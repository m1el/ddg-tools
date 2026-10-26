//! The components of §14.6, each timed on a pool of a given size, with a hash of its result to
//! check bitwise equality across thread counts.
#![allow(dead_code)]

use std::sync::Arc;
use std::time::Instant;

use ddg_repulsive_curves::metric::sobolev_gram;
use ddg_repulsive_curves::multigrid::HierarchyInput;
use ddg_repulsive_curves::solve::SaddleSolver;
use ddg_repulsive_curves::{
    BctOptions, BhOptions, CcdOptions, ConstraintSet, CurveNetwork, DenseSaddle, DofLayout,
    EnergyEval, Exponents, GramQuadrature, HierarchicalMetric, MgOptions, MultigridSaddle,
    default_constraints, energy, first_contact,
};

use super::P;

pub const COMPONENTS: &[&str] = &[
    "energy-exact",
    "energy-bh",
    "gram",
    "lu",
    "bct-build",
    "bct-product",
    "mg-build",
    "mg-solve",
    "ccd",
];

pub fn fnv(h: &mut u64, x: f64) {
    for b in x.to_bits().to_le_bytes() {
        *h ^= b as u64;
        *h = h.wrapping_mul(0x100000001b3);
    }
}

pub fn hash_slice(v: &[f64]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &x in v {
        fnv(&mut h, x);
    }
    h
}

fn flat(v: &[P]) -> Vec<f64> {
    v.iter().flat_map(|p| p.iter().copied()).collect()
}

pub fn pool(t: usize) -> Arc<rayon::ThreadPool> {
    Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build()
            .unwrap(),
    )
}

/// A fixed pseudo-random vector (splitmix64), the right-hand side of products and solves.
pub fn pseudo(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..n)
        .map(|_| {
            s = s.wrapping_add(0x9e3779b97f4a7c15);
            let mut z = s;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            z ^= z >> 31;
            (z >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        })
        .collect()
}

/// Inputs shared by the components (built once per curve, outside the timing).
pub struct Inputs<'a> {
    pub net: &'a CurveNetwork,
    pub pos: &'a [P],
    pub exps: Exponents,
    pub cset: ConstraintSet,
    pub jac: ddg_core::sparse::SparseMatrix,
    /// The exact (BH) differential, the right-hand side of the solves.
    pub grad: Vec<f64>,
    pub velocity: Vec<P>,
}

impl<'a> Inputs<'a> {
    pub fn new(net: &'a CurveNetwork, pos: &'a [P], exps: Exponents) -> Inputs<'a> {
        let cset = ConstraintSet::new(net, pos, &default_constraints(), None).unwrap();
        let jac = cset.jacobian(net, pos).unwrap();
        let r = energy(
            net,
            pos,
            exps,
            &EnergyEval::BarnesHut(BhOptions::default()),
            Some(8),
        )
        .unwrap();
        let grad = flat(&r.differential);
        // a typical step: 0.2 mean edge length along the descent direction
        let g = net.edge_geometry(pos);
        let mean_len = g.length.iter().sum::<f64>() / g.length.len() as f64;
        let mx = grad.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
        let velocity = r
            .differential
            .iter()
            .map(|d| d.map(|c| -c / mx * 0.2 * mean_len))
            .collect();
        Inputs {
            net,
            pos,
            exps,
            cset,
            jac,
            grad,
            velocity,
        }
    }
}

/// Dense data prepared outside the timing for `lu`.
pub struct DenseData {
    pub a: ddg_repulsive_curves::DenseMatrix,
}

pub fn dense_data(i: &Inputs) -> DenseData {
    let p = pool(8);
    DenseData {
        a: sobolev_gram(&p, i.net, i.pos, i.exps, GramQuadrature::Midpoint).a,
    }
}

/// Runs the component once on `t` threads; returns (seconds, result hash). The time is that of
/// the operation named, setup excluded (`bct-product` and `mg-solve` time the product / solve on
/// a structure built outside).
pub enum Prepared {
    None,
    Bct(Box<HierarchicalMetric>, Vec<f64>),
    Mg(Box<MultigridSaddle>),
}

pub fn prepare(c: &str, i: &Inputs, t: usize, dense: Option<&DenseData>) -> Prepared {
    let _ = dense;
    match c {
        "bct-product" => {
            let m = HierarchicalMetric::build(
                pool(t),
                i.net,
                i.pos,
                i.exps,
                GramQuadrature::Midpoint,
                DofLayout::new(None),
                BctOptions::default(),
            )
            .unwrap();
            let x = pseudo(m.n_dofs(), 7);
            Prepared::Bct(Box::new(m), x)
        }
        "mg-solve" => {
            let p = pool(t);
            let dp = pool(t.min(16));
            let input = HierarchyInput {
                network: i.net,
                positions: i.pos,
                exponents: i.exps,
                quadrature: GramQuadrature::Midpoint,
                layout: DofLayout::new(None),
                jacobian: &i.jac,
                constraints: Some(&i.cset),
            };
            Prepared::Mg(Box::new(
                MultigridSaddle::new(&p, &dp, &input, MgOptions::default()).unwrap(),
            ))
        }
        _ => Prepared::None,
    }
}

pub fn run(
    c: &str,
    i: &Inputs,
    t: usize,
    dense: Option<&DenseData>,
    prep: &mut Prepared,
) -> Result<(f64, u64), String> {
    let e = |e: ddg_repulsive_curves::RcError| e.to_string();
    match c {
        "energy-exact" | "energy-bh" => {
            let ev = if c == "energy-exact" {
                EnergyEval::Exact
            } else {
                EnergyEval::BarnesHut(BhOptions::default())
            };
            let t0 = Instant::now();
            let r = energy(i.net, i.pos, i.exps, &ev, Some(t)).map_err(e)?;
            let s = t0.elapsed().as_secs_f64();
            let mut h = hash_slice(&flat(&r.differential));
            fnv(&mut h, r.energy);
            Ok((s, h))
        }
        "gram" => {
            let p = pool(t);
            let t0 = Instant::now();
            let g = sobolev_gram(&p, i.net, i.pos, i.exps, GramQuadrature::Midpoint);
            let s = t0.elapsed().as_secs_f64();
            Ok((s, hash_slice(g.a.as_slice())))
        }
        "lu" => {
            let d = dense.ok_or("lu needs dense data")?;
            let p = pool(t);
            let t0 = Instant::now();
            let mut sad = DenseSaddle::new(p, &d.a, &i.jac, DofLayout::new(None)).map_err(e)?;
            let s = t0.elapsed().as_secs_f64();
            let k = i.jac.nrows();
            let sol = sad.solve(&i.grad, &vec![0.0; k]).map_err(e)?;
            Ok((s, hash_slice(&sol.x)))
        }
        "bct-build" => {
            let p = pool(t);
            let t0 = Instant::now();
            let m = HierarchicalMetric::build(
                p,
                i.net,
                i.pos,
                i.exps,
                GramQuadrature::Midpoint,
                DofLayout::new(None),
                BctOptions::default(),
            )
            .map_err(e)?;
            let s = t0.elapsed().as_secs_f64();
            let y = m.apply(&pseudo(m.n_dofs(), 7)).map_err(e)?;
            Ok((s, hash_slice(&y)))
        }
        "bct-product" => {
            let Prepared::Bct(m, x) = prep else {
                return Err("not prepared".into());
            };
            let t0 = Instant::now();
            let y = m.apply(x).map_err(e)?;
            Ok((t0.elapsed().as_secs_f64(), hash_slice(&y)))
        }
        "mg-build" => {
            let p = pool(t);
            let dp = pool(t.min(16));
            let input = HierarchyInput {
                network: i.net,
                positions: i.pos,
                exponents: i.exps,
                quadrature: GramQuadrature::Midpoint,
                layout: DofLayout::new(None),
                jacobian: &i.jac,
                constraints: Some(&i.cset),
            };
            let t0 = Instant::now();
            let mut m = MultigridSaddle::new(&p, &dp, &input, MgOptions::default()).map_err(e)?;
            let s = t0.elapsed().as_secs_f64();
            let sol = m.solve(&i.grad, &vec![0.0; i.jac.nrows()]).map_err(e)?;
            Ok((s, hash_slice(&sol.x)))
        }
        "mg-solve" => {
            let Prepared::Mg(m) = prep else {
                return Err("not prepared".into());
            };
            let t0 = Instant::now();
            let sol = m.solve(&i.grad, &vec![0.0; i.jac.nrows()]).map_err(e)?;
            let s = t0.elapsed().as_secs_f64();
            Ok((s, hash_slice(&sol.x)))
        }
        "ccd" => {
            let p = pool(t);
            let t0 = Instant::now();
            let r = first_contact(
                &p,
                i.net,
                i.pos,
                &i.velocity,
                1.0,
                &[],
                &CcdOptions::default(),
            )
            .map_err(e)?;
            let s = t0.elapsed().as_secs_f64();
            let mut h = hash_slice(&[r.contact.as_ref().map_or(-1.0, |c| c.tau)]);
            h ^= r.candidates.wrapping_mul(31) ^ r.advances;
            Ok((s, h))
        }
        _ => Err(format!("unknown component {c}")),
    }
}
