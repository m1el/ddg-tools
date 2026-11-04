//! Determinism (SPEC.md §5.4, §15.1-7): the whole pipeline is bitwise identical run to run and
//! across 1, 3, 8 and 32 threads (kNN, coefficients, σ, tori, SDF, gradient, Laplacian).

mod common;

use common::{bits, ellipsoid_cloud, lcg};
use ddg_points_as_tori::{EvalOut, FeatureOrder, Model, Neighborhood, PointsAsTori, TorusField};

fn run(cloud: &ddg_points_as_tori::PointCloud, model: &Model, threads: usize, order: FeatureOrder, nbh: Neighborhood) -> Vec<u64> {
    let mut opts = common::ref_build().with_threads(Some(threads));
    opts.predict.feature_order = order;
    opts.eval.neighborhood = nbh;
    let (field, rep) = PointsAsTori::build(cloud, model, &opts).unwrap();
    assert_eq!(rep.n, cloud.len());
    let mut out = Vec::new();
    for c in &field.tori().coefficients {
        out.extend(c.iter().map(|x| x.to_bits() as u64));
    }
    out.extend(bits(&field.tori().sigma));
    for t in &field.tori().tori {
        out.extend(bits(&t.center));
        out.extend(bits(&t.axis));
        out.extend(bits(&[t.major, t.minor]));
    }
    out.extend(eval_bits(&field, 3000));
    out
}

fn eval_bits(field: &TorusField, nq: usize) -> Vec<u64> {
    let mut s = 5;
    let qs: Vec<[f64; 3]> = (0..nq).map(|_| [1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8, 1.6 * lcg(&mut s) - 0.8]).collect();
    let mut d = vec![0.0; nq];
    let mut g = vec![[0.0; 3]; nq];
    let mut l = vec![0.0; nq];
    field.evaluate(&qs, EvalOut { distance: Some(&mut d), gradient: None, laplacian: None }).unwrap();
    let mut out = bits(&d);
    field.evaluate(&qs, EvalOut { distance: Some(&mut d), gradient: Some(&mut g), laplacian: Some(&mut l) }).unwrap();
    out.extend(bits(&d));
    out.extend(g.iter().flat_map(|x| bits(x)));
    out.extend(bits(&l));
    out
}

#[test]
fn bitwise_independent_of_threads() {
    let model = Model::reference().unwrap();
    let cloud = ellipsoid_cloud(if cfg!(debug_assertions) { 300 } else { 2500 }, 11, [0.6, 0.45, 0.3]);
    for (order, nbh) in [
        (FeatureOrder::DifferenceFirst, Neighborhood::Knn { k: 32 }),
        (FeatureOrder::Reference, Neighborhood::Radius { r: 0.1, fallback_k: 32 }),
    ] {
        let base = run(&cloud, &model, 1, order, nbh);
        let again = run(&cloud, &model, 1, order, nbh);
        assert!(base == again, "{order:?}: two runs on 1 thread differ");
        for t in [3, 8, 32] {
            let r = run(&cloud, &model, t, order, nbh);
            let first = base.iter().zip(&r).position(|(a, b)| a != b);
            assert!(first.is_none(), "{order:?}: {t} threads differ from 1 thread at word {first:?}");
        }
    }
}
