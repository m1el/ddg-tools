//! `spatial::block_tree::BlockTree` (WP17) against the curves crate's `BlockClusterTree`: the
//! same blocks and statistics, and **bitwise** equal products for any channel map, thread count,
//! leaf size and separation, with several right-hand sides through one pass.

use ddg_core::generators::{SplitMix64, icosphere, perturb, torus};
use ddg_repulsive_curves::spatial::bct::{BlockClusterTree, ClusterKernel};
use ddg_repulsive_curves::spatial::{Aabb, Bvh, Element};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::metric::hierarchical::separated;
use ddg_repulsive_surfaces::spatial::block_tree::BlockTree;
use ddg_repulsive_surfaces::spatial::farthest_corner_distance;

struct Pt {
    p: [f64; 3],
    m: f64,
}

impl Element<3> for Pt {
    fn point(&self) -> [f64; 3] {
        self.p
    }
    fn mass(&self) -> f64 {
        self.m
    }
}

struct Kern<'a> {
    bvh: &'a Bvh<3>,
    radius: Vec<f64>,
    boxes: Vec<Aabb<3>>,
    pts: &'a [[f64; 3]],
    chi: f64,
}

fn d2(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|k| (a[k] - b[k]) * (a[k] - b[k])).sum()
}

impl ClusterKernel<2> for Kern<'_> {
    fn admissible(&self, a: usize, b: usize) -> bool {
        separated(self.radius[a], &self.boxes[a], self.radius[b], &self.boxes[b], self.chi)
    }
    fn far(&self, a: usize, b: usize) -> [f64; 2] {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        let n = self.bvh.nodes();
        let r = d2(n[a].centroid(), n[b].centroid());
        [r.powf(-1.7), r.powf(-0.6) * 0.5]
    }
    fn near(&self, i: usize, j: usize) -> Option<[f64; 2]> {
        let (i, j) = if i <= j { (i, j) } else { (j, i) };
        let r = d2(self.pts[i], self.pts[j]);
        Some([r.powf(-1.7), r.powf(-0.6) * 0.5])
    }
    fn weight(&self, i: usize) -> f64 {
        1.0 + 0.1 * (i % 7) as f64
    }
}

fn points(which: usize) -> Vec<[f64; 3]> {
    let mut p = match which {
        0 => icosphere(4).unwrap().positions,
        _ => torus(60, 24, 1.5, 0.6).unwrap().positions,
    };
    perturb(&mut p, [0.01; 3], 3);
    p
}

fn bits<const C: usize>(v: &[[f64; C]]) -> Vec<u64> {
    v.iter().flatten().map(|x| x.to_bits()).collect()
}

#[test]
fn products_are_bitwise_those_of_the_curves_crate() {
    for which in 0..2 {
        let pts = points(which);
        let elems: Vec<Pt> = pts.iter().enumerate().map(|(i, &p)| Pt { p, m: 1.0 + (i % 5) as f64 }).collect();
        for leaf in [1usize, 3, 8] {
            let bvh = Bvh::build(&elems, leaf).unwrap();
            let boxes: Vec<Aabb<3>> = bvh.nodes().iter().map(|n| n.bbox).collect();
            let radius: Vec<f64> = bvh
                .nodes()
                .iter()
                .map(|n| farthest_corner_distance(&n.bbox, n.centroid()))
                .collect();
            for chi in [0.0, 0.25, 0.5, 1.0] {
                let k = Kern { bvh: &bvh, radius: radius.clone(), boxes: boxes.clone(), pts: &pts, chi };
                let mut rng = SplitMix64::new(11);
                let psi: Vec<[f64; 6]> =
                    (0..pts.len()).map(|_| std::array::from_fn(|_| rng.next_f64() - 0.5)).collect();
                let map = [0usize, 1, 0, 0, 1, 1];
                let p1 = make_pool(Some(1)).unwrap();
                let reference = BlockClusterTree::<2>::build(&p1, &bvh, &k, 64).unwrap();
                let want = reference.multiply(&p1, &map, &psi).unwrap();
                for threads in [1usize, 3, 8] {
                    let pool = make_pool(Some(threads)).unwrap();
                    let ours = BlockTree::<2>::build(&pool, &bvh, &k, 64).unwrap();
                    assert_eq!(ours.stats(), reference.stats(), "stats, leaf {leaf}, chi {chi}");
                    assert_eq!(bits(ours.row_sums()), bits(reference.row_sums()));
                    let got = ours.multiply(&pool, &map, &psi).unwrap();
                    assert_eq!(bits(&got), bits(&want), "leaf {leaf}, chi {chi}, {threads} threads");
                    // One channel at a time gives the same bits as the batch.
                    let one: Vec<[f64; 1]> = psi.iter().map(|x| [x[4]]).collect();
                    let g1 = ours.multiply(&pool, &[1], &one).unwrap();
                    assert!(g1.iter().zip(&got).all(|(a, b)| a[0].to_bits() == b[4].to_bits()));
                }
            }
        }
    }
}
