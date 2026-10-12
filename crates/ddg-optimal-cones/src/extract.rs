//! Extracting cones from the optimal measure (paper Sec. 6.1): thresholding, and the optional
//! merging of edge-connected clusters of cones of one sign into single cones.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use ddg_core::io::Cone;
use ddg_core::{Geometry, TriMesh};

use crate::problem::VertexData;

/// A cluster of cones that was replaced by one cone.
#[derive(Debug, Clone, PartialEq)]
pub struct MergedCluster {
    /// The cones of the optimal measure that formed the cluster, by ascending vertex. They all
    /// have the same sign.
    pub members: Vec<Cone>,
    /// The cone that replaced them: the sum of the members' angles at the cluster's Fréchet
    /// mean vertex.
    pub cone: Cone,
    /// The diameter of the cluster: the largest distance between two of its members, measured
    /// along mesh edges inside the cluster, in the length unit of the geometry.
    pub diameter: f64,
}

/// Why a cluster of cones was left as it is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnmergedReason {
    /// The cluster has more cones than
    /// [`max_merge_cluster_size`](crate::SolveOptions::max_merge_cluster_size).
    TooLarge,
    /// The cluster is wider than
    /// [`max_merge_cluster_diameter`](crate::SolveOptions::max_merge_cluster_diameter).
    TooWide,
    /// The summed angle violates the angle bounds of every vertex of the cluster.
    BoundsViolated,
}

/// A cluster of cones that was not merged.
#[derive(Debug, Clone, PartialEq)]
pub struct UnmergedCluster {
    /// The vertices of the cluster, ascending.
    pub vertices: Vec<u32>,
    /// The sum of its cone angles.
    pub angle: f64,
    /// Its diameter (see [`MergedCluster::diameter`]). Not computed for a cluster that was
    /// rejected for its number of cones: then it is NaN.
    pub diameter: f64,
    /// Why it was left alone.
    pub reason: UnmergedReason,
}

/// What cone extraction did to the optimal measure. Everything that makes the returned cones
/// differ from the raw measure is listed here.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtractionReport {
    /// The threshold that was applied: a vertex is a cone if `|µ_i|` exceeds it.
    pub threshold: f64,
    /// Number of vertices with a nonzero `|µ_i|` at or below the threshold; their angles were
    /// dropped.
    pub n_dropped: usize,
    /// The sum of the dropped angles.
    pub dropped_angle: f64,
    /// Whether merging of cone clusters was requested.
    pub merge_requested: bool,
    /// Number of clusters found: sets of two or more cones **of the same sign** that are
    /// connected through mesh edges between their members. Counted whether or not merging was
    /// requested. A positive and a negative cone next to each other are never one cluster.
    pub n_clusters: usize,
    /// The clusters that were merged.
    pub merged: Vec<MergedCluster>,
    /// The clusters that were found but not merged (only filled when merging was requested).
    pub unmerged: Vec<UnmergedCluster>,
    /// Whether the returned cones differ from the optimal measure (something was dropped or
    /// merged). If so the returned scale factor was recomputed for the returned cones.
    pub cones_differ_from_measure: bool,
    /// `distortion(returned cones) − distortion(optimal measure)`; zero if nothing changed.
    pub distortion_change: f64,
    /// `Σ |angle|` over the returned cones minus `Σ |µ_i|` of the optimal measure. Merging
    /// does not change the total (a cluster has one sign); dropping cones below the threshold
    /// lowers it.
    pub total_angle_change: f64,
    /// Number of returned cones whose angle is `≥ 2π` or `≤ −2π`. A cone angle of `2π` or more
    /// does not exist (the angle sum around the vertex would not be positive), and one of
    /// `−2π` or less has an angle sum of `4π` or more, which no flattening lays out without
    /// overlap. The linearized problem knows neither limit: such angles occur at a free cone
    /// site of a closed convex mesh, at needle-like spikes, and in merged clusters. Angle
    /// bounds inside `(−2π, 2π)` exclude them (a merged cone respects the bounds of its
    /// vertex).
    pub n_beyond_two_pi: usize,
}

#[derive(PartialEq)]
struct QueueEntry {
    distance: f64,
    node: usize,
}

impl Eq for QueueEntry {}

impl Ord for QueueEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        // Min-heap on the distance.
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.node.cmp(&self.node))
    }
}

impl PartialOrd for QueueEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// The Fréchet means of a cluster of cones, best first, and the diameter of the cluster.
///
/// The cluster is a set of vertices `c_1, …, c_m` that is connected through mesh edges between
/// its members, with cone angles `µ_1, …, µ_m`. Distances are measured inside the cluster:
/// `dist(v, c)` is the length of the shortest path from `v` to `c` along mesh edges whose two
/// ends both belong to the cluster. The Fréchet variance of a member `v` is
/// `Σ_k |µ_k| · dist(v, c_k)²`, and the Fréchet mean is the member that minimizes it (ties go to
/// the smaller vertex index). The diameter is the largest `dist` between two members.
///
/// Returns the members' positions in `members`, sorted by ascending variance, and the diameter.
///
/// [formula: cone-cluster-frechet-mean]
fn frechet_order(
    mesh: &TriMesh,
    geometry: &Geometry,
    members: &[u32],
    angles: &[f64],
) -> (Vec<usize>, f64) {
    let m = members.len();
    // Adjacency inside the cluster. `members` is ascending, so a binary search finds a vertex.
    let local = |v: usize| members.binary_search(&(v as u32)).ok();
    let mut neighbors: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
    for (a, &v) in members.iter().enumerate() {
        for w in mesh.vertex_neighbors(v as usize) {
            if let Some(b) = local(w) {
                let e = mesh
                    .find_edge(v as usize, w)
                    .expect("neighbors share an edge");
                neighbors[a].push((b, geometry.edge_lengths()[e]));
            }
        }
    }
    let mut variance = vec![0.0; m];
    let mut diameter = 0.0f64;
    let mut dist = vec![f64::INFINITY; m];
    let mut heap = BinaryHeap::new();
    for source in 0..m {
        dist.fill(f64::INFINITY);
        dist[source] = 0.0;
        heap.push(QueueEntry {
            distance: 0.0,
            node: source,
        });
        while let Some(QueueEntry { distance, node }) = heap.pop() {
            if distance > dist[node] {
                continue;
            }
            for &(next, length) in &neighbors[node] {
                let candidate = distance + length;
                if candidate < dist[next] {
                    dist[next] = candidate;
                    heap.push(QueueEntry {
                        distance: candidate,
                        node: next,
                    });
                }
            }
        }
        // dist[v] = dist(v, source) by symmetry.
        for v in 0..m {
            variance[v] += angles[source].abs() * dist[v] * dist[v];
            diameter = diameter.max(dist[v]);
        }
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| variance[a].total_cmp(&variance[b]).then(a.cmp(&b)));
    (order, diameter)
}

/// Options of [`extract_cones`].
pub(crate) struct ExtractOptions {
    pub threshold: f64,
    pub merge: bool,
    pub max_cluster_size: usize,
    pub max_cluster_diameter: f64,
}

/// Turns the optimal measure (per vertex) into a list of cones.
///
/// 1. A vertex is a cone if `|µ_i| > threshold` (the paper uses `10⁻¹²`).
/// 2. A *cluster* is a set of two or more cones **of the same sign** that is connected through
///    mesh edges between its members. Prescribed cones never belong to a cluster.
/// 3. If merging is requested, each cluster is replaced by one cone whose angle is the sum of
///    the cluster's angles, placed at the cluster's Fréchet mean vertex. The merged cone must
///    respect the angle bounds of the vertex it is placed at: if the mean vertex does not admit
///    the summed angle, the member with the next smallest Fréchet variance is taken; if no
///    member admits it, the cluster is left unmerged. A cluster with more than
///    `max_cluster_size` cones, or with a diameter above `max_cluster_diameter`, is left
///    unmerged.
///
/// Why one sign only. The paper merges because a single cone of the smooth optimum can show up
/// as a tiny ring of cones of that sign (Sec. 4.5, 6.1), and replacing the ring by its sum
/// moves little curvature over a short distance. A positive cone next to a negative one is
/// something else: a feature of the surface (the tip and the base of a spike, the two sides of
/// a crease). Summing them cancels cone angle that the optimum needs, and the distortion
/// returns: on the dataset's `camelshotnoise` that raised it from 3.5 to about 10⁴.
///
/// The sum of the returned angles equals the sum of the measure minus the reported
/// `dropped_angle` (up to rounding in the summation).
///
/// Returns the cones by ascending vertex, and the report with `distortion_change` and
/// `total_angle_change` still zero.
pub(crate) fn extract_cones(
    mesh: &TriMesh,
    geometry: &Geometry,
    data: &VertexData,
    prescribed: &[Cone],
    mu: &[f64],
    options: &ExtractOptions,
) -> (Vec<Cone>, ExtractionReport) {
    let nv = mesh.n_vertices();
    let mut report = ExtractionReport {
        threshold: options.threshold,
        n_dropped: 0,
        dropped_angle: 0.0,
        merge_requested: options.merge,
        n_clusters: 0,
        merged: Vec::new(),
        unmerged: Vec::new(),
        cones_differ_from_measure: false,
        distortion_change: 0.0,
        total_angle_change: 0.0,
        n_beyond_two_pi: 0,
    };
    let mut is_cone = vec![false; nv];
    for v in 0..nv {
        if mu[v].abs() > options.threshold {
            is_cone[v] = true;
        } else if mu[v] != 0.0 {
            report.n_dropped += 1;
            report.dropped_angle += mu[v];
        }
    }
    let mut is_prescribed = vec![false; nv];
    for cone in prescribed {
        is_prescribed[cone.vertex as usize] = true;
    }

    // Clusters: connected components of the cones (without the prescribed ones).
    let clusterable: Vec<bool> = (0..nv).map(|v| is_cone[v] && !is_prescribed[v]).collect();
    let mut angle = mu.to_vec();
    let mut visited = vec![false; nv];
    let mut stack = Vec::new();
    for seed in 0..nv {
        if visited[seed] || !clusterable[seed] {
            continue;
        }
        visited[seed] = true;
        stack.push(seed);
        let mut members: Vec<u32> = Vec::new();
        while let Some(v) = stack.pop() {
            members.push(v as u32);
            for w in mesh.vertex_neighbors(v) {
                if !visited[w] && clusterable[w] && (mu[w] > 0.0) == (mu[seed] > 0.0) {
                    visited[w] = true;
                    stack.push(w);
                }
            }
        }
        if members.len() < 2 {
            continue;
        }
        report.n_clusters += 1;
        if !options.merge {
            continue;
        }
        members.sort_unstable();
        let angles: Vec<f64> = members.iter().map(|&v| mu[v as usize]).collect();
        let total: f64 = angles.iter().sum();
        if members.len() > options.max_cluster_size {
            report.unmerged.push(UnmergedCluster {
                vertices: members,
                angle: total,
                diameter: f64::NAN,
                reason: UnmergedReason::TooLarge,
            });
            continue;
        }
        let (order, diameter) = frechet_order(mesh, geometry, &members, &angles);
        if diameter > options.max_cluster_diameter {
            report.unmerged.push(UnmergedCluster {
                vertices: members,
                angle: total,
                diameter,
                reason: UnmergedReason::TooWide,
            });
            continue;
        }
        let target = order
            .into_iter()
            .map(|k| members[k])
            .find(|&v| data.sites[v as usize].bound_violation(total) == 0.0);
        let Some(target) = target else {
            report.unmerged.push(UnmergedCluster {
                vertices: members,
                angle: total,
                diameter,
                reason: UnmergedReason::BoundsViolated,
            });
            continue;
        };
        // The members have one sign, so the sum is larger in magnitude than each of them: it
        // is a cone.
        for &v in &members {
            angle[v as usize] = 0.0;
            is_cone[v as usize] = false;
        }
        angle[target as usize] = total;
        is_cone[target as usize] = true;
        report.merged.push(MergedCluster {
            members: members
                .iter()
                .zip(&angles)
                .map(|(&vertex, &angle)| Cone { vertex, angle })
                .collect(),
            cone: Cone {
                vertex: target,
                angle: total,
            },
            diameter,
        });
    }
    report.cones_differ_from_measure = report.n_dropped > 0 || !report.merged.is_empty();
    let cones: Vec<Cone> = (0..nv)
        .filter(|&v| is_cone[v])
        .map(|v| Cone {
            vertex: v as u32,
            angle: angle[v],
        })
        .collect();
    report.n_beyond_two_pi = cones
        .iter()
        .filter(|c| c.angle.abs() >= 2.0 * std::f64::consts::PI)
        .count();
    (cones, report)
}
