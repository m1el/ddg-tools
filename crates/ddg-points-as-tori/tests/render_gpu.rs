//! The CUDA sphere tracer (WP11). All tests are `#[ignore]`d: they open the GPU, which is busy with
//! the training run until WP16. Run after training: `cargo test --release --features cuda --test
//! render_gpu -- --ignored --nocapture`.

#![cfg(feature = "cuda")]

mod common;

use common::ellipsoid_cloud;
use ddg_points_as_tori::render::{Camera, Tracer, TraceOptions};
use ddg_points_as_tori::{BuildOptions, Device, Model, PointsAsTori};

fn field(device: Device) -> ddg_points_as_tori::TorusField {
    let cloud = ellipsoid_cloud(400, 11, [1.0, 0.7, 0.5]);
    let mut bo = BuildOptions::default().with_threads(Some(4));
    bo.eval.device = device;
    PointsAsTori::build(&cloud, &Model::reference().unwrap(), &bo).unwrap().0
}

/// GPU image = CPU image (same traversal arithmetic; the blend's `exp` is the crate's own, so the
/// hit positions agree to rounding), and the GPU image is bitwise reproducible.
#[test]
#[ignore = "run after training (WP16)"]
fn gpu_trace_matches_cpu_and_repeats() {
    let cpu_field = field(Device::Cpu);
    let gpu_field = field(Device::Cuda { ordinal: 0 });
    let cam = Camera::look_at([0.0, 0.8, 4.0], [0.0; 3], [0.0, 1.0, 0.0], 30.0, 160, 120).unwrap();
    let cpu = Tracer::new(TraceOptions { threads: Some(8), ..Default::default() }).unwrap();
    let gpu = Tracer::new(TraceOptions { device: Device::Cuda { ordinal: 0 }, ..Default::default() }).unwrap();
    let (ic, rc) = cpu.render_torus(&cpu_field, &cam).unwrap();
    let (ig, rg) = gpu.render_torus(&gpu_field, &cam).unwrap();
    let (ig2, rg2) = gpu.render_torus(&gpu_field, &cam).unwrap();
    assert_eq!(ig, ig2, "GPU render differs run to run");
    assert_eq!((rg.hits, rg.misses, rg.steps_total), (rg2.hits, rg2.misses, rg2.steps_total));
    let diff = ic.rgb.iter().zip(&ig.rgb).filter(|(a, b)| a != b).count();
    eprintln!("CPU hits {} GPU hits {}; differing bytes {} of {}", rc.hits, rg.hits, diff, ic.rgb.len());
    assert!((rc.hits as i64 - rg.hits as i64).abs() <= 2, "{rc:?} {rg:?}");
    assert!(diff * 200 < ic.rgb.len(), "more than 0.5% of the bytes differ");
}

/// A CPU field cannot be traced on the GPU (an error before any device is opened).
#[test]
fn gpu_tracer_rejects_cpu_field() {
    let f = field(Device::Cpu);
    let cam = Camera::look_at([0.0, 0.8, 4.0], [0.0; 3], [0.0, 1.0, 0.0], 30.0, 16, 16).unwrap();
    let gpu = Tracer::new(TraceOptions { device: Device::Cuda { ordinal: 0 }, ..Default::default() }).unwrap();
    assert!(gpu.render_torus(&f, &cam).is_err());
}
