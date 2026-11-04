//! A Fig. 10-like picture of the 2D variant: input points with normals, the fitted circles, the
//! blended SDF with level lines, and the zero contour.
//!
//! `cargo run --release -p ddg-points-as-tori --example pat2d -- [out_dir] [points] [shape ...]`
//!
//! Shapes: circle ellipse box hexagon heart pill moon pie (default: heart box moon pie).
//! Writes `<out_dir>/pat2d-<shape>.png` and `.svg` (default `runs_tmp/pat-wp9a/`).

use ddg_points_as_tori::analytic::{Kind, Shape2d};
use ddg_points_as_tori::generate::d2::Instance2d;
use ddg_points_as_tori::generate::rng::Rng;
use ddg_points_as_tori::pat2d::*;
use ddg_points_as_tori::NormalPolicy;
use std::path::PathBuf;

fn shape(name: &str) -> Option<Shape2d> {
    let (k, p): (Kind, &[f64]) = match name {
        "circle" => (Kind::Circle, &[0.5]),
        "ellipse" => (Kind::Ellipse, &[0.6, 0.35]),
        "box" => (Kind::Box, &[0.5, 0.3]),
        "hexagon" => (Kind::Hexagon, &[0.5]),
        "heart" => (Kind::Heart, &[]),
        "pill" => (Kind::Pill, &[0.5, 0.2]),
        "moon" => (Kind::Moon, &[0.5, 0.6, 0.35]),
        "pie" => (Kind::Pie, &[0.6, 1.0]),
        _ => return None,
    };
    Shape2d::new(k, p).ok()
}

/// Blue inside, orange outside, light near zero.
fn color(v: f64, scale: f64) -> [u8; 3] {
    let t = (v / scale).clamp(-1.0, 1.0);
    let mix = |a: [f64; 3], b: [f64; 3], s: f64| [0, 1, 2].map(|i| (a[i] + (b[i] - a[i]) * s).round() as u8);
    let white = [250.0, 250.0, 250.0];
    if t < 0.0 { mix(white, [70.0, 120.0, 190.0], -t) } else { mix(white, [235.0, 150.0, 70.0], t) }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let out = PathBuf::from(args.first().cloned().unwrap_or_else(|| "runs_tmp/pat-wp9a".into()));
    let n: usize = args.get(1).map(|s| s.parse()).transpose()?.unwrap_or(70);
    let names: Vec<String> = if args.len() > 2 { args[2..].to_vec() } else { ["heart", "box", "moon", "pie"].map(String::from).to_vec() };
    for name in names {
        let Some(s) = shape(&name) else {
            eprintln!("unknown shape {name}");
            continue;
        };
        let c = Instance2d::Analytic(s.clone()).sample_uniform(n, &mut Rng::new(4))?;
        let (cloud, _) = Cloud2::new(c.positions, c.normals, NormalPolicy::Normalize)?;
        let (field, rep) = Pat2d::build(&cloud, &ClassicalFit::default(), &BuildOptions2::default())?;
        let bb = s.bbox();
        let w = (bb[1][0] - bb[0][0]).max(bb[1][1] - bb[0][1]);
        let (cx, cy) = (0.5 * (bb[0][0] + bb[1][0]), 0.5 * (bb[0][1] + bb[1][1]));
        let half = 0.5 * w * 1.35;
        let bbox = [[cx - half, cy - half], [cx + half, cy + half]];
        let px = 900;
        let g = evaluate_grid(&field, bbox, px + 1, px + 1)?;
        // exact comparison
        let nodes = g.nodes();
        let err: f64 = nodes.iter().zip(&g.values).map(|(p, v)| (v - s.distance(*p)).abs()).sum::<f64>() / nodes.len() as f64;
        let mut r = Raster::new(px, px, bbox, [255, 255, 255])?;
        let scale = 0.5 * w;
        let step = 0.04 * w;
        r.fill_with(|p| {
            let v = field.distance(p).unwrap_or(0.0);
            let mut c = color(v, scale);
            // level lines every `step`, thin
            let f = (v / step).rem_euclid(1.0);
            let d = f.min(1.0 - f) * step;
            if d < 0.0012 * w {
                for ch in &mut c {
                    *ch = (*ch as f64 * 0.7) as u8;
                }
            }
            c
        });
        // fitted circles: an arc of the extent of the neighborhood around each contact point
        let circles = field.circles();
        for (i, prim) in circles.primitives.iter().enumerate() {
            if let CirclePrimitive::Circle { center, radius, .. } = *prim
                && radius < 3.0 * w
            {
                let p = circles.contact[i];
                let dir = [p[0] - center[0], p[1] - center[1]];
                let arc = (0.12 * w / radius).min(1.2);
                r.circle_arc(center, radius, dir, arc, 1.2, [40, 140, 70]);
            }
        }
        for poly in chain_segments(&marching_squares(&g, 0.0)) {
            r.polyline(&poly, 2.5, [0, 0, 0]);
        }
        for (p, nrm) in cloud.positions().iter().zip(cloud.normals()) {
            r.line(*p, [p[0] + 0.05 * w * nrm[0], p[1] + 0.05 * w * nrm[1]], 1.2, [200, 30, 60]);
            r.dot(*p, 3.0, [200, 30, 60]);
        }
        write_png(&out.join(format!("pat2d-{name}.png")), &r)?;
        let polys = chain_segments(&marching_squares(&g, 0.0));
        write_svg_contours(&out.join(format!("pat2d-{name}.svg")), bbox, 700, &polys, cloud.positions())?;
        println!("{name}: n = {n}, circles {} lines {}, mean |phi - exact| on the picture {err:.5} (width {w:.2}), contours {}", rep.fit.circles, rep.fit.lines, polys.len());
    }
    Ok(())
}
