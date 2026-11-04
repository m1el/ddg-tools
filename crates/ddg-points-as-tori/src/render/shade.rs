//! Shading and colour maps, ported from the authors' demo shader (`demo/shaders_3d.py`, MIT,
//! Copyright (c) 2026 Nicole Feng): `diffuseShade`, `fresnelShade`, `normalShade` with its GGX
//! term (from shadertoy `MlB3DV`, `7tKfz1`), `shade`, `colormap_SDF` and `darken_color`.
//!
//! GLSL `pow` of a negative base is undefined; here the base is clamped at 0 where the demo's
//! arithmetic can make it negative (documented in `NOTES-WP11.md`). Everything is f64 and scalar
//! per pixel, so the result does not depend on the thread count.
//!
//! `[formula: pat-render-ggx]` `D = α²/(π((n·h)²(α²−1)+1)²)`, `F = F0 + (1−F0)(1−l·h)⁵`,
//! `V = 1/((n·l)(1−k)+k) · 1/((n·v)(1−k)+k)`, `k = α/2`, `α = roughness²`; the term is `D F V`.
//!
//! `[formula: pat-render-fresnel]` `(1 − |ω·n|)^80` (the rim term) and `(1 − n·v)₊^p`,
//! `p = mix(5, 3.5, metallic)` (the Schlick-type weight in `normalShade`).

use super::cubemap::CubeMap;

pub(crate) type V3 = [f64; 3];

pub(crate) fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(crate) fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(crate) fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(crate) fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(crate) fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(crate) fn length(a: V3) -> f64 {
    dot(a, a).sqrt()
}
pub(crate) fn normalize(a: V3) -> V3 {
    let l = length(a);
    scale(a, 1.0 / l)
}

fn mix1(x: f64, y: f64, a: f64) -> f64 {
    x * (1.0 - a) + y * a
}
fn mix3(x: V3, y: V3, a: f64) -> V3 {
    [
        mix1(x[0], y[0], a),
        mix1(x[1], y[1], a),
        mix1(x[2], y[2], a),
    ]
}
fn pow3(x: V3, p: f64) -> V3 {
    // GLSL pow(x, y) is undefined for x < 0; clamp at 0
    [
        x[0].max(0.0).powf(p),
        x[1].max(0.0).powf(p),
        x[2].max(0.0).powf(p),
    ]
}
fn reflect(i: V3, n: V3) -> V3 {
    sub(i, scale(n, 2.0 * dot(n, i)))
}

/// The lights of the demo, in world space.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Lights {
    /// The point light of the diffuse term and of the shadow rays.
    pub position: V3,
    /// The direction (not normalized) of the specular light.
    pub specular_dir: V3,
}

fn ggx(n: V3, v: V3, l: V3, roughness: f64, f0: f64) -> f64 {
    let alpha = roughness * roughness;
    let h = normalize(add(v, l));
    let clamp01 = |x: f64| x.clamp(0.0, 1.0);
    let dot_nl = clamp01(dot(n, l));
    let dot_nv = clamp01(dot(n, v));
    let dot_nh = clamp01(dot(n, h));
    let dot_lh = clamp01(dot(l, h));
    let alpha_sqr = alpha * alpha;
    let pi = std::f64::consts::PI;
    let denom = dot_nh * dot_nh * (alpha_sqr - 1.0) + 1.0;
    let d = alpha_sqr / (pi * denom * denom);
    let dot_lh5 = (1.0 - dot_lh).powf(5.0);
    let f = f0 + (1.0 - f0) * dot_lh5;
    let k = alpha / 2.0;
    let g1v = |x: f64| 1.0 / (x * (1.0 - k) + k);
    let vis = g1v(dot_nl) * g1v(dot_nv);
    d * f * vis
}

/// `normalShade`: IBL diffuse + GGX specular + reflection, gamma-corrected by 1/2.2.
fn normal_shade(env: &CubeMap, ray: V3, normal: V3, spec_pos: V3, material: V3) -> V3 {
    let metallic = 0.07;
    let roughness = 0.1;
    let fresnel_pow = mix1(5.0, 3.5, metallic);
    let sample = env.texture([1.0, 0.0, 0.0]);
    let light_color = pow3(scale(sample, 1.2), 2.2);

    let r = reflect(ray, normal);
    let ibl_diffuse = pow3(env.blurred(normal), 2.2);
    let ibl_reflection = pow3(env.blurred(r), 2.2);

    let fresnel = (1.0 - dot(normal, scale(ray, -1.0)))
        .max(0.0)
        .powf(fresnel_pow);

    let mut refl = pow3(env.texture(r), 2.2);
    refl = mix3(refl, ibl_reflection, (1.0 - fresnel) * roughness);
    refl = mix3(refl, ibl_reflection, roughness);

    let light = normalize(spec_pos);
    let g = ggx(normal, scale(ray, -1.0), light, roughness * 0.7, 0.2);
    let spec = scale(light_color, g);
    refl = sub(refl, spec);

    let diff0 = [
        ibl_diffuse[0] * material[0].max(0.0).powf(2.2),
        ibl_diffuse[1] * material[1].max(0.0).powf(2.2),
        ibl_diffuse[2] * material[2].max(0.0).powf(2.2),
    ];
    let diff = mix3(diff0, refl, fresnel);
    let m = mix3(diff, refl, metallic);
    let color = [
        (m[0] + spec[0]).min(1.0),
        (m[1] + spec[1]).min(1.0),
        (m[2] + spec[2]).min(1.0),
    ];
    pow3(color, 1.0 / 2.2)
}

/// `shade` followed by the `sqrt` the demo applies to surface colours. `normal` must face the
/// viewer.
pub(crate) fn shade_surface(
    env: &CubeMap,
    lights: &Lights,
    pos: V3,
    ray: V3,
    normal: V3,
    material: V3,
) -> V3 {
    let ndl = dot(normalize(sub(lights.position, pos)), normal).max(0.3);
    let diffuse = scale(material, ndl);
    let ns = normal_shade(env, ray, normal, lights.specular_dir, material);
    let fr = (1.0 - dot(ray, normal).abs()).powf(80.0);
    let mut c = [0.0; 3];
    for a in 0..3 {
        c[a] = (diffuse[a] + 0.3 * ns[a] + 0.6 * fr).clamp(0.0, 1.0).sqrt();
    }
    c
}

fn smoothstep(e0: f64, e1: f64, x: f64) -> f64 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// `colormap_SDF` on `t = (φ − min)/(max − min)`: white to Pasadena blue (`pow(2t, 7)`) on
/// `[0, ½)`, white to grapefruit on `[½, 1)`, black outside `[0, 1)` (the demo's `segment`).
pub(crate) fn colormap_sdf(t: f64) -> V3 {
    let white = [1.0, 1.0, 1.0];
    if (0.0..0.5).contains(&t) {
        let blue = [0.0627451, 0.517647, 0.94902];
        let tp = (2.0 * t).powf(7.0);
        [
            (1.0 - tp) * white[0] + tp * blue[0],
            (1.0 - tp) * white[1] + tp * blue[1],
            (1.0 - tp) * white[2] + tp * blue[2],
        ]
    } else if (0.5..1.0).contains(&t) {
        let grapefruit = [0.992157, 0.431373, 0.337255];
        mix3(white, grapefruit, 2.0 * (t - 0.5))
    } else {
        [0.0; 3]
    }
}

/// `redistribute_rgb`: keep the hue when a channel exceeds 1.
fn redistribute_rgb(c: V3) -> V3 {
    let m = c[0].max(c[1]).max(c[2]);
    if m <= 1.0 {
        return c;
    }
    let total = c[0] + c[1] + c[2];
    if total >= 3.0 {
        return [1.0; 3];
    }
    let x = (3.0 - total) / (3.0 * m - total);
    let gray = 1.0 - x * m;
    [gray + x * c[0], gray + x * c[1], gray + x * c[2]]
}

/// The colour of the cut plane at signed distance `phi` (already minus the level), as the demo's
/// `cutplaneColor`: colour map over `[min_d, max_d]` darkened by `0.9 (1 − 0.1 smoothstep(sin(−f φ)))`.
pub(crate) fn cutplane_color(phi: f64, min_d: f64, max_d: f64, isoband_frequency: f64) -> V3 {
    let col = colormap_sdf((phi - min_d) / (max_d - min_d));
    let dark = 0.9 * (1.0 - 0.1 * smoothstep(-0.01, 0.01, (-isoband_frequency * phi).sin()));
    redistribute_rgb(scale(col, dark))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colormap_endpoints() {
        assert_eq!(colormap_sdf(0.0), [1.0, 1.0, 1.0]);
        let near_zero_inside = colormap_sdf(0.4999999);
        assert!((near_zero_inside[0] - 0.0627451).abs() < 1e-3);
        assert_eq!(colormap_sdf(0.5), [1.0, 1.0, 1.0]);
        assert_eq!(colormap_sdf(1.5), [0.0; 3]);
    }

    #[test]
    fn shading_is_finite_and_bounded() {
        let env = CubeMap::forest().unwrap();
        let lights = Lights {
            position: [0.0, 10.0, 2.0],
            specular_dir: [-0.5, 1.0, 0.0],
        };
        for n in [
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            normalize([1.0, 1.0, -1.0]),
        ] {
            let ray = normalize([0.3, -0.2, -1.0]);
            let n = if dot(n, ray) > 0.0 { scale(n, -1.0) } else { n };
            let c = shade_surface(&env, &lights, [0.0; 3], ray, n, [0.27, 0.29, 0.57]);
            assert!(
                c.iter().all(|x| x.is_finite() && (0.0..=1.0).contains(x)),
                "{c:?}"
            );
        }
    }
}
