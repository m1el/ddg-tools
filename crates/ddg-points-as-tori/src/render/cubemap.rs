//! The environment of the IBL term: the demo's blurred forest cube map (or a constant colour).
//!
//! Ported from the authors' demo (`demo/shaders_3d.py`: `textureAVG`, `textureBlurred`;
//! `demo/demo_3d.py`: the face order and `GL_LINEAR` / `GL_CLAMP_TO_EDGE` sampling). The six PNG
//! faces in `assets/` are the authors' (MIT, Copyright (c) 2026 Nicole Feng; notice in
//! `assets/LICENSE-points-as-tori-MIT.txt`) and are embedded in the binary.
//!
//! The demo's `textureBlurred(dir)` averages 9 texture taps in each of six fixed directions on
//! every call; those six colours do not depend on `dir`, so they are computed once here
//! (the same arithmetic, in the same order).

use crate::PatError;

const FACE_BYTES: [&[u8]; 6] = [
    include_bytes!("../../assets/ForestBlurred0.png"),
    include_bytes!("../../assets/ForestBlurred1.png"),
    include_bytes!("../../assets/ForestBlurred2.png"),
    include_bytes!("../../assets/ForestBlurred3.png"),
    include_bytes!("../../assets/ForestBlurred4.png"),
    include_bytes!("../../assets/ForestBlurred5.png"),
];

/// A cube map (six faces in the order +X, -X, +Y, -Y, +Z, -Z) or a constant colour, with the
/// demo's six pre-averaged directions.
#[derive(Debug, Clone, PartialEq)]
pub struct CubeMap {
    size: usize,
    /// Per face, `size * size` RGB texels in `[0, 1]` (row 0 first, as the PNG file).
    faces: Vec<Vec<[f64; 3]>>,
    constant: Option<[f64; 3]>,
    /// `textureAVG` at +X, +Y, +Z, -X, -Y, -Z (the demo's `r t f l b a`).
    avg: [[f64; 3]; 6],
}

impl CubeMap {
    /// The demo's `ForestBlurred0..5.png` (embedded).
    pub fn forest() -> Result<CubeMap, PatError> {
        let mut size = 0usize;
        let mut faces = Vec::with_capacity(6);
        for (i, bytes) in FACE_BYTES.iter().enumerate() {
            let (w, h, data) = decode_rgb8(bytes).map_err(|e| PatError::Format { reason: format!("cube map face {i}: {e}") })?;
            if w != h || (i > 0 && w != size) {
                return Err(PatError::Format { reason: format!("cube map face {i} is {w}x{h}, expected square faces of one size") });
            }
            size = w;
            faces.push(data);
        }
        Ok(CubeMap::from_faces(size, faces, None))
    }

    /// A constant environment: every lookup returns `c` (so the IBL terms see a uniform sky).
    pub fn constant(c: [f64; 3]) -> Result<CubeMap, PatError> {
        if !c.iter().all(|x| x.is_finite() && *x >= 0.0) {
            return Err(PatError::InvalidOption { what: format!("constant environment colour {c:?} (finite and >= 0)") });
        }
        Ok(CubeMap::from_faces(1, vec![vec![c]; 6], Some(c)))
    }

    /// Faces as RGB texels in `[0, 1]` (`size * size` each, +X -X +Y -Y +Z -Z).
    pub fn from_faces(size: usize, faces: Vec<Vec<[f64; 3]>>, constant: Option<[f64; 3]>) -> CubeMap {
        let mut m = CubeMap { size, faces, constant, avg: [[0.0; 3]; 6] };
        let dirs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
        for (a, d) in m.avg.iter_mut().zip(dirs) {
            *a = Self::average_at(&m.faces, m.size, m.constant, d);
        }
        m
    }

    fn average_at(faces: &[Vec<[f64; 3]>], size: usize, constant: Option<[f64; 3]>, tc: [f64; 3]) -> [f64; 3] {
        const D0: f64 = 0.35;
        const D1: f64 = 0.12;
        let off = |t: [f64; 3], o: [f64; 3]| [t[0] + o[0], t[1] + o[1], t[2] + o[2]];
        let offsets = [
            [0.0, 0.0, 0.0],
            [D0, D0, D0],
            [-D0, -D0, -D0],
            [-D0, D0, -D0],
            [D0, -D0, D0],
            [D1, D1, D1],
            [-D1, -D1, -D1],
            [-D1, D1, -D1],
            [D1, -D1, D1],
        ];
        let mut s = [0.0; 3];
        for o in offsets {
            let c = sample(faces, size, constant, off(tc, o));
            for a in 0..3 {
                s[a] += c[a];
            }
        }
        [s[0] * 0.111111111, s[1] * 0.111111111, s[2] * 0.111111111]
    }

    /// `texture(cubemap, dir)`: bilinear, clamped to the face edge.
    pub fn texture(&self, dir: [f64; 3]) -> [f64; 3] {
        sample(&self.faces, self.size, self.constant, dir)
    }

    /// The demo's `textureBlurred(cubemap, tc)`: the six averaged directions weighted by
    /// `((tc·axis)/2 + 1/2)^4`.
    pub fn blurred(&self, tc: [f64; 3]) -> [f64; 3] {
        // order in `avg`: +X +Y +Z -X -Y -Z (r t f l b a)
        let [r, t, f, l, b, a] = self.avg;
        let kr = tc[0] * 0.5 + 0.5;
        let kt = tc[1] * 0.5 + 0.5;
        let kf = tc[2] * 0.5 + 0.5;
        let (kl, kb, ka) = (1.0 - kr, 1.0 - kt, 1.0 - kf);
        let somestep = |x: f64| x.powi(4);
        let (kr, kt, kf, kl, kb, ka) = (somestep(kr), somestep(kt), somestep(kf), somestep(kl), somestep(kb), somestep(ka));
        let mut ret = [0.0; 3];
        let mut d;
        // the demo's accumulation order: f, a, l, r, t, b
        for (c, k) in [(f, kf), (a, ka), (l, kl), (r, kr), (t, kt), (b, kb)] {
            for i in 0..3 {
                ret[i] += c[i] * k;
            }
        }
        d = kf;
        d += ka;
        d += kl;
        d += kr;
        d += kt;
        d += kb;
        [ret[0] / d, ret[1] / d, ret[2] / d]
    }
}

/// OpenGL cube-face selection and bilinear lookup with `GL_CLAMP_TO_EDGE`.
fn sample(faces: &[Vec<[f64; 3]>], size: usize, constant: Option<[f64; 3]>, d: [f64; 3]) -> [f64; 3] {
    if let Some(c) = constant {
        return c;
    }
    let (ax, ay, az) = (d[0].abs(), d[1].abs(), d[2].abs());
    // (face, sc, tc, ma) from the OpenGL specification, table 8.19
    let (face, sc, tc, ma) = if ax >= ay && ax >= az {
        if d[0] > 0.0 { (0, -d[2], -d[1], ax) } else { (1, d[2], -d[1], ax) }
    } else if ay >= az {
        if d[1] > 0.0 { (2, d[0], d[2], ay) } else { (3, d[0], -d[2], ay) }
    } else if d[2] > 0.0 {
        (4, d[0], -d[1], az)
    } else {
        (5, -d[0], -d[1], az)
    };
    if ma == 0.0 || !ma.is_finite() {
        return [0.0; 3];
    }
    let s = 0.5 * (sc / ma + 1.0);
    let t = 0.5 * (tc / ma + 1.0);
    let n = size as f64;
    let (u, v) = (s * n - 0.5, t * n - 0.5);
    let (u0, v0) = (u.floor(), v.floor());
    let (fu, fv) = (u - u0, v - v0);
    let at = |i: f64, j: f64| {
        let ii = (i.max(0.0) as usize).min(size - 1);
        let jj = (j.max(0.0) as usize).min(size - 1);
        faces[face][jj * size + ii]
    };
    let (c00, c10, c01, c11) = (at(u0, v0), at(u0 + 1.0, v0), at(u0, v0 + 1.0), at(u0 + 1.0, v0 + 1.0));
    let mut out = [0.0; 3];
    for a in 0..3 {
        let top = c00[a] * (1.0 - fu) + c10[a] * fu;
        let bot = c01[a] * (1.0 - fu) + c11[a] * fu;
        out[a] = top * (1.0 - fv) + bot * fv;
    }
    out
}

/// Decode an 8-bit RGB or RGBA PNG into `[0, 1]` texels.
fn decode_rgb8(bytes: &[u8]) -> Result<(usize, usize, Vec<[f64; 3]>), String> {
    let mut dec = png::Decoder::new(std::io::Cursor::new(bytes));
    dec.set_transformations(png::Transformations::EXPAND);
    let mut reader = dec.read_info().map_err(|e| e.to_string())?;
    let mut buf = vec![0u8; reader.output_buffer_size().ok_or("image too large")?];
    let info = reader.next_frame(&mut buf).map_err(|e| e.to_string())?;
    if info.bit_depth != png::BitDepth::Eight {
        return Err(format!("bit depth {:?}, expected 8", info.bit_depth));
    }
    let ch = match info.color_type {
        png::ColorType::Rgb => 3,
        png::ColorType::Rgba => 4,
        other => return Err(format!("colour type {other:?}, expected RGB or RGBA")),
    };
    let (w, h) = (info.width as usize, info.height as usize);
    let texels = (0..w * h).map(|i| [0, 1, 2].map(|a| buf[i * ch + a] as f64 / 255.0)).collect();
    Ok((w, h, texels))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forest_loads_and_samples_in_range() {
        let m = CubeMap::forest().unwrap();
        assert_eq!(m.size, 64);
        for d in [[1.0, 0.0, 0.0], [0.3, -0.7, 0.2], [0.0, 0.0, -1.0], [-0.5, 0.5, 0.5]] {
            let c = m.texture(d);
            let b = m.blurred(d);
            assert!(c.iter().chain(b.iter()).all(|x| (0.0..=1.0).contains(x)), "{c:?} {b:?}");
        }
    }

    #[test]
    fn constant_environment_is_constant() {
        let m = CubeMap::constant([0.25, 0.5, 0.75]).unwrap();
        for d in [[1.0, 0.0, 0.0], [0.3, -0.7, 0.2]] {
            assert_eq!(m.texture(d), [0.25, 0.5, 0.75]);
            let b = m.blurred(d);
            for (x, y) in b.iter().zip([0.25, 0.5, 0.75]) {
                assert!((x - y).abs() < 1e-8, "{b:?}");
            }
        }
        assert!(CubeMap::constant([f64::NAN, 0.0, 0.0]).is_err());
    }
}
