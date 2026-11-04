//! Perlin's improved noise (K. Perlin, "Improving Noise", SIGGRAPH 2002) in 2D and 3D, written
//! from the paper: the quintic fade `6t⁵ − 15t⁴ + 10t³`, hashed lattice gradients (the 12 edge
//! vectors of the cube in 3D, 8 directions in 2D), trilinear/bilinear blending, and fractal sums
//! with octaves, persistence and lacunarity.
//!
//! The permutation table is not Perlin's published one: it is a Fisher-Yates shuffle of `0..256`
//! drawn from [`Rng`] seeded by `base`, so each `base` gives a different but reproducible field
//! (the reference uses `base = seed % 256`; NOTES-REFERENCE §11). Output lies in `[-1, 1]`.

use super::rng::Rng;

/// A noise field (one permutation table).
#[derive(Debug, Clone)]
pub struct Perlin {
    perm: [u8; 512],
}

fn fade(t: f64) -> f64 {
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

fn lerp(t: f64, a: f64, b: f64) -> f64 {
    a + t * (b - a)
}

/// Gradient of Perlin's improved noise: `hash & 15` selects one of 12 edge directions
/// (four of them repeated).
fn grad3(h: u8, x: f64, y: f64, z: f64) -> f64 {
    let h = h & 15;
    let u = if h < 8 { x } else { y };
    let v = if h < 4 {
        y
    } else if h == 12 || h == 14 {
        x
    } else {
        z
    };
    (if h & 1 == 0 { u } else { -u }) + (if h & 2 == 0 { v } else { -v })
}

fn grad2(h: u8, x: f64, y: f64) -> f64 {
    const S: f64 = std::f64::consts::FRAC_1_SQRT_2;
    let (gx, gy) = match h & 7 {
        0 => (1.0, 0.0),
        1 => (-1.0, 0.0),
        2 => (0.0, 1.0),
        3 => (0.0, -1.0),
        4 => (S, S),
        5 => (-S, S),
        6 => (S, -S),
        _ => (-S, -S),
    };
    gx * x + gy * y
}

impl Perlin {
    /// The field number `base` (any integer; reduced mod 256 as in the reference).
    pub fn new(base: u64) -> Self {
        let mut p: Vec<u8> = (0..=255u8).collect();
        let mut rng = Rng::new(0x5045_524C ^ (base % 256));
        rng.shuffle(&mut p);
        let mut perm = [0u8; 512];
        for i in 0..512 {
            perm[i] = p[i & 255];
        }
        Perlin { perm }
    }

    /// One octave of 3D noise, in `[-1, 1]`, zero at lattice points.
    pub fn noise3(&self, x: f64, y: f64, z: f64) -> f64 {
        let (fx, fy, fz) = (x.floor(), y.floor(), z.floor());
        let (xi, yi, zi) = ((fx as i64 & 255) as usize, (fy as i64 & 255) as usize, (fz as i64 & 255) as usize);
        let (x, y, z) = (x - fx, y - fy, z - fz);
        let (u, v, w) = (fade(x), fade(y), fade(z));
        let p = &self.perm;
        let a = p[xi] as usize + yi;
        let aa = p[a] as usize + zi;
        let ab = p[a + 1] as usize + zi;
        let b = p[xi + 1] as usize + yi;
        let ba = p[b] as usize + zi;
        let bb = p[b + 1] as usize + zi;
        lerp(
            w,
            lerp(
                v,
                lerp(u, grad3(p[aa], x, y, z), grad3(p[ba], x - 1.0, y, z)),
                lerp(u, grad3(p[ab], x, y - 1.0, z), grad3(p[bb], x - 1.0, y - 1.0, z)),
            ),
            lerp(
                v,
                lerp(u, grad3(p[aa + 1], x, y, z - 1.0), grad3(p[ba + 1], x - 1.0, y, z - 1.0)),
                lerp(u, grad3(p[ab + 1], x, y - 1.0, z - 1.0), grad3(p[bb + 1], x - 1.0, y - 1.0, z - 1.0)),
            ),
        )
    }

    /// One octave of 2D noise, in `[-1, 1]`.
    pub fn noise2(&self, x: f64, y: f64) -> f64 {
        let (fx, fy) = (x.floor(), y.floor());
        let (xi, yi) = ((fx as i64 & 255) as usize, (fy as i64 & 255) as usize);
        let (x, y) = (x - fx, y - fy);
        let (u, v) = (fade(x), fade(y));
        let p = &self.perm;
        let a = p[xi] as usize + yi;
        let b = p[xi + 1] as usize + yi;
        lerp(
            v,
            lerp(u, grad2(p[a], x, y), grad2(p[b], x - 1.0, y)),
            lerp(u, grad2(p[a + 1], x, y - 1.0), grad2(p[b + 1], x - 1.0, y - 1.0)),
        )
    }

    /// Fractal 3D noise: `Σ persistenceⁱ noise(lacunarityⁱ x) / Σ persistenceⁱ` over `octaves ≥ 1`.
    pub fn fbm3(&self, p: [f64; 3], octaves: usize, persistence: f64, lacunarity: f64) -> f64 {
        let (mut freq, mut amp, mut sum, mut norm) = (1.0, 1.0, 0.0, 0.0);
        for _ in 0..octaves.max(1) {
            sum += amp * self.noise3(p[0] * freq, p[1] * freq, p[2] * freq);
            norm += amp;
            freq *= lacunarity;
            amp *= persistence;
        }
        sum / norm
    }

    /// Fractal 2D noise (see [`Perlin::fbm3`]).
    pub fn fbm2(&self, p: [f64; 2], octaves: usize, persistence: f64, lacunarity: f64) -> f64 {
        let (mut freq, mut amp, mut sum, mut norm) = (1.0, 1.0, 0.0, 0.0);
        for _ in 0..octaves.max(1) {
            sum += amp * self.noise2(p[0] * freq, p[1] * freq);
            norm += amp;
            freq *= lacunarity;
            amp *= persistence;
        }
        sum / norm
    }
}
