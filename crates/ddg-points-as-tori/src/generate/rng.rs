//! A small deterministic random generator (xoshiro256** seeded by SplitMix64; both public domain
//! algorithms by Blackman and Vigna). Every generator of this module takes a seed and draws from
//! its own stream, so results never depend on thread counts or on the order of other items
//! (the reference shares one `mt19937` across OpenMP threads, a data race; NOTES-REFERENCE §11).

/// SplitMix64 step (also used as a 64-bit hash).
pub fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// xoshiro256** generator.
#[derive(Debug, Clone, PartialEq)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    /// A generator from a seed.
    pub fn new(seed: u64) -> Self {
        let mut st = seed;
        let s = [splitmix64(&mut st), splitmix64(&mut st), splitmix64(&mut st), splitmix64(&mut st)];
        Rng { s }
    }

    /// An independent stream for item `index` of the run `seed`.
    pub fn for_item(seed: u64, index: u64) -> Self {
        let mut st = seed ^ 0xD1B5_4A32_D192_ED03;
        let a = splitmix64(&mut st);
        let mut st2 = a ^ index.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let b = splitmix64(&mut st2);
        Rng::new(b)
    }

    /// Next 64 random bits.
    pub fn next_u64(&mut self) -> u64 {
        let r = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        r
    }

    /// Uniform in `[0, 1)` (53 bits).
    pub fn f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform in `[a, b)` (`a` and `b` may be given in either order).
    pub fn uniform(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.f64()
    }

    /// Uniform integer in `0..n` (`n = 0` gives 0).
    pub fn below(&mut self, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        // Lemire's multiply-shift without the rejection step (bias < n / 2^64).
        (((self.next_u64() as u128) * (n as u128)) >> 64) as usize
    }

    /// Standard normal (Box-Muller).
    pub fn normal(&mut self) -> f64 {
        let u1 = 1.0 - self.f64();
        let u2 = self.f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// A uniformly random direction on the circle.
    pub fn unit_circle(&mut self) -> [f64; 2] {
        let a = self.uniform(0.0, std::f64::consts::TAU);
        [a.cos(), a.sin()]
    }

    /// A uniformly random direction on the sphere.
    pub fn unit_sphere(&mut self) -> [f64; 3] {
        let z = self.uniform(-1.0, 1.0);
        let a = self.uniform(0.0, std::f64::consts::TAU);
        let r = (1.0 - z * z).max(0.0).sqrt();
        [r * a.cos(), r * a.sin(), z]
    }

    /// A uniformly random rotation matrix (rows), from a uniform unit quaternion (Shoemake).
    pub fn rotation(&mut self) -> [[f64; 3]; 3] {
        let (u1, u2, u3) = (self.f64(), self.f64(), self.f64());
        let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
        let t = std::f64::consts::TAU;
        let (x, y, z, w) = (a * (t * u2).sin(), a * (t * u2).cos(), b * (t * u3).sin(), b * (t * u3).cos());
        [
            [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
            [2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
            [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)],
        ]
    }

    /// Fisher-Yates shuffle.
    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = self.below(i + 1);
            v.swap(i, j);
        }
    }
}
