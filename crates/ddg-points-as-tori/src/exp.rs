//! A crate-owned `exp` for the blend (SPEC.md §18, WP4 item 3 / WP4b).
//!
//! The blend's weights `e^{−λ(ρ−s)}` are the only transcendental function of the f64 evaluation.
//! glibc's `exp` and CUDA's `exp` differ by up to 1 ulp on about 6 % of arguments, which made the
//! GPU blend level 1 (tolerance) instead of level 2 (bitwise = CPU). This `exp` is one fixed
//! sequence of IEEE operations (`+ − × ÷`, `rint`, `floor`-free integer casts) with no fused
//! multiply-add: the same text is in `gpu/kernels.cu` (`pat_exp`, NVRTC with `--fmad=false`), so
//! CPU and GPU return the same bits.
//!
//! **Method.** `k = rint(x/ln 2)`, `r = (x − k·ln2_hi) − k·ln2_lo` (Cody–Waite; `ln2_hi` has 21
//! trailing zero bits so `k·ln2_hi` is exact for `|k| ≤ 1075`), `|r| ≤ ½ ln 2`;
//! `e^r = 1 + (r + r²·P(r))` with `P` the Taylor polynomial `Σ_{n=2}^{13} rⁿ⁻²/n!` in Horner form
//! (truncation error `< 5e-18`, below 0.05 ulp); then multiplication by `2^k` (exact for normal
//! results, a single rounding in the subnormal range).
//!
//! **Accuracy** (measured in `examples/exp_check.rs` against glibc, whose own error is `< 0.52`
//! ulp, over 10⁸ arguments in `[−745, 709]` and the blend's range `[−800, 0]` plus edge cases):
//! see `NOTES-WP15.md`. Edge cases: NaN → NaN, `x > 709.782712893384` → +∞,
//! `x < −745.1332191019412` → +0, `−∞` → 0, `+∞` → +∞, `exp(0) = 1`.

/// `1/ln 2`.
const LOG2E: f64 = 1.442_695_040_888_963_4;
/// High part of `ln 2` (fdlibm's `ln2HI`: 21 trailing zero bits).
const LN2_HI: f64 = 6.931_471_803_691_238_164_90e-1;
/// Low part: `ln 2 − LN2_HI`.
const LN2_LO: f64 = 1.908_214_929_270_587_700_02e-10;
/// Largest argument with a finite result.
const X_MAX: f64 = 709.782_712_893_384;
/// Smallest argument with a non-zero result (below it `exp` rounds to 0).
const X_MIN: f64 = -745.133_219_101_941_2;

#[inline(always)]
fn pow2(k: i32) -> f64 {
    f64::from_bits(((k + 1023) as u64) << 52)
}

/// `e^x`, bitwise identical to `pat_exp` in `gpu/kernels.cu`.
#[inline]
pub fn exp(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x > X_MAX {
        return f64::INFINITY;
    }
    if x < X_MIN {
        return 0.0;
    }
    let kf = (x * LOG2E).round_ties_even();
    let r = (x - kf * LN2_HI) - kf * LN2_LO;
    let mut p = 1.605_904_383_682_161_3e-10;
    p = p * r + 2.087_675_698_786_81e-9;
    p = p * r + 2.505_210_838_544_172e-8;
    p = p * r + 2.755_731_922_398_589e-7;
    p = p * r + 2.755_731_922_398_589_3e-6;
    p = p * r + 2.480_158_730_158_73e-5;
    p = p * r + 1.984_126_984_126_984e-4;
    p = p * r + 1.388_888_888_888_889e-3;
    p = p * r + 8.333_333_333_333_333e-3;
    p = p * r + 4.166_666_666_666_666_4e-2;
    p = p * r + 0.166_666_666_666_666_66;
    p = p * r + 0.5;
    let y = 1.0 + (r + r * r * p);
    let k = kf as i32;
    if k > 1023 {
        (y * 2.0) * pow2(1023)
    } else if k < -1021 {
        (y * pow2(k + 1000)) * pow2(-1000)
    } else {
        y * pow2(k)
    }
}

#[cfg(test)]
mod tests {
    use super::exp;

    fn ulps(a: f64, b: f64) -> f64 {
        if a == b {
            return 0.0;
        }
        (a.to_bits() as i64 - b.to_bits() as i64).unsigned_abs() as f64
    }

    #[test]
    fn edge_cases() {
        assert_eq!(exp(0.0), 1.0);
        assert_eq!(exp(-0.0), 1.0);
        assert!(exp(f64::NAN).is_nan());
        assert_eq!(exp(f64::INFINITY), f64::INFINITY);
        assert_eq!(exp(f64::NEG_INFINITY), 0.0);
        assert_eq!(exp(710.0), f64::INFINITY);
        assert_eq!(exp(-746.0), 0.0);
        assert!(exp(709.78).is_finite());
        assert!(exp(-745.0) > 0.0);
        assert!(exp(-1e-300) == 1.0);
    }

    #[test]
    fn close_to_libm() {
        let mut worst = 0.0f64;
        let mut s = 12345u64;
        for _ in 0..200_000 {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let x = ((s >> 11) as f64 / (1u64 << 53) as f64) * 800.0 - 750.0;
            if x > 700.0 {
                continue;
            }
            worst = worst.max(ulps(exp(x), x.exp()));
        }
        assert!(worst <= 2.0, "worst {worst} ulp");
    }
}
