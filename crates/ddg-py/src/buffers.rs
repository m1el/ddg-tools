//! Decoding the memory of a Python buffer into typed vectors, without Python.
//!
//! The PyO3 glue hands every array argument over as a [`BufferView`]: the format string, item
//! size, shape and bytes that the buffer protocol reports. Everything else — which element types
//! a parameter accepts, how rows are counted, which values are representable — is decided here, so
//! it can be tested with `cargo test`.
//!
//! # Accepted element types
//!
//! | parameter kind | accepted `struct` formats |
//! |----------------|---------------------------|
//! | real numbers ([`decode_reals`]) | `f` (float32), `d` (float64) |
//! | indices ([`decode_indices`]) | `h H i I l L q Q n N` with an item size of 2, 4 or 8 bytes |
//! | either | `B`, `b`, `c` (`bytes`, `bytearray`, uint8 arrays): **raw bytes** of native-endian float32 (real numbers) or int32 (indices) |
//!
//! A format may carry the byte-order prefix `@`, `=` or the one of the machine (`<` on x86-64 and
//! ARM64). Buffers in the other byte order are rejected. One-byte integers are never read as
//! index values: a one-byte buffer always means raw bytes, so that `bytes` works and a `uint8`
//! array is not silently misread.
//!
//! # Shapes
//!
//! A parameter with `width` items per row (3 per vertex position, 2 per seam edge, ...) accepts
//! a flat buffer whose length is a multiple of `width`, or a two-dimensional one of shape
//! `(n, width)`. Buffers must be C-contiguous; the glue checks that before building the view.

use crate::error::InputError;

/// A borrowed, C-contiguous Python buffer as the buffer protocol describes it.
#[derive(Debug, Clone, Copy)]
pub struct BufferView<'a> {
    /// Name of the Python parameter, for error messages.
    pub name: &'static str,
    /// `struct`-module format string of one item (`"f"`, `"<i"`, `"B"`, ...).
    pub format: &'a str,
    /// Size of one item in bytes.
    pub item_size: usize,
    /// Shape of the buffer; one entry for a flat buffer.
    pub shape: &'a [usize],
    /// The memory, `item_size * product(shape)` bytes in C order.
    pub bytes: &'a [u8],
}

/// The element types this module can read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Element {
    F32,
    F64,
    Signed(usize),
    Unsigned(usize),
    RawBytes,
}

const REALS_EXPECTED: &str =
    "float32 ('f') or float64 ('d'), or raw bytes of native-endian float32";
const INDICES_EXPECTED: &str = "a 16-, 32- or 64-bit integer type ('h', 'H', 'i', 'I', 'l', 'L', 'q', 'Q'), \
     or raw bytes of native-endian int32";

/// Parses a single-item `struct` format. `None` means "not an element type this module reads".
fn parse_format(format: &str, item_size: usize) -> Option<Element> {
    let bytes = format.as_bytes();
    let native_order = if cfg!(target_endian = "little") {
        b'<'
    } else {
        b'>'
    };
    let code = match bytes {
        [c] => *c,
        [b'@' | b'=', c] => *c,
        [order, c] if *order == native_order => *c,
        _ => return None,
    };
    match (code, item_size) {
        (b'f', 4) => Some(Element::F32),
        (b'd', 8) => Some(Element::F64),
        (b'B' | b'b' | b'c', 1) => Some(Element::RawBytes),
        (b'h' | b'i' | b'l' | b'q' | b'n', 2 | 4 | 8) => Some(Element::Signed(item_size)),
        (b'H' | b'I' | b'L' | b'Q' | b'N', 2 | 4 | 8) => Some(Element::Unsigned(item_size)),
        _ => None,
    }
}

/// Checks the shape of `view` against `width` items per row. Returns the number of items.
fn check_shape(
    view: &BufferView<'_>,
    element_size: usize,
    width: usize,
    row: &'static str,
) -> Result<usize, InputError> {
    debug_assert!(width > 0 && element_size > 0);
    let len = view.bytes.len() / element_size;
    match view.shape {
        [] | [_] => {}
        [_, w] if *w == width => {}
        _ => {
            return Err(InputError::BadShape {
                name: view.name,
                shape: view.shape.to_vec(),
                width,
            });
        }
    }
    if !len.is_multiple_of(width) {
        return Err(InputError::LengthNotMultiple {
            name: view.name,
            len,
            width,
            row,
        });
    }
    Ok(len)
}

fn unsupported(view: &BufferView<'_>, expected: &'static str) -> InputError {
    InputError::UnsupportedFormat {
        name: view.name,
        format: view.format.to_owned(),
        item_size: view.item_size,
        expected,
    }
}

/// For raw byte buffers: the bytes stand for 4-byte elements and must be flat.
fn check_raw_bytes(view: &BufferView<'_>, element: &'static str) -> Result<(), InputError> {
    if view.shape.len() > 1 {
        return Err(InputError::BadShape {
            name: view.name,
            shape: view.shape.to_vec(),
            width: 1,
        });
    }
    if !view.bytes.len().is_multiple_of(4) {
        return Err(InputError::RawByteLength {
            name: view.name,
            len_bytes: view.bytes.len(),
            element,
        });
    }
    Ok(())
}

/// Reads a buffer of real numbers with `width` items per row as `f64`.
///
/// `row` names one row in error messages ("vertex", "cone"). Values are converted exactly
/// (float32 widens to float64); finiteness is not checked here.
///
/// # Errors
/// [`InputError::UnsupportedFormat`], [`InputError::BadShape`],
/// [`InputError::LengthNotMultiple`], [`InputError::RawByteLength`].
pub fn decode_reals(
    view: &BufferView<'_>,
    width: usize,
    row: &'static str,
) -> Result<Vec<f64>, InputError> {
    let element = parse_format(view.format, view.item_size)
        .ok_or_else(|| unsupported(view, REALS_EXPECTED))?;
    match element {
        Element::F64 => {
            check_shape(view, 8, width, row)?;
            Ok(view
                .bytes
                .chunks_exact(8)
                .map(|c| f64::from_ne_bytes(c.try_into().expect("chunk of 8 bytes")))
                .collect())
        }
        Element::F32 | Element::RawBytes => {
            if element == Element::RawBytes {
                check_raw_bytes(view, "float32")?;
            }
            check_shape(view, 4, width, row)?;
            Ok(view
                .bytes
                .chunks_exact(4)
                .map(|c| f64::from(f32::from_ne_bytes(c.try_into().expect("chunk of 4 bytes"))))
                .collect())
        }
        Element::Signed(_) | Element::Unsigned(_) => Err(unsupported(view, REALS_EXPECTED)),
    }
}

/// Reads a buffer of indices with `width` items per row as `u32`.
///
/// Every value must be non-negative and at most `u32::MAX − 1` (`u32::MAX` is the Rust side's
/// "no element" sentinel). Range checks against the mesh happen later, in
/// [`crate::contract::UnwrapInput::from_buffers`].
///
/// # Errors
/// [`InputError::UnsupportedFormat`], [`InputError::BadShape`],
/// [`InputError::LengthNotMultiple`], [`InputError::RawByteLength`],
/// [`InputError::NegativeIndex`], [`InputError::IndexTooLarge`].
pub fn decode_indices(
    view: &BufferView<'_>,
    width: usize,
    row: &'static str,
) -> Result<Vec<u32>, InputError> {
    let element = parse_format(view.format, view.item_size)
        .ok_or_else(|| unsupported(view, INDICES_EXPECTED))?;
    let (size, signed) = match element {
        Element::Signed(size) => (size, true),
        Element::Unsigned(size) => (size, false),
        Element::RawBytes => {
            check_raw_bytes(view, "int32")?;
            (4, true)
        }
        Element::F32 | Element::F64 => return Err(unsupported(view, INDICES_EXPECTED)),
    };
    let len = check_shape(view, size, width, row)?;
    let mut out = Vec::with_capacity(len);
    for (position, chunk) in view.bytes.chunks_exact(size).enumerate() {
        // Widen to i128 so that every accepted type has one code path for the range checks.
        let value: i128 = match (size, signed) {
            (2, true) => i128::from(i16::from_ne_bytes(chunk.try_into().expect("2 bytes"))),
            (2, false) => i128::from(u16::from_ne_bytes(chunk.try_into().expect("2 bytes"))),
            (4, true) => i128::from(i32::from_ne_bytes(chunk.try_into().expect("4 bytes"))),
            (4, false) => i128::from(u32::from_ne_bytes(chunk.try_into().expect("4 bytes"))),
            (8, true) => i128::from(i64::from_ne_bytes(chunk.try_into().expect("8 bytes"))),
            (_, _) => i128::from(u64::from_ne_bytes(chunk.try_into().expect("8 bytes"))),
        };
        if value < 0 {
            return Err(InputError::NegativeIndex {
                name: view.name,
                position,
                value: value as i64,
            });
        }
        if value >= i128::from(u32::MAX) {
            return Err(InputError::IndexTooLarge {
                name: view.name,
                position,
                value: value as u64,
            });
        }
        out.push(value as u32);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn view<'a>(
        format: &'a str,
        item_size: usize,
        shape: &'a [usize],
        bytes: &'a [u8],
    ) -> BufferView<'a> {
        BufferView {
            name: "arg",
            format,
            item_size,
            shape,
            bytes,
        }
    }

    fn bytes_of<T: Copy, const N: usize>(values: &[T], to: fn(T) -> [u8; N]) -> Vec<u8> {
        values.iter().flat_map(|&v| to(v)).collect()
    }

    #[test]
    fn reals_from_f32_f64_and_raw_bytes() {
        let f32s = bytes_of(&[1.5f32, -2.25, 0.1], f32::to_ne_bytes);
        let f64s = bytes_of(&[1.5f64, -2.25, 0.1], f64::to_ne_bytes);
        let expect32 = vec![1.5, -2.25, f64::from(0.1f32)];
        assert_eq!(
            decode_reals(&view("f", 4, &[3], &f32s), 3, "vertex").unwrap(),
            expect32
        );
        assert_eq!(
            decode_reals(&view("<f", 4, &[3], &f32s), 3, "vertex").unwrap(),
            expect32
        );
        assert_eq!(
            decode_reals(&view("d", 8, &[3], &f64s), 3, "vertex").unwrap(),
            vec![1.5, -2.25, 0.1]
        );
        // `bytes`: format B, shape = byte count.
        assert_eq!(
            decode_reals(&view("B", 1, &[12], &f32s), 3, "vertex").unwrap(),
            expect32
        );
        // Two-dimensional (1, 3).
        assert_eq!(
            decode_reals(&view("d", 8, &[1, 3], &f64s), 3, "vertex").unwrap(),
            vec![1.5, -2.25, 0.1]
        );
    }

    #[test]
    fn reals_reject_wrong_types_and_shapes() {
        let f32s = bytes_of(&[1.0f32, 2.0, 3.0, 4.0], f32::to_ne_bytes);
        let e = decode_reals(&view("f", 4, &[4], &f32s), 3, "vertex").unwrap_err();
        assert_eq!(
            e.to_string(),
            "arg: 4 items is not a multiple of 3 (3 per vertex)"
        );
        let e = decode_reals(&view("f", 4, &[2, 2], &f32s), 3, "vertex").unwrap_err();
        assert_eq!(
            e.to_string(),
            "arg: buffer has shape (2, 2); expected a flat buffer or one of shape (n, 3)"
        );
        let e = decode_reals(&view("i", 4, &[4], &f32s), 1, "angle").unwrap_err();
        assert!(matches!(e, InputError::UnsupportedFormat { .. }), "{e}");
        assert!(e.to_string().contains("format \"i\""), "{e}");
        // float16 and big-endian (on little-endian machines) are rejected.
        assert!(decode_reals(&view("e", 2, &[8], &f32s), 1, "angle").is_err());
        let other_order = if cfg!(target_endian = "little") {
            ">f"
        } else {
            "<f"
        };
        assert!(decode_reals(&view(other_order, 4, &[4], &f32s), 1, "angle").is_err());
        // Raw bytes whose length is not a multiple of 4.
        let e = decode_reals(&view("B", 1, &[7], &f32s[..7]), 1, "angle").unwrap_err();
        assert_eq!(
            e.to_string(),
            "arg: a byte buffer is read as native-endian float32, but its length 7 is not a \
             multiple of 4"
        );
    }

    #[test]
    fn indices_from_every_integer_width() {
        let expect = vec![0u32, 7, 65_000];
        let cases: Vec<(&str, usize, Vec<u8>)> = vec![
            ("i", 4, bytes_of(&[0i32, 7, 65_000], i32::to_ne_bytes)),
            ("I", 4, bytes_of(&[0u32, 7, 65_000], u32::to_ne_bytes)),
            ("l", 8, bytes_of(&[0i64, 7, 65_000], i64::to_ne_bytes)),
            ("q", 8, bytes_of(&[0i64, 7, 65_000], i64::to_ne_bytes)),
            ("Q", 8, bytes_of(&[0u64, 7, 65_000], u64::to_ne_bytes)),
            ("H", 2, bytes_of(&[0u16, 7, 65_000], u16::to_ne_bytes)),
            ("=i", 4, bytes_of(&[0i32, 7, 65_000], i32::to_ne_bytes)),
            // `bytes` holding int32.
            ("B", 1, bytes_of(&[0i32, 7, 65_000], i32::to_ne_bytes)),
        ];
        for (format, size, bytes) in &cases {
            let n = bytes.len() / size;
            let shape = [n];
            assert_eq!(
                decode_indices(&view(format, *size, &shape, bytes), 3, "triangle").unwrap(),
                expect,
                "format {format}"
            );
        }
    }

    #[test]
    fn indices_reject_negative_huge_and_non_integer() {
        let neg = bytes_of(&[3i32, -1, 2], i32::to_ne_bytes);
        let e = decode_indices(&view("i", 4, &[3], &neg), 1, "vertex").unwrap_err();
        assert_eq!(e.to_string(), "arg[1] = -1 is negative");
        let huge = bytes_of(&[1u64, u64::from(u32::MAX)], u64::to_ne_bytes);
        let e = decode_indices(&view("Q", 8, &[2], &huge), 1, "vertex").unwrap_err();
        assert_eq!(
            e.to_string(),
            "arg[1] = 4294967295 does not fit a 32-bit index (limit 4294967294)"
        );
        let floats = bytes_of(&[1.0f32], f32::to_ne_bytes);
        assert!(matches!(
            decode_indices(&view("f", 4, &[1], &floats), 1, "vertex").unwrap_err(),
            InputError::UnsupportedFormat { .. }
        ));
        // A bool buffer is not an index list.
        assert!(decode_indices(&view("?", 1, &[4], &floats), 1, "vertex").is_err());
        // Shape (2, 3) for a width-2 parameter.
        let six = bytes_of(&[0i32, 1, 2, 3, 4, 5], i32::to_ne_bytes);
        let e = decode_indices(&view("i", 4, &[2, 3], &six), 2, "edge").unwrap_err();
        assert_eq!(
            e.to_string(),
            "arg: buffer has shape (2, 3); expected a flat buffer or one of shape (n, 2)"
        );
        // Three dimensions are never accepted.
        assert!(decode_indices(&view("i", 4, &[1, 2, 3], &six), 3, "triangle").is_err());
    }

    #[test]
    fn empty_buffers_decode_to_empty_vectors() {
        assert!(
            decode_reals(&view("f", 4, &[0], &[]), 3, "vertex")
                .unwrap()
                .is_empty()
        );
        assert!(
            decode_indices(&view("i", 4, &[0], &[]), 2, "edge")
                .unwrap()
                .is_empty()
        );
        assert!(
            decode_indices(&view("B", 1, &[0], &[]), 2, "edge")
                .unwrap()
                .is_empty()
        );
    }
}
