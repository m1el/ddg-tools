//! Images of 2D slices: a small raster with drawing helpers, an uncompressed PNG writer, and SVG
//! contours. No dependencies.

use super::grid::Polyline;
use crate::PatError;
use std::io::Write;
use std::path::Path;

/// An RGB raster whose pixels map to a world box (`y` up).
#[derive(Debug, Clone, PartialEq)]
pub struct Raster {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// `[[xmin, ymin], [xmax, ymax]]` covered by the pixel centers' outer edges.
    pub bbox: [[f64; 2]; 2],
    /// `3 * width * height` bytes, top row first.
    pub rgb: Vec<u8>,
}

impl Raster {
    /// A raster filled with `bg`.
    pub fn new(width: usize, height: usize, bbox: [[f64; 2]; 2], bg: [u8; 3]) -> Result<Raster, PatError> {
        if width == 0 || height == 0 || width.saturating_mul(height) > 1 << 28 || bbox.iter().flatten().any(|c| !c.is_finite()) || bbox[1][0] <= bbox[0][0] || bbox[1][1] <= bbox[0][1] {
            return Err(PatError::InvalidOption { what: format!("Raster {width} x {height} on {bbox:?}") });
        }
        Ok(Raster { width, height, bbox, rgb: bg.iter().copied().cycle().take(3 * width * height).collect() })
    }

    /// World position of the center of pixel `(px, py)` (`py = 0` is the top row).
    pub fn pixel_center(&self, px: usize, py: usize) -> [f64; 2] {
        let fx = (px as f64 + 0.5) / self.width as f64;
        let fy = 1.0 - (py as f64 + 0.5) / self.height as f64;
        [self.bbox[0][0] + fx * (self.bbox[1][0] - self.bbox[0][0]), self.bbox[0][1] + fy * (self.bbox[1][1] - self.bbox[0][1])]
    }

    /// Pixel coordinates (continuous) of a world point.
    pub fn to_pixel(&self, p: [f64; 2]) -> [f64; 2] {
        [
            (p[0] - self.bbox[0][0]) / (self.bbox[1][0] - self.bbox[0][0]) * self.width as f64,
            (1.0 - (p[1] - self.bbox[0][1]) / (self.bbox[1][1] - self.bbox[0][1])) * self.height as f64,
        ]
    }

    /// Set pixel colors from a function of the world position.
    pub fn fill_with(&mut self, f: impl Fn([f64; 2]) -> [u8; 3]) {
        for py in 0..self.height {
            for px in 0..self.width {
                let c = f(self.pixel_center(px, py));
                self.set(px as i64, py as i64, c, 1.0);
            }
        }
    }

    /// Blend a color into a pixel (`alpha` in `[0, 1]`); outside pixels are ignored.
    pub fn set(&mut self, px: i64, py: i64, c: [u8; 3], alpha: f64) {
        if px < 0 || py < 0 || px as usize >= self.width || py as usize >= self.height {
            return;
        }
        let o = 3 * (py as usize * self.width + px as usize);
        for (dst, &src) in self.rgb[o..o + 3].iter_mut().zip(&c) {
            *dst = (*dst as f64 * (1.0 - alpha) + src as f64 * alpha).round() as u8;
        }
    }

    /// A filled disc of `radius` pixels at a world point.
    pub fn dot(&mut self, p: [f64; 2], radius: f64, c: [u8; 3]) {
        let q = self.to_pixel(p);
        let r = radius.max(0.5);
        for py in (q[1] - r - 1.0).floor() as i64..=(q[1] + r + 1.0).ceil() as i64 {
            for px in (q[0] - r - 1.0).floor() as i64..=(q[0] + r + 1.0).ceil() as i64 {
                let d = ((px as f64 + 0.5 - q[0]).powi(2) + (py as f64 + 0.5 - q[1]).powi(2)).sqrt();
                let a = (r + 0.5 - d).clamp(0.0, 1.0);
                if a > 0.0 {
                    self.set(px, py, c, a);
                }
            }
        }
    }

    /// An anti-aliased line of `width` pixels between world points (clipped to the raster).
    pub fn line(&mut self, a: [f64; 2], b: [f64; 2], width: f64, c: [u8; 3]) {
        let (pa, pb) = (self.to_pixel(a), self.to_pixel(b));
        let len = ((pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2)).sqrt();
        if !len.is_finite() || len > 1e6 {
            return;
        }
        let steps = (len * 2.0).ceil().max(1.0) as usize;
        let half = 0.5 * width;
        for s in 0..=steps {
            let t = s as f64 / steps as f64;
            let q = [pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1])];
            for py in (q[1] - half - 1.0).floor() as i64..=(q[1] + half + 1.0).ceil() as i64 {
                for px in (q[0] - half - 1.0).floor() as i64..=(q[0] + half + 1.0).ceil() as i64 {
                    let d = ((px as f64 + 0.5 - q[0]).powi(2) + (py as f64 + 0.5 - q[1]).powi(2)).sqrt();
                    let al = (half + 0.5 - d).clamp(0.0, 1.0) * 0.5;
                    if al > 0.0 {
                        self.set(px, py, c, al);
                    }
                }
            }
        }
    }

    /// A polyline.
    pub fn polyline(&mut self, p: &Polyline, width: f64, c: [u8; 3]) {
        for w in p.points.windows(2) {
            self.line(w[0], w[1], width, c);
        }
        if p.closed && p.points.len() > 2 {
            self.line(p.points[p.points.len() - 1], p.points[0], width, c);
        }
    }

    /// A circle outline (64 to 720 chords by size), optionally only the arc of half-angle `half_arc`
    /// around the direction `dir`.
    pub fn circle_arc(&mut self, center: [f64; 2], radius: f64, dir: [f64; 2], half_arc: f64, width: f64, c: [u8; 3]) {
        let rp = radius / (self.bbox[1][0] - self.bbox[0][0]) * self.width as f64;
        let n = ((rp * half_arc.min(std::f64::consts::PI)).ceil() as usize).clamp(16, 720);
        let a0 = dir[1].atan2(dir[0]);
        let pt = |t: f64| [center[0] + radius * (a0 + t).cos(), center[1] + radius * (a0 + t).sin()];
        for s in 0..n {
            let t0 = -half_arc + 2.0 * half_arc * s as f64 / n as f64;
            let t1 = -half_arc + 2.0 * half_arc * (s + 1) as f64 / n as f64;
            self.line(pt(t0), pt(t1), width, c);
        }
    }
}

fn crc32(data: &[u8], table: &[u32; 256], mut crc: u32) -> u32 {
    crc = !crc;
    for &b in data {
        crc = table[((crc ^ b as u32) & 0xff) as usize] ^ (crc >> 8);
    }
    !crc
}

/// Write an RGB image as a PNG (stored, i.e. uncompressed, deflate blocks).
pub fn write_png(path: &Path, raster: &Raster) -> Result<(), PatError> {
    let io = |e: std::io::Error| PatError::Format { reason: format!("{}: {e}", path.display()) };
    let table: [u32; 256] = std::array::from_fn(|n| {
        let mut c = n as u32;
        for _ in 0..8 {
            c = if c & 1 == 1 { 0xEDB8_8320 ^ (c >> 1) } else { c >> 1 };
        }
        c
    });
    let mut raw = Vec::with_capacity((3 * raster.width + 1) * raster.height);
    for row in raster.rgb.chunks(3 * raster.width) {
        raw.push(0);
        raw.extend_from_slice(row);
    }
    let mut z = vec![0x78, 0x01];
    let mut chunks = raw.chunks(65535).peekable();
    while let Some(ch) = chunks.next() {
        z.push(u8::from(chunks.peek().is_none()));
        z.extend_from_slice(&(ch.len() as u16).to_le_bytes());
        z.extend_from_slice(&(!(ch.len() as u16)).to_le_bytes());
        z.extend_from_slice(ch);
    }
    let (mut a, mut b) = (1u32, 0u32);
    for &x in &raw {
        a = (a + x as u32) % 65521;
        b = (b + a) % 65521;
    }
    z.extend_from_slice(&((b << 16) | a).to_be_bytes());
    let mut out = vec![0x89, b'P', b'N', b'G', 0x0d, 0x0a, 0x1a, 0x0a];
    let mut chunk = |kind: &[u8; 4], data: &[u8]| {
        out.extend_from_slice(&(data.len() as u32).to_be_bytes());
        out.extend_from_slice(kind);
        out.extend_from_slice(data);
        let c = crc32(data, &table, crc32(kind, &table, 0));
        out.extend_from_slice(&c.to_be_bytes());
    };
    let mut ihdr = Vec::new();
    ihdr.extend_from_slice(&(raster.width as u32).to_be_bytes());
    ihdr.extend_from_slice(&(raster.height as u32).to_be_bytes());
    ihdr.extend_from_slice(&[8, 2, 0, 0, 0]);
    chunk(b"IHDR", &ihdr);
    chunk(b"IDAT", &z);
    chunk(b"IEND", &[]);
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(io)?;
    }
    std::fs::File::create(path).and_then(|mut f| f.write_all(&out)).map_err(io)
}

/// Write contours (and optional points) as an SVG of `width` pixels covering `bbox`.
pub fn write_svg_contours(path: &Path, bbox: [[f64; 2]; 2], width: usize, contours: &[Polyline], points: &[[f64; 2]]) -> Result<(), PatError> {
    let io = |e: std::io::Error| PatError::Format { reason: format!("{}: {e}", path.display()) };
    if bbox.iter().flatten().any(|c| !c.is_finite()) || bbox[1][0] <= bbox[0][0] || bbox[1][1] <= bbox[0][1] || width == 0 {
        return Err(PatError::InvalidOption { what: format!("write_svg_contours: box {bbox:?}, width {width}") });
    }
    let sx = width as f64 / (bbox[1][0] - bbox[0][0]);
    let height = ((bbox[1][1] - bbox[0][1]) * sx).ceil();
    let tp = |p: [f64; 2]| format!("{:.3},{:.3}", (p[0] - bbox[0][0]) * sx, (bbox[1][1] - p[1]) * sx);
    let mut s = format!("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\">\n<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");
    for c in contours {
        let pts: Vec<String> = c.points.iter().map(|&p| tp(p)).collect();
        s.push_str(&format!(
            "<{} points=\"{}\" fill=\"none\" stroke=\"black\" stroke-width=\"1.5\"/>\n",
            if c.closed { "polygon" } else { "polyline" },
            pts.join(" ")
        ));
    }
    for &p in points {
        let q: Vec<String> = tp(p).split(',').map(String::from).collect();
        s.push_str(&format!("<circle cx=\"{}\" cy=\"{}\" r=\"2\" fill=\"crimson\"/>\n", q[0], q[1]));
    }
    s.push_str("</svg>\n");
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(io)?;
    }
    std::fs::write(path, s).map_err(io)
}
