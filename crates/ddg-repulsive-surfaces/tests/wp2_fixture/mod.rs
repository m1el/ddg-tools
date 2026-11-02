//! A small reader of the `rs` golden fixtures of the Lean oracle (`lean/oracle/fixtures/rs/`,
//! protocol in `lean/oracle/README.md`, "Repulsive surfaces") for the WP2 checks of `df`,
//! `metric`, `constraints`, `saddle` and `precond`, shared by the replay of every task in
//! `tests/oracle.rs` (WP11) and by the per-task tests. The parsers take text, so that the fuzz
//! test can parse the streams of `ddg-oracle rs-gen` and `ddg-oracle rs`.

#![allow(dead_code)]

use std::collections::HashMap;
use std::path::PathBuf;

use ddg_repulsive_surfaces::constraints::Constraint;

/// The fixture directory.
pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/rs")
}

/// A hex float (16 digits, the bit pattern).
pub fn hex(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap_or_else(|_| panic!("bad hex float {s}")))
}

/// One input case.
#[derive(Debug, Clone, Default)]
pub struct Case {
    pub name: String,
    pub task: String,
    pub exponents: Option<(f64, f64)>,
    pub tau: Option<f64>,
    /// `(w_TP, c_lo, w_W)`.
    pub weights: Option<(f64, f64, f64)>,
    pub matrices_a_only: bool,
    pub eps_rel: Option<f64>,
    pub r: Vec<[f64; 3]>,
    pub u: Vec<Vec<f64>>,
    pub constraints: Vec<Constraint>,
    pub positions: Vec<[f64; 3]>,
    pub faces: Vec<[u32; 3]>,
    /// `remesh-ops`: `op NAME A B`.
    pub op: Option<(String, u32, u32)>,
    /// `remesh-ops`: pinned vertices.
    pub pinned: Vec<u32>,
    pub delaunay: bool,
    pub foldover: bool,
    /// Every line of the case as tokens (for the `step` options).
    pub lines: Vec<Vec<String>>,
}

impl Case {
    /// The tokens of the first line with key `k`.
    pub fn line(&self, k: &str) -> Option<&[String]> {
        self.lines.iter().find(|l| l[0] == k).map(|l| l.as_slice())
    }
}

fn comp(tok: &str) -> Option<u32> {
    if tok == "all" {
        None
    } else {
        Some(tok.parse().unwrap())
    }
}

/// Parses `NAME.in`.
pub fn read_case(name: &str) -> Case {
    parse_case(&std::fs::read_to_string(dir().join(format!("{name}.in"))).unwrap())
}

/// Splits a stream of `case NAME … end` blocks into `(name, text)` (the text keeps the
/// `case` line and drops `end`).
pub fn blocks(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut cur: Option<(String, String)> = None;
    for line in text.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        match t.first().copied() {
            Some("case") => cur = Some((t[1].to_string(), format!("{line}\n"))),
            Some("end") => out.extend(cur.take()),
            Some(_) => {
                if let Some((_, s)) = cur.as_mut() {
                    s.push_str(line);
                    s.push('\n');
                }
            }
            None => {}
        }
    }
    out
}

/// Parses the text of one input case.
pub fn parse_case(text: &str) -> Case {
    let mut c = Case::default();
    for line in text.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.is_empty() {
            continue;
        }
        c.lines.push(t.iter().map(|s| s.to_string()).collect());
        match t[0] {
            "op" => c.op = Some((t[1].into(), t[2].parse().unwrap(), t[3].parse().unwrap())),
            "pinned" => c.pinned = t[1..].iter().map(|s| s.parse().unwrap()).collect(),
            "legality" => {
                c.delaunay = t.contains(&"delaunay");
                c.foldover = t.contains(&"foldover");
            }
            "case" => c.name = t[1].into(),
            "task" => c.task = t[1].into(),
            "exponents" => c.exponents = Some((hex(t[1]), hex(t[2]))),
            "tau" => c.tau = Some(hex(t[1])),
            "weights" => c.weights = Some((hex(t[1]), hex(t[2]), hex(t[3]))),
            "matrices" => c.matrices_a_only = t[1] == "a",
            "eps-rel" => c.eps_rel = Some(hex(t[1])),
            "r" => c.r.push([hex(t[2]), hex(t[3]), hex(t[4])]),
            "u" => c.u.push(t[1..].iter().map(|s| hex(s)).collect()),
            "p" => c.positions.push([hex(t[1]), hex(t[2]), hex(t[3])]),
            "f" => c
                .faces
                .push([t[1].parse().unwrap(), t[2].parse().unwrap(), t[3].parse().unwrap()]),
            "c" => {
                let k = match t[1] {
                    "barycenter" => Constraint::Barycenter {
                        component: comp(t[2]),
                        target: (t.len() > 3).then(|| [hex(t[3]), hex(t[4]), hex(t[5])]),
                    },
                    "pins" => {
                        let n: usize = t[2].parse().unwrap();
                        let vertices: Vec<u32> =
                            t[3..3 + n].iter().map(|s| s.parse().unwrap()).collect();
                        let rest = &t[3 + n..];
                        let targets = (!rest.is_empty()).then(|| {
                            rest.chunks(3)
                                .map(|q| [hex(q[0]), hex(q[1]), hex(q[2])])
                                .collect()
                        });
                        Constraint::Pins { vertices, targets }
                    }
                    "area" => Constraint::Area {
                        component: comp(t[2]),
                        target: t.get(3).map(|s| hex(s)),
                    },
                    "volume" => Constraint::Volume {
                        component: comp(t[2]),
                        target: t.get(3).map(|s| hex(s)),
                    },
                    other => panic!("unknown constraint {other}"),
                };
                c.constraints.push(k);
            }
            _ => {}
        }
    }
    c
}

/// The output of a case: every line split into tokens, grouped by the first token.
#[derive(Debug, Clone, Default)]
pub struct Output {
    pub status: String,
    pub lines: HashMap<String, Vec<Vec<String>>>,
}

impl Output {
    /// All lines with key `k` (tokens after the key).
    pub fn get(&self, k: &str) -> &[Vec<String>] {
        self.lines.get(k).map_or(&[], |v| v.as_slice())
    }

    /// `error KIND …`: the kind of an error status.
    pub fn error_kind(&self) -> Option<&str> {
        self.status.strip_prefix("error ").map(|r| r.split_whitespace().next().unwrap_or(""))
    }

    /// Lines `p X Y Z` (the mesh of a `remesh-ops` answer).
    pub fn mesh_positions(&self) -> Vec<[f64; 3]> {
        self.get("p").iter().map(|l| [hex(&l[0]), hex(&l[1]), hex(&l[2])]).collect()
    }

    /// Lines `f A B C`.
    pub fn mesh_faces(&self) -> Vec<[u32; 3]> {
        self.get("f")
            .iter()
            .map(|l| [l[0].parse().unwrap(), l[1].parse().unwrap(), l[2].parse().unwrap()])
            .collect()
    }

    /// A one-line float value.
    pub fn scalar(&self, k: &str) -> f64 {
        hex(&self.get(k)[0][0])
    }

    /// Lines `KEY I X0 X1 …` as a dense row-major matrix (rows indexed by `I`).
    pub fn matrix(&self, k: &str) -> Vec<Vec<f64>> {
        let mut rows: Vec<(usize, Vec<f64>)> = self
            .get(k)
            .iter()
            .map(|l| (l[0].parse().unwrap(), l[1..].iter().map(|s| hex(s)).collect()))
            .collect();
        rows.sort_by_key(|r| r.0);
        rows.into_iter().map(|r| r.1).collect()
    }

    /// Lines `KEY V X Y Z` as per-vertex 3-vectors, flattened `3 v + axis`.
    pub fn vectors(&self, k: &str) -> Vec<f64> {
        self.matrix(k).into_iter().flatten().collect()
    }
}

/// Parses `NAME.out`.
pub fn read_output(name: &str) -> Output {
    parse_output(&std::fs::read_to_string(dir().join(format!("{name}.out"))).unwrap())
}

/// Parses the text of one output case.
pub fn parse_output(text: &str) -> Output {
    let mut o = Output::default();
    for line in text.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.is_empty() {
            continue;
        }
        if t[0] == "status" {
            o.status = t[1..].join(" ");
            continue;
        }
        o.lines
            .entry(t[0].to_string())
            .or_default()
            .push(t[1..].iter().map(|s| s.to_string()).collect());
    }
    o
}

/// The cases of `MANIFEST.tsv` with task `task` and status `ok`.
pub fn cases(task: &str) -> Vec<String> {
    let text = std::fs::read_to_string(dir().join("MANIFEST.tsv")).unwrap();
    text.lines()
        .filter_map(|l| {
            let t: Vec<&str> = l.split('\t').collect();
            (t.len() >= 3 && t[1] == task && t[2] == "ok").then(|| t[0].to_string())
        })
        .collect()
}

/// `max |a − b| / max(max |b|, floor)`.
pub fn rel_err(a: &[f64], b: &[f64], floor: f64) -> f64 {
    assert_eq!(a.len(), b.len(), "length mismatch");
    let scale = b.iter().fold(floor, |m, x| m.max(x.abs()));
    a.iter()
        .zip(b)
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()))
        / scale
}

/// Every case of `MANIFEST.tsv` as `(name, task, is_ok)`.
pub fn manifest() -> Vec<(String, String, bool)> {
    let text = std::fs::read_to_string(dir().join("MANIFEST.tsv")).unwrap();
    text.lines()
        .filter_map(|l| {
            let t: Vec<&str> = l.split('\t').collect();
            (t.len() >= 3).then(|| (t[0].to_string(), t[1].to_string(), t[2] == "ok"))
        })
        .collect()
}
