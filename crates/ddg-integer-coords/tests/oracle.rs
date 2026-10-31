//! Differential tests against the committed fixtures of the Lean oracle
//! (`lean/oracle/fixtures/ic`, `ddg-oracle ic`; format in `lean/oracle/README.md`): the tasks
//! `corner`, `trace`, `tracecurve`, `identify`, `splitedge` and `flip` (Eq. 3, 4 and the flip
//! back), compared by exact integer equality, line by line. The fixtures are read, never regenerated.

use std::path::PathBuf;

use ddg_integer_coords::formulas::{
    Step, StepError, corner_coordinate, emanating, face_corner_numerator, face_emanating,
    face_valid, flip_doubled, flip_normal, identify_index, n_minus, roundabout_update,
    split_edge_normal, trace_step,
};
use ddg_integer_coords::{Crossing, IntegerError, NormalGluing, trace_with};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/ic")
}

/// One case: its name, task and lines (tokens), without `case`, `task`, `status` and `end`.
struct Case {
    name: String,
    task: String,
    status: Option<String>,
    lines: Vec<Vec<String>>,
}

fn read_cases(path: &PathBuf) -> Vec<Case> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut cases = Vec::new();
    let mut current: Option<Case> = None;
    for line in text.lines() {
        let tokens: Vec<String> = line.split_whitespace().map(str::to_string).collect();
        if tokens.is_empty() || tokens[0].starts_with('#') {
            continue;
        }
        match tokens[0].as_str() {
            "case" => {
                current = Some(Case {
                    name: tokens[1].clone(),
                    task: String::new(),
                    status: None,
                    lines: Vec::new(),
                })
            }
            "task" => current.as_mut().unwrap().task = tokens[1].clone(),
            "status" => current.as_mut().unwrap().status = Some(tokens[1..].join(" ")),
            "end" => cases.push(current.take().unwrap()),
            _ => current.as_mut().unwrap().lines.push(tokens),
        }
    }
    cases
}

fn int(t: &str) -> i64 {
    t.parse().unwrap_or_else(|_| panic!("not an integer: {t}"))
}

/// The pairs of `.in` / `.out` files of a task.
fn pairs(task: &str) -> Vec<(Case, Case)> {
    let mut out = Vec::new();
    let mut names: Vec<PathBuf> = std::fs::read_dir(fixtures())
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "in"))
        .collect();
    names.sort();
    for input in names {
        let output = input.with_extension("out");
        let ins = read_cases(&input);
        let outs = read_cases(&output);
        assert_eq!(ins.len(), outs.len(), "{}", input.display());
        for (i, o) in ins.into_iter().zip(outs) {
            assert_eq!(i.name, o.name);
            if i.task == task {
                assert_eq!(o.status.as_deref(), Some("ok"), "{}", i.name);
                out.push((i, o));
            }
        }
    }
    assert!(!out.is_empty(), "no fixture for task {task}");
    out
}

#[test]
fn corner() {
    let mut n = 0;
    for (input, output) in pairs("corner") {
        assert_eq!(input.lines.len(), output.lines.len());
        for (i, o) in input.lines.iter().zip(&output.lines) {
            let f = [int(&i[1]), int(&i[2]), int(&i[3])];
            let mut mine = vec![
                "face".to_string(),
                f[0].to_string(),
                f[1].to_string(),
                f[2].to_string(),
                (face_valid(f) as u8).to_string(),
            ];
            for k in 0..3 {
                mine.push(face_corner_numerator(f, k).to_string());
            }
            for k in 0..3 {
                mine.push(
                    corner_coordinate(f[(k + 1) % 3], f[k], f[(k + 2) % 3])
                        .map_or("odd".to_string(), |c| c.to_string()),
                );
            }
            for k in 0..3 {
                mine.push(face_emanating(f, k).to_string());
            }
            assert_eq!(&mine, o, "{}", input.name);
            n += 1;
        }
    }
    eprintln!("corner: {n} faces agree");
}

#[test]
fn trace() {
    let mut n = 0;
    for (input, output) in pairs("trace") {
        assert_eq!(input.lines.len(), output.lines.len());
        for (i, o) in input.lines.iter().zip(&output.lines) {
            let (nij, nik, nkj, p) = (int(&i[1]), int(&i[2]), int(&i[3]), int(&i[4]));
            let tail = match trace_step(nij, nik, nkj, p) {
                Ok(Step::ToIk(q)) => format!("step ik {q}"),
                Ok(Step::ToKj(q)) => format!("step kj {q}"),
                Ok(Step::Vertex) => "step vertex".to_string(),
                Err(StepError::InvalidFace) => "error face not valid".to_string(),
                Err(StepError::OutOfRange) => "error crossing out of range".to_string(),
            };
            let mine = format!("{} {tail}", i.join(" "));
            assert_eq!(mine, o.join(" "), "{}", input.name);
            n += 1;
        }
    }
    eprintln!("trace: {n} steps agree");
}

/// The hand-built Δ-complexes of the `tracecurve` fixtures.
struct Gluing {
    vert: Vec<usize>,
    twin: Vec<i64>,
    n: Vec<i64>,
}

impl NormalGluing for Gluing {
    fn n_halfedges(&self) -> usize {
        self.vert.len()
    }
    fn twin(&self, h: usize) -> Option<usize> {
        (self.twin[h] >= 0).then_some(self.twin[h] as usize)
    }
    fn tail(&self, h: usize) -> usize {
        self.vert[h]
    }
    fn normal(&self, h: usize) -> i64 {
        self.n[h]
    }
}

#[test]
fn tracecurve() {
    let mut n_curves = 0;
    let mut n_errors = 0;
    for (input, output) in pairs("tracecurve") {
        let mut g = Gluing {
            vert: Vec::new(),
            twin: Vec::new(),
            n: Vec::new(),
        };
        let mut starts = Vec::new();
        let mut bound = None;
        for l in &input.lines {
            match l[0].as_str() {
                "faces" => {
                    let nf = int(&l[1]) as usize;
                    g.vert = vec![0; 3 * nf];
                    g.twin = vec![-2; 3 * nf];
                    g.n = vec![0; 3 * nf];
                }
                "f" | "tw" | "nf" => {
                    let f = int(&l[1]) as usize;
                    for k in 0..3 {
                        let x = int(&l[2 + k]);
                        match l[0].as_str() {
                            "f" => g.vert[3 * f + k] = x as usize,
                            "tw" => g.twin[3 * f + k] = x,
                            _ => g.n[3 * f + k] = x,
                        }
                    }
                }
                "start" => starts.push((int(&l[1]) as usize, int(&l[2]))),
                "bound" => bound = Some(int(&l[1]) as u64),
                other => panic!("unexpected line {other}"),
            }
        }
        // The oracle's default bound, and its loop runs `bound + 1` times.
        let bound = bound.unwrap_or_else(|| {
            g.n.iter().map(|&x| x.max(0) as u64).sum::<u64>() + g.n.len() as u64
        }) + 1;
        let mut mine: Vec<String> = Vec::new();
        for (i, &(h, p)) in starts.iter().enumerate() {
            mine.push(format!("curve {i} {h} {p}"));
            let mut visited: Vec<Crossing> = Vec::new();
            let mut seen = std::collections::HashSet::new();
            let result = trace_with(&g, Crossing::new(h, p), Some(bound), &mut |z| {
                visited.push(z);
                let key = match z.reversed(&g) {
                    Some(r) if r.halfedge < z.halfedge => (r.halfedge, r.index),
                    _ => (z.halfedge, z.index),
                };
                if seen.insert(key) {
                    Ok(())
                } else {
                    Err(IntegerError::Conversion {
                        edge: 0,
                        reason: "closed".into(),
                    })
                }
            });
            for z in &visited {
                mine.push(format!("x {} {}", z.halfedge, z.index));
            }
            match result {
                Ok((end, steps)) => {
                    assert_eq!(steps as usize, visited.len());
                    mine.push(format!(
                        "stop vertex {} {} {} {}",
                        end.vertex,
                        end.corner / 3,
                        end.corner % 3,
                        steps
                    ));
                }
                Err(_) => {
                    n_errors += 1;
                    mine.push("stop error".to_string())
                }
            }
            n_curves += 1;
        }
        let theirs: Vec<String> = output
            .lines
            .iter()
            .map(|l| {
                if l[0] == "stop" && l[1] == "error" {
                    "stop error".to_string()
                } else {
                    l.join(" ")
                }
            })
            .collect();
        assert_eq!(mine, theirs, "{}", input.name);
    }
    eprintln!("tracecurve: {n_curves} curves agree ({n_errors} closed)");
    assert!(n_errors > 0);
}

#[test]
fn identify() {
    let (mut n, mut invalid) = (0, 0);
    for (input, output) in pairs("identify") {
        for (i, o) in input.lines.iter().zip(&output.lines) {
            let v: Vec<i64> = i[1..].iter().map(|t| int(t)).collect();
            let (r, nai, nij, nja, p, deg) = (v[0], v[1], v[2], v[3], v[4], v[5] as u32);
            let theirs = int(&o[8]);
            assert_eq!(o[7], "index");
            match identify_index(r, nai, nij, nja, p, deg) {
                Some(k) => {
                    assert_eq!(k, theirs, "{}: {}", input.name, i.join(" "));
                    n += 1;
                }
                None => {
                    // The face ija is not valid: our formula refuses to round Eq. 1 (I14).
                    assert!(!face_valid([nij, nja, nai]), "{}", i.join(" "));
                    invalid += 1;
                }
            }
        }
    }
    eprintln!("identify: {n} indices agree, {invalid} arguments with an invalid face refused");
}

#[test]
fn splitedge() {
    let mut n = 0;
    for (input, output) in pairs("splitedge") {
        for (i, o) in input.lines.iter().zip(&output.lines) {
            let mine = split_edge_normal(int(&i[1]), int(&i[2]));
            assert_eq!(mine, int(&o[3]), "{}", i.join(" "));
            n += 1;
        }
    }
    eprintln!("splitedge: {n} agree");
}

/// One block of the `flip` task: the quad, the optional `roundabout` arguments and `again`.
struct FlipBlock {
    quad: [i64; 5],
    roundabout: Option<[i64; 4]>,
    again: bool,
}

#[test]
fn flip() {
    let (mut n_quads, mut n_invalid, mut n_roundabouts, mut n_again) = (0, 0, 0, 0);
    for (input, output) in pairs("flip") {
        let mut blocks: Vec<FlipBlock> = Vec::new();
        for l in &input.lines {
            match l[0].as_str() {
                "quad" => blocks.push(FlipBlock {
                    quad: [1, 2, 3, 4, 5].map(|k| int(&l[k])),
                    roundabout: None,
                    again: false,
                }),
                "roundabout" => {
                    blocks.last_mut().unwrap().roundabout = Some([1, 2, 3, 4].map(|k| int(&l[k])))
                }
                "again" => blocks.last_mut().unwrap().again = true,
                x => panic!("{}: unexpected input line {x}", input.name),
            }
        }
        // The output lines of each block: from its `quad` echo to the next.
        let mut out: Vec<Vec<&Vec<String>>> = Vec::new();
        for l in &output.lines {
            if l[0] == "quad" {
                out.push(Vec::new());
            }
            out.last_mut().unwrap().push(l);
        }
        assert_eq!(blocks.len(), out.len(), "{}", input.name);
        for (b, lines) in blocks.iter().zip(&out) {
            let [nij, njk, nki, nil, nlj] = b.quad;
            let ctx = format!("{} quad {:?}", input.name, b.quad);
            assert_eq!(lines[0][1..].iter().map(|t| int(t)).collect::<Vec<_>>(), b.quad.to_vec());
            let v1 = face_valid([nij, njk, nki]);
            let v2 = face_valid([nij, nil, nlj]);
            let faces = lines[1];
            assert_eq!(faces[0], "faces", "{ctx}");
            assert_eq!((int(&faces[1]) == 1, int(&faces[2]) == 1), (v1, v2), "{ctx}");
            n_quads += 1;
            if faces.len() == 4 {
                assert_eq!(faces[3], "invalid", "{ctx}");
                assert!(!(v1 && v2));
                assert_eq!(flip_doubled(nij, njk, nki, nil, nlj), None, "{ctx}");
                assert!(flip_normal(nij, njk, nki, nil, nlj).is_err(), "{ctx}");
                assert_eq!(lines.len(), 2, "{ctx}");
                n_invalid += 1;
                continue;
            }
            let mut rest = lines[2..].iter();
            let l = rest.next().unwrap();
            assert_eq!(l[0], "flip", "{ctx}");
            let d = flip_doubled(nij, njk, nki, nil, nlj).unwrap();
            let nkl = flip_normal(nij, njk, nki, nil, nlj).unwrap();
            assert_eq!((d, nkl), (int(&l[1]) as i128, int(&l[2])), "{ctx}");
            let l = rest.next().unwrap();
            assert_eq!(l[0], "newfaces", "{ctx}");
            assert_eq!(
                (face_valid([nki, nil, nkl]), face_valid([nlj, njk, nkl])),
                (int(&l[1]) == 1, int(&l[2]) == 1),
                "{ctx}"
            );
            let l = rest.next().unwrap();
            assert_eq!(l[0], "emanating", "{ctx}");
            let (ek, el) = (emanating(nil, nki, nkl), emanating(njk, nlj, nkl));
            assert_eq!((ek, el), (int(&l[1]), int(&l[2])), "{ctx}");
            if let Some([rki, dk, rlj, dl]) = b.roundabout {
                let l = rest.next().unwrap();
                assert_eq!(l[0], "roundabout", "{ctx}");
                let rkl = roundabout_update(rki, ek, n_minus(nki), dk as u32).unwrap();
                let rlk = roundabout_update(rlj, el, n_minus(nlj), dl as u32).unwrap();
                assert_eq!((rkl, rlk), (int(&l[1]), int(&l[2])), "{ctx}");
                n_roundabouts += 1;
            }
            if b.again {
                let l = rest.next().unwrap();
                assert_eq!(l[0], "again", "{ctx}");
                let back = flip_normal(nkl, nki, nil, nlj, njk).unwrap();
                assert_eq!(back, int(&l[1]), "{ctx}");
                assert_eq!(back, nij, "{ctx}: ic-flip-involution");
                n_again += 1;
            }
            assert!(rest.next().is_none(), "{ctx}");
        }
    }
    eprintln!(
        "flip: {n_quads} quads ({n_invalid} with an invalid face), {n_roundabouts} roundabout \
         updates, {n_again} flips back: all agree"
    );
    assert!(n_quads >= 300 && n_roundabouts >= 300 && n_again >= 300);
}
