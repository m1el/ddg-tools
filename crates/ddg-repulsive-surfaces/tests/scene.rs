//! Scene reader tests (WP12): grammar, errors, round trips, every reference scene.

use std::path::{Path, PathBuf};

use ddg_repulsive_surfaces::Penalty;
use ddg_repulsive_surfaces::RsError;
use ddg_repulsive_surfaces::constraints::Constraint;
use ddg_repulsive_surfaces::scene::*;

fn p() -> &'static Path {
    Path::new("test.txt")
}

fn parse(t: &str) -> Result<SceneFile, RsError> {
    SceneFile::parse(t, p())
}

fn line_of(e: RsError) -> usize {
    match e {
        RsError::Parse { line, .. } => line,
        o => panic!("expected Parse, got {o}"),
    }
}

fn strip(s: &mut SceneFile) {
    for c in &mut s.constraints {
        c.line = 0;
    }
    for c in &mut s.potentials {
        c.line = 0;
    }
    for c in &mut s.obstacles {
        c.line = 0;
    }
    for c in &mut s.implicits {
        c.line = 0;
    }
}

#[test]
fn grammar_and_values() {
    let s = parse(
        "# comment\n\nrepel_mesh a.obj 6 12\nconstrain barycenter\nconstrain area 3 200 1.5\nconstrain volume\n\
         constrain vertices 1 2 3 move 0 0 .5 10\nconstrain boundary_vertices\nconstrain boundary_normals\n\
         minimize area 2\nminimize boundary_length 3 1.5\nminimize willmore\n\
         obstacle b.obj 0.1 centered\npoint_cloud_obstacle c.obj\n\
         implicit torus attract 2 1 1 .8 0 0 0\nimplicit plane repel 2 20 0 0 0 0 1 0\n\
         iteration_limit 0\ntime_limit 500\nallow_barycenter_shift\ndisable_near_field\nmethod l2\nlog out.csv\nautotarget_volume\n",
    )
    .unwrap();
    assert_eq!((s.alpha, s.beta), (6.0, 12.0));
    assert_eq!(s.constraints.len(), 6);
    assert_eq!(
        s.constraints[1].value,
        ConstrainStmt::Area(Some(Growth {
            multiplier: 3.0,
            iterations: 200,
            addition: 1.5
        }))
    );
    assert_eq!(s.iteration_limit, None);
    assert_eq!(s.time_limit_ms, Some(500));
    assert_eq!(s.potentials[2].value.weight, 1.0);
    assert_eq!(s.potentials[1].value.target, 1.5);
    assert!(s.obstacles[0].value.centered && s.obstacles[1].value.point_cloud);
    assert_eq!(s.autotarget_volume, Some(1.0));
    assert_eq!(s.method, SceneMethod::L2);
    let mut s1 = s.clone();
    let mut s2 = parse(&s.to_text()).unwrap();
    strip(&mut s1);
    strip(&mut s2);
    assert_eq!(s1, s2);
}

#[test]
fn errors_are_located() {
    let cases: [(&str, usize); 11] = [
        ("repel_mesh a.obj 6 12\nfrobnicate 3\n", 2),
        ("repel_mesh a.obj\n", 1),
        ("repel_mesh a.obj 6 12\nminimzie area\n", 2),
        ("repel_mesh a.obj 6 12\nconstrain area 3\n", 2),
        ("repel_mesh a.obj 6 12\nconstrain foo\n", 2),
        ("repel_mesh a.obj 6 12\nimplicit sphere repel 2 1 1\n", 2),
        ("repel_mesh a.obj 6 12\nmethod nope\n", 2),
        ("repel_mesh a.obj 1 2\n", 1),
        ("repel_mesh a.obj 6 12\nrepel_mesh b.obj 6 12\n", 2),
        ("constrain barycenter\n", 0),
        ("repel_mesh a.obj 6 12\nconstrain vertices 1 x\n", 2),
    ];
    for (text, line) in cases {
        assert_eq!(line_of(parse(text).unwrap_err()), line, "{text}");
    }
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("rs-scene-test-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn resolve_small_scene() {
    let d = scratch("resolve");
    // octahedron with UVs: u > 0 on the corners of the first face only (vertices 0 and 2... see below)
    let obj = "v 1 0 0\nv -1 0 0\nv 0 1 0\nv 0 -1 0\nv 0 0 1\nv 0 0 -1\n\
        vt 1 0\nvt 0 0\nf 1/1 3/2 5/2\nf 3 2 5\nf 2 4 5\nf 4 1 5\nf 3 1 6\nf 2 3 6\nf 4 2 6\nf 1 4 6\n";
    std::fs::write(d.join("m.obj"), obj).unwrap();
    std::fs::write(
        d.join("cube.obj"),
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\nf 1 3 2\nf 1 2 4\nf 1 4 3\nf 2 3 4\n",
    )
    .unwrap();
    let scene = "repel_mesh m.obj 6 12\nconstrain barycenter\nconstrain volume 2 10\nconstrain vertices 0 1\n\
        minimize area 3\nminimize area_deviation 100\nminimize squared_error\n\
        obstacle cube.obj 0.5 centered\nimplicit sphere attract 2 1 3 0 0 0\nimplicit cylinder repel 2 1 .1 0 0 0 0 0 1\n\
        autotarget_volume 2\nlog x.csv\n";
    std::fs::write(d.join("scene.txt"), scene).unwrap();
    let s = SceneSpec::load(d.join("scene.txt")).unwrap();
    assert_eq!(s.surface.n_faces(), 8);
    assert_eq!(s.constraints.len(), 3);
    assert!(matches!(s.constraints[2].constraint, Constraint::Pins { .. }));
    assert_eq!(s.penalties.len(), 6);
    match &s.penalties[0] {
        Penalty::TotalArea { weight, power } => assert_eq!((*weight, *power), (3.0, 2.0)),
        o => panic!("{o:?}"),
    }
    match &s.penalties[1] {
        Penalty::AreaDeviation { target, .. } => {
            assert!((*target - 4.0 * 3f64.sqrt()).abs() < 1e-12)
        }
        o => panic!("{o:?}"),
    }
    match &s.penalties[4] {
        Penalty::ImplicitAttractor { vertices, .. } => {
            assert_eq!(vertices.as_deref(), Some(&[0u32][..]))
        }
        o => panic!("{o:?}"),
    }
    assert!((s.obstacle_volume.abs() - 1.0 / 6.0).abs() < 1e-12);
    assert_eq!(s.log.as_ref().unwrap(), &d.join("x.csv"));
    std::fs::write(d.join("bad.txt"), "repel_mesh m.obj 6 12\nconstrain vertices 99\n").unwrap();
    assert_eq!(line_of(SceneSpec::load(d.join("bad.txt")).unwrap_err()), 2);
    std::fs::write(
        d.join("bad2.txt"),
        "repel_mesh m.obj 6 12\nconstrain boundary_vertices\n",
    )
    .unwrap();
    assert_eq!(line_of(SceneSpec::load(d.join("bad2.txt")).unwrap_err()), 2);
    std::fs::write(d.join("bad3.txt"), "repel_mesh nothere.obj 6 12\n").unwrap();
    match SceneSpec::load(d.join("bad3.txt")).unwrap_err() {
        RsError::Io { path, .. } => assert!(path.ends_with("nothere.obj")),
        o => panic!("{o}"),
    }
    std::fs::remove_dir_all(d).ok();
}

fn scenes_dir() -> Option<PathBuf> {
    let here = Path::new(env!("CARGO_MANIFEST_DIR"));
    let cands = [
        std::env::var("RS_SCENES").ok().map(PathBuf::from),
        Some(here.join("../../data/repulsive-surfaces/scenes")),
        Some(PathBuf::from(
            "/var/data/discrete-geometry/RepulsiveSurfaces/repulsive-surfaces/scenes",
        )),
    ];
    cands.into_iter().flatten().find(|d| d.is_dir())
}

/// Every scene of the reference parses; only `Boundary/scene-tall.txt` (no `α β`) is an error;
/// meshes that do not exist are listed; scenes whose files exist resolve.
#[test]
fn reference_scenes_parse() {
    let Some(dir) = scenes_dir() else {
        eprintln!("no scenes directory; run scripts/rs-fetch-scenes.sh");
        return;
    };
    let mut parsed = 0;
    let mut missing = Vec::new();
    let mut resolved = 0;
    let mut errors = Vec::new();
    for d in std::fs::read_dir(&dir).unwrap().flatten() {
        if !d.path().is_dir() {
            continue;
        }
        for f in std::fs::read_dir(d.path()).unwrap().flatten() {
            let name = f.file_name().to_string_lossy().to_string();
            if !name.ends_with(".txt") {
                continue;
            }
            let label = format!("{}/{name}", d.file_name().to_string_lossy());
            let sf = match SceneFile::read(f.path()) {
                Ok(s) => s,
                Err(e) => {
                    errors.push((label, e.to_string()));
                    continue;
                }
            };
            parsed += 1;
            let again = SceneFile::parse(&sf.to_text(), &sf.path).unwrap();
            let (mut a, mut b) = (sf.clone(), again);
            strip(&mut a);
            strip(&mut b);
            a.path = PathBuf::new();
            b.path = PathBuf::new();
            assert_eq!(a, b, "{label}");
            let gone: Vec<_> = sf.mesh_files().into_iter().filter(|m| !m.exists()).collect();
            if !gone.is_empty() {
                for g in gone {
                    missing.push(format!(
                        "{label}: {}",
                        g.file_name().unwrap().to_string_lossy()
                    ));
                }
                continue;
            }
            match SceneSpec::resolve(sf) {
                Ok(_) => resolved += 1,
                Err(e) => errors.push((label, e.to_string())),
            }
        }
    }
    missing.sort();
    eprintln!(
        "parsed {parsed}, resolved {resolved}\nmissing meshes:\n  {}\nerrors: {errors:?}",
        missing.join("\n  ")
    );
    assert_eq!(errors.len(), 1, "{errors:?}");
    assert!(
        errors[0].0 == "Boundary/scene-tall.txt" && errors[0].1.contains("exponents"),
        "{errors:?}"
    );
    assert!(parsed >= 34);
}
