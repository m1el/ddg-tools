//! Bad data is an error naming the place, never a panic and never a huge allocation (SPEC §3.1,
//! I2): extreme normal coordinates, out-of-range ids, hand-built segments.

use ddg_core::Geometry;
use ddg_core::generators::icosphere;
use ddg_integer_coords::{
    Crossing, CurveSegment, IntegerError, IntegerTriangulation, NO_ROUNDABOUT, Parts, Place,
    ValidityItem,
};
use ddg_intrinsic::SurfacePoint;

fn sphere() -> IntegerTriangulation {
    let g = icosphere(1).unwrap();
    let mesh = ddg_core::TriMesh::new(g.positions.len(), g.faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &g.positions).unwrap();
    IntegerTriangulation::new(&mesh, &geometry).unwrap()
}

fn with(parts: Parts) -> Result<ddg_integer_coords::CertifyReport, IntegerError> {
    IntegerTriangulation::from_parts(parts)?.certify()
}

#[test]
fn extreme_normal_coordinates_are_errors() {
    let t = sphere();
    // One edge at i64::MAX.
    let mut p = t.clone().into_parts();
    p.normal[0] = i64::MAX;
    let e = with(p).unwrap_err();
    assert!(e.validity_item().is_some(), "{e}");
    // Every edge at an even huge value: every face is valid (all corners cut, none emanates)
    // and nothing runs along an input edge; the roundabouts around a vertex then disagree with
    // the counts, which is found before anything of the size of Σ n⁺ is allocated.
    let mut p = t.clone().into_parts();
    p.normal.iter_mut().for_each(|n| *n = 1 << 40);
    let e = with(p).unwrap_err();
    assert_eq!(e.validity_item(), Some(ValidityItem::Extraction), "{e}");
    // Below −1.
    let mut p = t.clone().into_parts();
    p.normal[3] = -2;
    let e = with(p).unwrap_err();
    assert_eq!(e.validity_item(), Some(ValidityItem::Ranges));
    assert!(e.to_string().contains("edge 3"), "{e}");
    // A huge emanating count at a vertex: reported, not enumerated.
    let mut p = t.clone().into_parts();
    let c = p.complex.clone();
    let h = 0usize;
    p.normal[c.edge(c.next(h))] = i64::MAX - 1;
    let m = IntegerTriangulation::from_parts(p).unwrap();
    assert!(m.curve_ends(c.from_vertex(h)).is_err());
    assert!(m.input_directions(c.from_vertex(h)).is_err());
    assert!(m.certify().is_err());
}

#[test]
fn ranges_and_ids_are_checked() {
    let t = sphere();
    let mut p = t.clone().into_parts();
    p.roundabouts[5] = NO_ROUNDABOUT;
    let e = with(p).unwrap_err();
    assert_eq!(e.validity_item(), Some(ValidityItem::Ranges));
    match &e {
        IntegerError::InvalidCoordinates(x) => assert_eq!(x.place, Place::Halfedge(5)),
        other => panic!("{other:?}"),
    }
    let mut p = t.clone().into_parts();
    p.locations[2] = SurfacePoint::Edge { edge: 0, t: 0.5 };
    assert_eq!(
        with(p).unwrap_err().validity_item(),
        Some(ValidityItem::Ranges)
    );
    let mut p = t.clone().into_parts();
    p.normal.pop();
    assert!(IntegerTriangulation::from_parts(p).is_err());

    let nh = t.complex().n_halfedges();
    assert!(matches!(
        t.trace_from(Crossing::new(nh + 7, 0)),
        Err(IntegerError::IndexOutOfRange { .. })
    ));
    // Every edge runs along an input edge: there is no crossing 0 anywhere.
    let e = t.trace_from(Crossing::new(0, 0)).unwrap_err();
    assert_eq!(e.validity_item(), Some(ValidityItem::Extraction));
    assert!(t.extract_curve(Crossing::new(0, -3)).is_err());
    assert!(matches!(
        t.extract_edge(usize::MAX),
        Err(IntegerError::IndexOutOfRange { .. })
    ));
    assert!(t.input_directions(t.complex().n_vertices()).is_err());
    assert!(t.curve_ends(usize::MAX).is_err());
    let empty = CurveSegment {
        start: 0,
        start_vertex: 0,
        crossings: Vec::new(),
        end: ddg_integer_coords::Arrival {
            corner: 1,
            vertex: 1,
            back: Crossing::new(0, 0),
        },
    };
    assert!(t.identify_traced_edge(&empty).is_err());
    let wrong = CurveSegment {
        crossings: vec![Crossing::new(nh, 0)],
        ..empty
    };
    assert!(t.identify_traced_edge(&wrong).is_err());
}
