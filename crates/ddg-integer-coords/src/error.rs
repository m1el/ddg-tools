//! Errors of this crate.

use std::fmt;

use ddg_core::GeometryError;
use ddg_core::mollify::MollifyError;
use ddg_intrinsic::IntrinsicError;

/// Which part of the validity of a state (SPEC §2.4) a check belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ValidityItem {
    /// The ranges of the stored quantities (SPEC §2.2): `n ≥ −1`, `n = −1` on boundary edges,
    /// a roundabout in `[0, deg₀)` exactly for halfedges whose tail is an input vertex, the
    /// place `q⁰` of every vertex, array lengths.
    Ranges,
    /// Item 1: the numerators of Eq. 1 are even in every face.
    FaceParity,
    /// Item 3: an edge with `n = −1` runs along an input edge, and the roundabouts of its
    /// halfedges leaving input vertices name that input edge.
    AlongInputEdge,
    /// Item 4: every input edge is extracted (Supp. Alg. 4) as a curve from its tail to its head
    /// through inserted vertices on it only, with matching roundabout labels at both ends, and
    /// every crossing belongs to exactly one such curve. This includes the agreement of the
    /// roundabouts around every input vertex with the emanating counts (the window of Alg. 4).
    Extraction,
}

impl fmt::Display for ValidityItem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            ValidityItem::Ranges => "ranges of the stored quantities",
            ValidityItem::FaceParity => "item 1 (Eq. 1 integral)",
            ValidityItem::AlongInputEdge => "item 3 (edges along input edges)",
            ValidityItem::Extraction => "item 4 (extraction of the input edges)",
        })
    }
}

/// Where invalid integer data was found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Place {
    /// A face of the triangulation.
    Face(u32),
    /// An edge of the triangulation.
    Edge(u32),
    /// A halfedge of the triangulation (its roundabout, or its corner).
    Halfedge(u32),
    /// A vertex of the triangulation.
    Vertex(u32),
    /// A combinatorial crossing `(halfedge, index)`: the curve leaves the face of `halfedge`
    /// through it at the `index`-th crossing from its tail.
    Crossing {
        /// The halfedge.
        halfedge: u32,
        /// The crossing index from the halfedge's tail.
        index: i64,
        /// The face of the halfedge.
        face: u32,
    },
    /// An edge of the input mesh.
    InputEdge(u32),
}

impl fmt::Display for Place {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Place::Face(x) => write!(f, "face {x}"),
            Place::Edge(x) => write!(f, "edge {x}"),
            Place::Halfedge(x) => write!(f, "halfedge {x}"),
            Place::Vertex(x) => write!(f, "vertex {x}"),
            Place::Crossing {
                halfedge,
                index,
                face,
            } => write!(f, "crossing ({halfedge}, {index}) of face {face}"),
            Place::InputEdge(x) => write!(f, "input edge {x}"),
        }
    }
}

/// Integer data that is not that of any realization (SPEC §2.4, §3.1): what is wrong, which
/// item of the validity it violates and where.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvalidCoordinates {
    /// The item of SPEC §2.4.
    pub item: ValidityItem,
    /// Where.
    pub place: Place,
    /// What, in words, with the values involved.
    pub detail: String,
}

impl fmt::Display for InvalidCoordinates {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid integer coordinates at {} ({}): {}",
            self.place, self.item, self.detail
        )
    }
}

/// The errors of this crate. Nothing is changed when an operation returns one.
#[derive(Debug, Clone, PartialEq)]
pub enum IntegerError {
    /// An error of the shared intrinsic-triangulation code (connectivity, input mesh, lengths,
    /// refusals of flips and splits, traces of the signpost representation).
    Intrinsic(IntrinsicError),
    /// The geometry was rejected or not built for the mesh.
    Geometry(GeometryError),
    /// Mollification refused its input.
    Mollify(MollifyError),
    /// An element id is out of range.
    IndexOutOfRange {
        /// What kind of element.
        what: &'static str,
        /// The id.
        index: usize,
        /// The exclusive upper bound it violates.
        bound: usize,
    },
    /// An integer computation would leave the range of `i64` (SPEC I2); nothing wraps.
    Overflow {
        /// The edge whose normal coordinate was being computed or read.
        edge: u32,
    },
    /// The certificate marks every crossing in a bitmap of `Σ n⁺` bits; this state has more
    /// crossings than [`MAX_CERTIFIED_CROSSINGS`](crate::certify::MAX_CERTIFIED_CROSSINGS).
    TooManyCrossings {
        /// `Σ n⁺`.
        total: u128,
        /// The limit.
        limit: u64,
    },
    /// The integers are not those of any realization; see [`InvalidCoordinates`].
    InvalidCoordinates(Box<InvalidCoordinates>),
    /// An inserted vertex could not be removed (SPEC §4.6, I8): no edge at it can be flipped
    /// to lower its degree, and its degree is not three (or its star cannot be merged). The
    /// state is restored to what it was before the removal.
    RemovalStuck {
        /// The vertex.
        vertex: u32,
        /// Its degree when no flip was left.
        degree: usize,
    },
    /// A face split was refused because the point is closer to an arc of the face (an input
    /// edge crossing it) than [`FaceSplitOptions::too_close`](crate::insert::FaceSplitOptions)
    /// allows (SPEC §4.4, I6). Nothing is changed.
    TooCloseToCurve {
        /// The face of the triangulation.
        face: u32,
        /// The distance to the nearest arc, in units of the face's shortest side.
        distance: f64,
        /// The option.
        limit: f64,
    },
    /// `RecoverBarycentric` (App. A, SPEC I16) gave no point of the input face of the region:
    /// the corners of the region do not span it, a corner is not on that input face, or the
    /// result lies outside it by more than rounding. Nothing is changed.
    RecoverFailed {
        /// The face of the triangulation.
        face: u32,
        /// The input face of the region.
        input_face: u32,
        /// What went wrong, with the values.
        detail: String,
    },
    /// An operation this representation does not provide (yet), e.g. the removal of an
    /// inserted boundary vertex (SPEC I7). Nothing is changed.
    Unsupported {
        /// The operation.
        operation: &'static str,
        /// Why, or what is missing.
        reason: String,
    },
    /// [`IntegerTriangulation::from_intrinsic`](crate::IntegerTriangulation::from_intrinsic)
    /// could not read integers off the signpost representation (floating point decides it).
    Conversion {
        /// The edge of the signpost triangulation.
        edge: u32,
        /// What went wrong.
        reason: String,
    },
}

impl IntegerError {
    /// An [`IntegerError::InvalidCoordinates`].
    pub(crate) fn invalid(item: ValidityItem, place: Place, detail: String) -> IntegerError {
        IntegerError::InvalidCoordinates(Box::new(InvalidCoordinates {
            item,
            place,
            detail,
        }))
    }

    /// The validity item an [`IntegerError::InvalidCoordinates`] names.
    pub fn validity_item(&self) -> Option<ValidityItem> {
        match self {
            IntegerError::InvalidCoordinates(x) => Some(x.item),
            _ => None,
        }
    }
}

impl fmt::Display for IntegerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IntegerError::Intrinsic(e) => write!(f, "{e}"),
            IntegerError::Geometry(e) => write!(f, "{e}"),
            IntegerError::Mollify(e) => write!(f, "mollification: {e}"),
            IntegerError::IndexOutOfRange { what, index, bound } => {
                write!(f, "{what} {index} out of range (there are {bound})")
            }
            IntegerError::Overflow { edge } => write!(
                f,
                "an integer computation at edge {edge} would leave the range of i64"
            ),
            IntegerError::TooManyCrossings { total, limit } => write!(
                f,
                "the state has {total} crossings, more than the certificate handles ({limit})"
            ),
            IntegerError::InvalidCoordinates(x) => write!(f, "{x}"),
            IntegerError::RemovalStuck { vertex, degree } => write!(
                f,
                "vertex {vertex} could not be removed: no edge at it can be flipped to lower its \
                 degree {degree}"
            ),
            IntegerError::TooCloseToCurve {
                face,
                distance,
                limit,
            } => write!(
                f,
                "the point is {distance:e} (in units of the shortest side) from an input edge \
                 crossing face {face}, closer than the limit {limit:e}"
            ),
            IntegerError::RecoverFailed {
                face,
                input_face,
                detail,
            } => write!(
                f,
                "the place on input face {input_face} of a point of face {face} could not be \
                 recovered: {detail}"
            ),
            IntegerError::Unsupported { operation, reason } => {
                write!(f, "{operation} is not supported by integer coordinates: {reason}")
            }
            IntegerError::Conversion { edge, reason } => write!(
                f,
                "edge {edge} of the signpost triangulation could not be converted: {reason}"
            ),
        }
    }
}

impl std::error::Error for IntegerError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            IntegerError::Intrinsic(e) => Some(e),
            IntegerError::Geometry(e) => Some(e),
            IntegerError::Mollify(e) => Some(e),
            _ => None,
        }
    }
}

impl From<IntrinsicError> for IntegerError {
    fn from(e: IntrinsicError) -> Self {
        IntegerError::Intrinsic(e)
    }
}

impl From<GeometryError> for IntegerError {
    fn from(e: GeometryError) -> Self {
        IntegerError::Geometry(e)
    }
}

impl From<MollifyError> for IntegerError {
    fn from(e: MollifyError) -> Self {
        IntegerError::Mollify(e)
    }
}
