"""Type stubs of the extension module ``ddg_tools`` (built from ``crates/ddg-py``).

Array arguments accept any C-contiguous object with the buffer protocol; see README.md for the
element types. Results are ``array.array`` objects.
"""

from array import array
from collections.abc import Buffer, Callable, Sequence
from typing import Final, Literal, TypeAlias, TypedDict, final

__version__: Final[str]

PIPELINE: Final[str]
"""Name of the unwrap pipeline in this build, ``"ddg-uv"``."""

IS_PLACEHOLDER: Final[bool]
"""False: the pipeline is the real unwrap (``ddg-uv``). The flag exists so that a build that
links a stand-in pipeline would say so; no build does."""

class DdgError(Exception):
    """Base class of the errors raised by ddg_tools."""

class InputError(DdgError):
    """An argument is malformed. The message starts with the parameter name."""

class MeshError(DdgError):
    """The mesh cannot be unwrapped as it is (non-manifold, inconsistently oriented, ...)."""

class SolveError(DdgError):
    """The mesh was accepted but no valid UV map could be produced (the cone solver or the
    flattening failed, or a result did not pass the pipeline's own checks)."""

@final
class UnwrapResult:
    @property
    def uvs(self) -> array[float]:
        """``array('f')``, 2 per polygon corner in corner ("loop") order. Not packed: the
        charts are laid out where the flattening put them and overlap each other until
        ``pack_charts``."""
    @property
    def polygon_charts(self) -> array[int]:
        """``array('i')``, the chart of every polygon."""
    @property
    def cone_vertices(self) -> array[int]:
        """``array('i')``, the vertices that carry a cone, ascending."""
    @property
    def cone_angles(self) -> array[float]:
        """``array('d')``, cone angles in radians, parallel to ``cone_vertices``."""
    @property
    def cut_edges(self) -> array[int]:
        """``array('i')``, 2 vertex indices per edge the surface was cut along (seams
        included), sorted, ``min, max`` within an edge."""
    @property
    def report(self) -> dict[str, bool | int | float | str | list[str]]:
        """The pipeline's report: counts, cone statistics, ``area_distortion``, timings
        (``seconds_*``), the distortion measured by the binding (``measured_*``) and
        ``warnings`` (``list[str]``). Keys are listed in the crate README."""

def unwrap(
    positions: Buffer,
    corner_vertices: Buffer,
    polygon_starts: Buffer,
    triangles: Buffer,
    *,
    seam_edges: Buffer | None = None,
    pinned_vertices: Buffer | None = None,
    prescribed_vertices: Buffer | None = None,
    prescribed_angles: Buffer | None = None,
    forbidden_vertices: Buffer | None = None,
    lambda_: float = 0.1,
    positive_cones_only: bool = True,
    optimize_boundary: bool = False,
    split_overlaps: bool = False,
    split_max: int = 32,
    split_threshold: float = 0.001,
    cone_feedback: bool = True,
    feedback_threshold: float = 20.0,
    feedback_max_rounds: int = 4,
) -> UnwrapResult: ...
def pack_charts(
    uvs: Buffer,
    polygon_starts: Buffer,
    polygon_charts: Buffer,
    *,
    margin: float = 0.0,
) -> tuple[array[float], dict[str, int | float]]: ...

@final
class ReoptimizeResult:
    @property
    def cone_vertices(self) -> array[int]:
        """``array('i')``, the cone vertices as given."""
    @property
    def cone_angles(self) -> array[float]:
        """``array('d')``, the optimal cone angles in radians, parallel to ``cone_vertices``.
        An angle may be 0."""
    @property
    def log_scale(self) -> array[float]:
        """``array('d')``, the log conformal scale factor ``u``, one per vertex."""
    @property
    def boundary_curvature(self) -> array[float]:
        """``array('d')``, the boundary curvature of the flattened boundary, one per vertex
        (0 away from an optimized boundary)."""
    @property
    def report(self) -> dict[str, bool | int | float | str | None]:
        """``status`` (``"converged"``, ``"not_converged"``, ``"linear_solver_failed"``),
        ``method`` (``"direct_kkt"`` or ``"solve"``), ``distortion``, ``objective``,
        ``total_cone_angle``, ``n_components``, ``certificate_holds`` (at tolerance 1e-9),
        ``state_residual``, ``duality_gap``, ``seconds``, ``message`` (``str`` or ``None``) and
        ``seconds_total`` (including decoding and encoding)."""

def reoptimize_cones(
    positions: Buffer,
    triangles: Buffer,
    cone_vertices: Buffer,
    *,
    boundary: Literal["fixed", "free", "convex", "penalized"] = "fixed",
    distortion_weights: Buffer | None = None,
    lower_bounds: Buffer | None = None,
    upper_bounds: Buffer | None = None,
    prescribed_vertices: Buffer | None = None,
    prescribed_angles: Buffer | None = None,
    method: Literal["auto", "direct", "solve"] = "auto",
    lambda_: float = 0.1,
) -> ReoptimizeResult:
    """The least-distortion cone angles for cones at the given vertices.

    ``positions`` has 3 reals per vertex, ``triangles`` 3 vertex indices per triangle (a
    manifold, consistently oriented mesh in which every vertex is used; several components are
    allowed). Cones are placed at ``cone_vertices`` (interior vertices, each once) and nowhere
    else; their angles are optimized. ``boundary`` says what happens at the boundary: ``"fixed"``
    keeps ``u = 0`` there. ``distortion_weights``, ``lower_bounds`` (<= 0) and ``upper_bounds``
    (>= 0) are per vertex; with ``"penalized"``, ``lambda_`` (finite, >= 0) is the cost of the
    boundary curvature. ``prescribed_vertices`` and ``prescribed_angles`` fix some cone
    angles. Raises ``InputError`` (bad arguments, a cone vertex out of range, repeated, on the
    boundary or prescribed, a closed component whose cones cannot carry 2*pi*chi) and
    ``MeshError`` (not a valid manifold triangle mesh). A result that did not converge is not an
    error: see ``report["status"]``."""


@final
class SphereResult:
    @property
    def positions(self) -> array[float]:
        """``array('d')``, 3 per vertex: a point of the unit sphere."""
    @property
    def report(self) -> dict[str, bool | int | float | str | None]:
        """``pole``, ``is_injective``, ``flipped_triangles``, ``degenerate_triangles``,
        ``covering_degree``, ``quasi_conformal_mean``, ``quasi_conformal_max``,
        ``area_distortion``, ``log_scale_min``, ``log_scale_max``, ``max_radius_error``,
        ``centering_status`` (``"converged"``, ``"max_iterations"``, ``"stalled"``,
        ``"singular_jacobian"`` or ``None`` without centering), ``centering_iterations``,
        ``center_of_mass`` (``float`` or ``None``), ``seconds`` and ``seconds_total``."""

def map_to_sphere(
    positions: Buffer,
    triangles: Buffer,
    *,
    pole: int | None = None,
    center: bool = True,
) -> SphereResult:
    """A conformal map of a connected closed genus-0 triangle mesh to the unit sphere
    (Boundary First Flattening, Sec. 7.5: remove one vertex star, flatten the rest to a round
    disk, project stereographically, Moebius-center). ``pole`` is the removed vertex (``None``:
    the vertex with the most regular star); ``center=False`` leaves it at ``(0, 0, -1)``. The
    result is determined up to a rotation of the sphere. Raises ``InputError`` (bad arguments, a
    pole that is not a vertex) and ``MeshError`` (not a valid manifold triangle mesh, or not
    connected, closed and of genus 0)."""

# ----------------------------------------------------------------------------------------------
# Developability of triangle meshes (Stein, Grinspun, Crane 2018).
# ----------------------------------------------------------------------------------------------

ENERGIES: Final[list[str]]
"""The six energies: ``combinatorial``, ``combinatorial_max``, ``covariance``, ``covariance_max``,
``intrinsic``, ``intrinsic_max``."""

DEFAULT_DEFECT_TOLERANCE: Final[float]
"""Default ``defect_tolerance`` of the cut functions, radians."""

Energy = Literal["combinatorial", "combinatorial_max", "covariance", "covariance_max", "intrinsic", "intrinsic_max"]

@final
class DevelopabilityResult:
    @property
    def vertex_energy(self) -> array[float]:
        """``array('d')``, one per vertex."""
    @property
    def counted(self) -> array[int]:
        """``array('B')``: 1 where the vertex counts toward the total."""
    @property
    def energy(self) -> float: ...
    @property
    def report(self) -> dict[str, int | float]:
        """``near_ties``, ``n_counted``, ``max_vertex_energy``, ``seconds``, ``seconds_total``."""

def developability(
    positions: Buffer,
    corner_vertices: Buffer,
    polygon_starts: Buffer,
    triangles: Buffer,
    *,
    energy: Energy = "covariance",
    include_boundary: bool = True,
    omit_valence_three: bool = True,
    threads: int | None = None,
) -> DevelopabilityResult:
    """The developability energy per vertex of a polygon mesh with Blender's triangulation.
    Raises ``InputError`` and ``MeshError``."""

@final
class MakeDevelopableResult:
    @property
    def positions(self) -> array[float]:
        """``array('d')``, 3 per vertex: of the input mesh if ``faces`` is ``None``, else of the
        new mesh."""
    @property
    def faces(self) -> array[int] | None:
        """``array('i')``, 3 vertices per triangle; ``None`` if the topology did not change."""
    @property
    def face_polygon(self) -> array[int] | None:
        """``array('i')``: the polygon of the input each triangle of ``faces`` descends from
        (exact: tracked through subdivision, flips and collapses)."""
    @property
    def vertex_source(self) -> array[int]:
        """``array('i')``: the input vertex each vertex descends from, -1 if created by
        subdivision."""
    @property
    def report(self) -> dict[str, object]:
        """See the README: ``energy_initial``, ``energy_final``, ``cancelled``,
        ``topology_changed``, ``nonplanarity_*``, ``min_angle``, ``bad_faces``,
        ``inadmissible_trials``, ``stopped`` (``None`` or a dict ``round``, ``after_subdivision``,
        ``vertex``, ``reason``: the flow ended before that round because the energy is undefined
        there), ``rounds`` (list of dicts), ``warnings``, ..."""

def make_developable(
    positions: Buffer,
    corner_vertices: Buffer,
    polygon_starts: Buffer,
    triangles: Buffer,
    *,
    energy: Energy = "covariance",
    include_boundary: bool = True,
    omit_valence_three: bool = True,
    method: Literal["lbfgs", "gradient_descent", "fixed_step"] = "lbfgs",
    memory: int = 10,
    step: float = 1e-3,
    line_search: Literal["wolfe", "backtracking"] = "wolfe",
    first_trial: Literal["unit", "grow"] = "unit",
    first_trial_growth: float = 2.0,
    first_trial_max: float = 100.0,
    iterations: int = 500,
    subdivisions: int = 2,
    coarse_energy: Energy | None = "combinatorial",
    max_displacement: float = 0.5,
    gradient_tolerance: float | None = None,
    energy_tolerance: float | None = None,
    stall_iterations: int | None = None,
    initial_step: float | None = None,
    min_face_angle: float | None = None,
    pinned_vertices: Buffer | None = None,
    remesh: bool = True,
    small_angle: float = 0.017453292519943295,
    collapse_boundary_edges: bool = True,
    fallback_operations: bool = True,
    return_faces: bool = False,
    threads: int | None = None,
    progress: "Callable[[int, int, float, float], object] | None" = None,
    progress_interval: float = 0.1,
) -> MakeDevelopableResult:
    """The developability flow. ``progress(round, iteration, energy, gradient_norm)`` may return
    ``False`` to stop. Raises ``InputError`` and ``MeshError``, ``SolveError``."""

@final
class DevelopableSeamsResult:
    @property
    def seam_edges(self) -> array[int]: ...
    @property
    def terminals(self) -> array[int]: ...
    @property
    def vertex_energy(self) -> array[float]: ...
    @property
    def vertex_defect(self) -> array[float]: ...
    @property
    def star_embeddedness(self) -> array[float]: ...
    @property
    def report(self) -> dict[str, object]: ...

@final
class DevelopableUnwrapResult:
    @property
    def corner_uvs(self) -> array[float]:
        """``array('f')``, 2 per polygon corner, packed."""
    @property
    def seam_edges(self) -> array[int]: ...
    @property
    def polygon_chart(self) -> array[int]: ...
    @property
    def corner_log_scale(self) -> array[float]: ...
    @property
    def terminals(self) -> array[int]: ...
    @property
    def vertex_energy(self) -> array[float]: ...
    @property
    def vertex_defect(self) -> array[float]: ...
    @property
    def star_embeddedness(self) -> array[float]: ...
    @property
    def report(self) -> dict[str, object]: ...

def developable_seams(
    positions: Buffer,
    corner_vertices: Buffer,
    polygon_starts: Buffer,
    triangles: Buffer,
    *,
    seam_edges: Buffer | None = None,
    extra_terminals: Buffer | None = None,
    energy: Energy | None = None,
    epsilon: float = 1e-3,
    defect_tolerance: float | None = DEFAULT_DEFECT_TOLERANCE,
    length_weight: float = 1.0,
    alignment_weight: float = 0.0,
    extension: Literal["harmonic", "conjugate"] = "harmonic",
    overlap_samples_per_triangle: float | None = 0.0,
    threads: int | None = None,
) -> DevelopableSeamsResult:
    """Where to cut a polygon mesh (polygon edges only). Feed ``seam_edges`` to ``unwrap``."""

def developable_unwrap(
    positions: Buffer,
    corner_vertices: Buffer,
    polygon_starts: Buffer,
    triangles: Buffer,
    *,
    seam_edges: Buffer | None = None,
    extra_terminals: Buffer | None = None,
    energy: Energy | None = None,
    epsilon: float = 1e-3,
    defect_tolerance: float | None = DEFAULT_DEFECT_TOLERANCE,
    length_weight: float = 1.0,
    alignment_weight: float = 0.0,
    extension: Literal["harmonic", "conjugate"] = "harmonic",
    overlap_samples_per_triangle: float | None = None,
    margin: float = 1 / 256,
    threads: int | None = None,
) -> DevelopableUnwrapResult:
    """The paper's cut and flatten: cut through the high-energy vertices and flatten each chart
    with BFF keeping the boundary lengths; charts packed. A chart may overlap itself."""

# ----------------------------------------------------------------------------------------------
# The cut flow: "Variational Surface Cutting" (Sharp, Crane 2018) on a polygon mesh.
# ----------------------------------------------------------------------------------------------

class StaleStateError(DdgError):
    """A stored state (``CutSession.snapshot()``) belongs to another triangulation: the mesh or
    its seams changed since it was made. Start again."""

class NotAvailableError(DdgError):
    """The request is valid but not implemented yet: a cut session with more than two regions."""

@final
class CutDescription:
    """The cut on the caller's polygons: points on polygon edges and, per polygon, the paths of
    the cut through it. All arrays are flat; see README.md."""

    @property
    def n_edge_points(self) -> int: ...
    @property
    def n_paths(self) -> int: ...
    @property
    def edge_vertices(self) -> array[int]:
        """``array('i')``, 2 per edge point: the vertices ``[min, max]`` of its polygon edge."""
    @property
    def edge_t(self) -> array[float]:
        """``array('d')``: the point is ``(1 - t) x[min] + t x[max]``, ``0 < t < 1``."""
    @property
    def edge_positions(self) -> array[float]:
        """``array('d')``, 3 per edge point."""
    @property
    def edge_side(self) -> array[int]:
        """``array('i')``: on a seam, the polygon whose side of the seam the point is on, else -1."""
    @property
    def path_polygon(self) -> array[int]:
        """``array('i')``: the polygon of every path, ascending."""
    @property
    def path_start(self) -> array[int]:
        """``array('i')``: an edge point index, or ``-(j + 1)`` for junction ``j`` (three paths
        end at every junction)."""
    @property
    def path_end(self) -> array[int]:
        """``array('i')``: like ``path_start``."""
    @property
    def path_regions(self) -> array[int]:
        """``array('i')``, 2 per path: the regions on its left and on its right."""
    @property
    def path_interior_starts(self) -> array[int]:
        """``array('i')``, ``n_paths + 1`` offsets into the interior points."""
    @property
    def interior_corners(self) -> array[int]:
        """``array('i')``, 2 per interior point: the polygon corners of the diagonal it is on."""
    @property
    def interior_t(self) -> array[float]:
        """``array('d')``: parameter from the first corner."""
    @property
    def interior_positions(self) -> array[float]:
        """``array('d')``, 3 per interior point."""
    @property
    def junction_positions(self) -> array[float]:
        """``array('d')``, 3 per junction (empty for two regions)."""
    @property
    def junction_polygon(self) -> array[int]:
        """``array('i')``: the polygon of every junction (a triple point of the network)."""
    @property
    def junction_triangle(self) -> array[int]:
        """``array('i')``: the triangle of the supplied triangulation that contains it."""
    @property
    def junction_corners(self) -> array[int]:
        """``array('i')``, 3 per junction: the polygon corners of that triangle."""
    @property
    def junction_weights(self) -> array[float]:
        """``array('d')``, 3 per junction: the junction is ``sum(w[k] * x[corner k])``."""
    def lines(self) -> tuple[array[float], array[int], list[bool], array[int]]:
        """``(positions, starts, closed, regions)``: the paths joined where exactly two meet."""
    @property
    def n_clamped_crossings(self) -> int: ...
    @property
    def report(self) -> dict[str, int | float]:
        """``n_merged_nodes``, ``max_displacement``, ``max_displacement_relative``,
        ``n_moved_points``, ``n_removed_excursions`` (all zero on the caller's triangulation, empty for ``unwrap``)."""

@final
class InsertedMesh:
    """The polygon mesh with a cut inserted (the Rust reference of what an editor must produce)."""

    @property
    def positions(self) -> array[float]: ...
    @property
    def corner_vertices(self) -> array[int]: ...
    @property
    def polygon_starts(self) -> array[int]: ...
    @property
    def triangles(self) -> array[int]: ...
    @property
    def new_seams(self) -> array[int]:
        """``array('i')``, 2 per edge of the cut, ``[min, max]``, ascending."""
    @property
    def polygon_parent(self) -> array[int]: ...
    @property
    def vertex_source(self) -> array[int]:
        """``array('i')``, 3 per new vertex: kind (0 old, 1 edge point, 2 interior point,
        3 junction), then the old vertex / edge point / path / junction, and for kind 2 the
        position in the path's points."""
    @property
    def edge_point_vertex(self) -> array[int]: ...
    @property
    def junction_vertex(self) -> array[int]:
        """``array('i')``: the vertex of every junction of the description."""
    @property
    def polygon_region(self) -> array[int]:
        """``array('i')``: the region of every new polygon where the description says (a part
        next to a piece of the cut), else -1."""
    @property
    def changed_polygons(self) -> array[int]:
        """New polygons (more than three corners) that were split or got extra vertices: the
        ones whose triangulation an editor should keep. With ``changed_vertex_starts`` /
        ``changed_vertices`` (their vertices in order) and ``changed_triangle_starts`` /
        ``changed_triangles`` (3 vertices per triangle)."""
    @property
    def changed_vertex_starts(self) -> array[int]: ...
    @property
    def changed_vertices(self) -> array[int]: ...
    @property
    def changed_triangle_starts(self) -> array[int]: ...
    @property
    def changed_triangles(self) -> array[int]: ...
    @property
    def n_old_vertices(self) -> int: ...

@final
class CutUnwrapResult:
    @property
    def cut(self) -> CutDescription: ...
    @property
    def inserted(self) -> InsertedMesh: ...
    @property
    def seams(self) -> array[int]:
        """``array('i')``, 2 per seam of the new mesh: the old seams and the cut."""
    @property
    def uvs(self) -> array[float]:
        """``array('f')``, 2 per corner of the new mesh. Not packed."""
    @property
    def polygon_charts(self) -> array[int]: ...
    @property
    def cone_vertices(self) -> array[int]: ...
    @property
    def cone_angles(self) -> array[float]: ...
    @property
    def cut_edges(self) -> array[int]: ...
    @property
    def report(self) -> dict[str, bool | int | float | str | list[str]]: ...

Report = dict[str, bool | int | float | str | list["Report"]]

@final
class CutSession:
    """The cut flow on a polygon mesh. Weights, penalties, step sizes and the ripple wavelength
    refer to the surface scaled to unit area (``length_scale``); the state and the positions are in
    the unit of the mesh. A pointwise energy given as a callable is not offered."""

    def __init__(
        self,
        positions: Buffer,
        corner_vertices: Buffer,
        polygon_starts: Buffer,
        triangles: Buffer,
        *,
        n_regions: int | None = None,
        triangulation: Literal["intrinsic", "caller", "refined"] = "intrinsic",
        refine_bound: float | None = None,
        optimal_delaunay: bool = False,
        representation: Literal["integer", "signposts"] | None = None,
        min_offset: float | None = None,
        merge_tolerance: float | None = None,
        state: Buffer | None = None,
        seam_edges: Buffer | None = None,
        pinned_vertices: Buffer | None = None,
        forbidden_vertices: Buffer | None = None,
        prescribed_vertices: Buffer | None = None,
        prescribed_angles: Buffer | None = None,
        circle_vertex: int | None = None,
        circle_radius: float | None = None,
        polygon_labels: Buffer | None = None,
        corner_values: Buffer | None = None,
        cluster_iterations: int | None = None,
    ) -> None:
        """``triangulation``: the flow runs on an intrinsic triangulation of the surface flipped
        to Delaunay (``"intrinsic"``, the default: usually no vertex is added), on the same with
        refinement to ``refine_bound`` degrees (``"refined"``, default 25: vertices are added), or
        on the caller's triangulation (``"caller"``, as before). ``optimal_delaunay``: also build the
        intrinsic optimal Delaunay triangulation (the library's default options) on top of
        ``"intrinsic"`` or ``"refined"``: more even triangle areas in the intrinsic triangulation,
        much slower setup (seconds on 10k vertices, tens of seconds on 50k) and flow steps about 2.5
        times as slow;
        ``ValueError`` (``InputError``) with ``"caller"`` or together with ``state`` (the state
        carries its own flag). ``representation``: how the intrinsic triangulation is stored,
        ``"integer"`` (``None``, the default: integer coordinates, Gillespie, Sharp, Crane 2021) or
        ``"signposts"`` (Sharp, Soliman, Crane 2019); both make the same flips, splits and
        refinement (the signposts' rules), so the cut agrees to rounding; ``InputError`` with
        ``"caller"`` or together with ``state`` (the state carries its own). ``state``: bytes from
        ``snapshot()``, instead of a start (the triangulation and ``n_regions`` come from it;
        ``StaleStateError`` if the mesh or seams changed). ``min_offset`` (intrinsic and refined
        only; default the library's 1e-6): every point of the cut keeps at least this fraction of
        a side from the vertices and sides of the caller's triangle it lies in. ``merge_tolerance`` (intrinsic and refined only; default the library's 0.1 of an edge,
        0 switches it off): bends of the returned cut that deviate from the straight line between
        their neighbours by less than this fraction of the longest edge of their triangle are
        removed (the described cut can move by up to that much); ``CutDescription.report`` says
        how many. Exactly one start: ``circle_vertex`` with ``circle_radius`` (two regions only),
        ``polygon_labels`` (one label per polygon below ``n_regions``), ``corner_values``
        (``(n_regions, n_corners)``) or ``cluster_iterations`` (the vertex normals clustered into
        ``n_regions`` clusters in at most that many rounds: the start of the atlas). Any number
        of regions is accepted; ``NotAvailableError`` stays as a class but is not raised for
        this."""
    @property
    def triangulation(self) -> Literal["intrinsic", "caller", "refined"]: ...
    @property
    def representation(self) -> Literal["integer", "signposts"] | None:
        """The representation of the intrinsic triangulation (the state's after ``state=``;
        None for ``"caller"``)."""
    @property
    def optimal_delaunay(self) -> bool:
        """Whether the intrinsic optimal Delaunay triangulation is on (the state's flag after
        ``state=``; always False for ``"caller"``)."""
    @property
    def n_patches(self) -> int:
        """Connected pieces of the surface in one region (corners joined along polygons and
        across edges that are not seams)."""
    @property
    def state_is_per_corner(self) -> bool:
        """True when ``state`` (one value per polygon corner) is the whole state; False when the
        triangulation added vertices: keep ``snapshot()`` too."""
    @property
    def fingerprint(self) -> str:
        """Hex fingerprint of the triangulation the flow runs on."""
    @property
    def intrinsic_report(self) -> dict[str, bool | int | float | str] | None:
        """What the construction of the intrinsic triangulation did (``n_added``, flips,
        ``min_angle_degrees``, ``max_angle_degrees``, ``is_delaunay``, ``seconds``, ``n_optimized`` (vertices inserted by the optimal Delaunay
        triangulation, 0 if off), and with it on ``optimal_delaunay_status``,
        ``optimal_delaunay_iterations``, ``optimal_delaunay_last_motion``, ...); ``None`` for
        ``"caller"``."""
    def snapshot(self) -> bytes:
        """The whole state: levels on the vertices of the triangulation, its fingerprint, the
        number of regions, the triangulation option."""
    def restore(self, state: Buffer) -> None:
        """Sets the state from ``snapshot()`` bytes; ``StaleStateError`` for another mesh."""
    @property
    def n_regions(self) -> int: ...
    @property
    def n_vertices(self) -> int: ...
    @property
    def n_corners(self) -> int: ...
    @property
    def n_polygons(self) -> int: ...
    @property
    def n_triangles(self) -> int: ...
    @property
    def length_scale(self) -> float: ...
    @property
    def total_area(self) -> float: ...
    @property
    def time(self) -> float: ...
    @property
    def n_steps(self) -> int: ...
    @property
    def warnings(self) -> list[str]: ...
    @property
    def seams_cut(self) -> array[int]: ...
    @property
    def state(self) -> memoryview:
        """float64, shape ``(n_regions, n_corners)``; assign any float buffer of that shape or
        flat row-major."""
    @state.setter
    def state(self, value: Buffer) -> None: ...
    def state_row(self, region: int) -> array[float]: ...
    def corner_regions(self) -> array[int]: ...
    def restart(
        self,
        *,
        circle_vertex: int | None = None,
        circle_radius: float | None = None,
        polygon_labels: Buffer | None = None,
        corner_values: Buffer | None = None,
        cluster_iterations: int | None = None,
    ) -> None: ...
    def roughening_rate(self, **energy_options: object) -> float:
        """The roughening rate of the state; takes the keywords of ``step``."""
    def step(
        self,
        *,
        energy: Literal["dirichlet", "scaled_dirichlet", "hencky"] = "hencky",
        scale_weights: Buffer | None = None,
        feature_r_eps: float = 1.0,
        feature_smoothing: float = 0.02,
        length_weight: float | None = None,
        ripple_wavelength: float | None = None,
        penalty: Buffer | None = None,
        penalty_rule: Literal["gradient", "paper_diagonal"] | None = None,
        distortion_weights: Buffer | Sequence[float] | None = None,
        area_weights: Buffer | Sequence[float] | None = None,
        symmetry: Buffer | None = None,
        curvature_uniform: float | None = None,
        curvature_per_triangle: Buffer | None = None,
        stages: Literal["default", "paper", "recommended"] | None = None,
        flux_recovery: Literal["nodal", "smoothed"] | None = None,
        flux_smoothing: float | None = None,
        sliver_aspect_ratio: float | None = None,
        smoothing_passes: int | None = None,
        length_step: Literal["on_cut", "on_mesh"] | None = None,
        extension: Literal["local", "harmonic"] | None = None,
        speed_scaling: Literal["mean_gradient", "unit"] | None = None,
        redistance: Literal["keep_adjacent", "off"] | None = None,
        step_rule: Literal["adaptive", "capped", "fixed"] = "adaptive",
        time_scale: float | None = None,
        max_motion: float | None = None,
        tau_max: float | None = None,
        tau: float | None = None,
        min_crossing_offset: float | None = None,
    ) -> dict[str, bool | int | float | str]:
        """One flow step. ``None`` means the library's default. Neither ``length_weight`` nor
        ``ripple_wavelength``: ``ripple_wavelength=0.25``."""
    def run(
        self,
        *,
        energy: Literal["dirichlet", "scaled_dirichlet", "hencky"] = "hencky",
        scale_weights: Buffer | None = None,
        feature_r_eps: float = 1.0,
        feature_smoothing: float = 0.02,
        length_weight: float | None = None,
        ripple_wavelength: float | None = None,
        penalty: Buffer | None = None,
        penalty_rule: Literal["gradient", "paper_diagonal"] | None = None,
        distortion_weights: Buffer | Sequence[float] | None = None,
        area_weights: Buffer | Sequence[float] | None = None,
        symmetry: Buffer | None = None,
        curvature_uniform: float | None = None,
        curvature_per_triangle: Buffer | None = None,
        stages: Literal["default", "paper", "recommended"] | None = None,
        flux_recovery: Literal["nodal", "smoothed"] | None = None,
        flux_smoothing: float | None = None,
        sliver_aspect_ratio: float | None = None,
        smoothing_passes: int | None = None,
        length_step: Literal["on_cut", "on_mesh"] | None = None,
        extension: Literal["local", "harmonic"] | None = None,
        speed_scaling: Literal["mean_gradient", "unit"] | None = None,
        redistance: Literal["keep_adjacent", "off"] | None = None,
        step_rule: Literal["adaptive", "capped", "fixed"] = "adaptive",
        time_scale: float | None = None,
        max_motion: float | None = None,
        tau_max: float | None = None,
        tau: float | None = None,
        min_crossing_offset: float | None = None,
        max_steps: int = 1000,
        motion_tolerance: float = 1e-5,
        energy_tolerance: float = 1e-4,
        energy_window: int = 25,
    ) -> Report:
        """Steps until the stop rule says so; a tolerance of 0 never stops for that reason.
        ``steps`` in the result lists the dicts that ``step`` returns."""
    def preview(
        self, *, min_crossing_offset: float | None = None
    ) -> tuple[array[float], array[int], list[bool], array[int]]:
        """``(positions, starts, closed, regions)``: 3 floats per point, ``n_polylines + 1``
        offsets into the points, whether each polyline is closed, and 2 regions per polyline
        (``[lower, higher]``). For more than two regions the polylines are the curves of the
        network between its junctions (a curve that ends at a junction ends at its position)."""
    def lines(
        self, *, min_crossing_offset: float | None = None
    ) -> tuple[array[float], array[int], list[bool], array[int]]:
        """The cut as lines from its description (``CutDescription.lines()``)."""
    def describe_cut(self, *, min_crossing_offset: float | None = None) -> CutDescription: ...
    def insert(
        self, cut: CutDescription | None = None, *, min_crossing_offset: float | None = None
    ) -> InsertedMesh: ...
    def unwrap(
        self,
        *,
        min_crossing_offset: float | None = None,
        lambda_: float = 0.1,
        positive_cones_only: bool = True,
        optimize_boundary: bool = False,
        cones: bool = False,
    ) -> CutUnwrapResult: ...

def visibility_penalty(
    positions: Buffer,
    corner_vertices: Buffer,
    polygon_starts: Buffer,
    triangles: Buffer,
    camera: Sequence[float],
    *,
    visibility: Buffer | None = None,
) -> array[float]:
    """``eta_i = nu_i * max(n_i . e_i, 0)`` (paper Sec. 7.1), without a unit; no rays are cast."""

def length_weight_for_ripples(roughening_rate: float, wavelength: float) -> float:
    """``alpha_L = omega * wavelength / (2 pi)``."""


# ----------------------------------------------------------------------------------------------
# Applications of intrinsic triangulations: "Navigating Intrinsic Triangulations" (Sharp,
# Soliman, Crane 2019), Sec. 5.1, 5.3, 5.4 and 5.5 (crate ddg-intrinsic-apps).

@final
class GeodesicDistanceResult:
    @property
    def distances(self) -> array[float]:
        """``array('d')``, the geodesic distance at every vertex (0 at the sources)."""
    @property
    def report(self) -> dict[str, int | float]:
        """``n_triangulation_vertices``, ``time`` (the time step), ``seconds_factor``,
        ``seconds_solve``."""

def geodesic_distance(
    positions: Buffer,
    triangles: Buffer,
    sources: Buffer,
    *,
    triangulation: Literal["refined", "delaunay", "input"] = "refined",
    time_factor: float = 1.0,
    boundary: Literal["neumann", "dirichlet", "average"] = "neumann",
    representation: Literal["signposts", "integer"] = "integer",
) -> GeodesicDistanceResult:
    """Geodesic distance from the vertices ``sources`` by the heat method (Crane, Weischedel,
    Wardetzky 2013) on the intrinsic Delaunay refinement of the mesh (``"refined"``), its
    intrinsic Delaunay triangulation (``"delaunay"``) or the mesh itself (``"input"``). The
    time step is ``time_factor * h**2``, ``h`` the mean edge length of that triangulation.
    Raises ``InputError`` (bad arguments; a source on the boundary with ``"dirichlet"``),
    ``MeshError`` (not a valid manifold triangle mesh) and ``SolveError``."""

@final
class SteinerTreeResult:
    @property
    def points(self) -> array[float]:
        """``array('d')``, 3 per point of the polylines; every polyline lies on the surface."""
    @property
    def offsets(self) -> array[int]:
        """``array('i')``, where the polyline of every edge of the tree starts in ``points``
        (in points), and the number of points last."""
    @property
    def report(self) -> dict[str, bool | int | float]:
        """``length``, ``mst_length`` (the minimal spanning tree of the same triangulation),
        ``n_steiner_points``, ``is_disk`` (cutting along the tree and the boundary leaves one
        disk), ``max_interior_defect`` (0 when that disk is flat and unfolds)."""

def steiner_tree(
    positions: Buffer,
    triangles: Buffer,
    *,
    intrinsic: bool = True,
    steiner_points: bool = True,
    representation: Literal["signposts", "integer"] = "integer",
) -> SteinerTreeResult:
    """A short cut through all vertices: the minimal spanning tree of the intrinsic Delaunay
    triangulation (``intrinsic=True``) or of the mesh's edges, shortened by Fermat points
    (``steiner_points=True``) with the heuristic of Smith, Lee and Liebman (1981). Raises
    ``InputError``, ``MeshError`` and ``SolveError``."""

@final
class DirectionFieldResult:
    @property
    def vectors(self) -> array[float]:
        """``array('d')``, 3 per vertex: one of the field's ``n`` directions, a unit vector in
        the plane of an input face at the vertex (zero where the field vanishes)."""
    @property
    def singularity_points(self) -> array[float]:
        """``array('d')``, 3 per singularity: its place on the surface."""
    @property
    def singularity_indices(self) -> array[int]:
        """``array('i')``, per singularity: its index in units of ``1/n``."""
    @property
    def report(self) -> dict[str, str | int | float]:
        """``status`` (``"converged"`` or ``"iteration_limit"``), ``iterations``,
        ``eigenvalue``, ``energy`` (eigenvalue times area), ``residual``,
        ``n_triangulation_vertices``, ``index_sum`` (``n`` times the Euler characteristic on a
        closed surface), ``euler_characteristic``, ``seconds``."""

def smoothest_direction_field(
    positions: Buffer,
    triangles: Buffer,
    *,
    n: int = 1,
    triangulation: Literal["delaunay", "refined", "input"] = "delaunay",
    laplacian: Literal["fem", "fd"] = "fem",
    min_angle_degrees: float = 25.0,
    representation: Literal["integer", "signposts"] = "integer",
) -> DirectionFieldResult:
    """The smoothest ``n``-direction field (``n`` = 1 vectors, 2 lines, 4 crosses; Knöppel,
    Crane, Pinkall, Schröder 2013) on the intrinsic Delaunay triangulation of the mesh
    (``"delaunay"``), its intrinsic Delaunay refinement to ``min_angle_degrees``
    (``"refined"``) or the mesh itself (``"input"``), with the finite element (``"fem"``) or
    finite difference (``"fd"``) connection Laplacian. ``representation``: integer coordinates
    (the default; input-aligned frames computed from the integers, Integer Coordinates I26) or
    signposts (the stored angles); the field agrees to rounding, the singularities exactly.
    ``report["representation"]`` says which ran. Raises ``InputError``, ``MeshError`` and
    ``SolveError``."""

@final
class AdaptiveResult:
    @property
    def values(self) -> array[float]:
        """``array('d')``, ``n_components`` per vertex (interleaved)."""
    @property
    def n_components(self) -> int:
        """1 (``"green"``, ``"heat"``) or 2 (``"disk"``: coordinates in the unit disk)."""
    @property
    def report(self) -> dict[str, str | int | float]:
        """``status`` (``"tolerance"``, ``"vertex_limit"``, ``"iteration_limit"``,
        ``"stalled"``), ``iterations``, ``n_triangulation_vertices``, ``estimate``,
        ``oscillation``, ``seconds``."""

def adaptive_solve(
    positions: Buffer,
    triangles: Buffer,
    problem: Literal["green", "heat", "disk"],
    *,
    source: int = 0,
    time: float | None = None,
    max_vertices: int = 20000,
    bulk: float = 0.5,
    min_angle_degrees: float = 25.0,
    representation: Literal["signposts", "integer"] = "integer",
) -> AdaptiveResult:
    """Solve a model problem with adaptive refinement of an intrinsic triangulation (residual
    error estimator, bulk marking, intrinsic Delaunay refinement): the harmonic Green's
    function of ``source`` (``"green"``, mean zero), one backward Euler step of length ``time``
    of the heat equation from a unit point source (``"heat"``; default time the squared mean
    edge length), or a harmonic map of a disk-shaped mesh onto the unit disk (``"disk"``).
    Returns the solution at the mesh's vertices. Raises ``InputError``, ``MeshError`` (also a
    mesh that is not a disk for ``"disk"``) and ``SolveError``."""


# ---------------------------------------------------------------------------------------------
# Repulsive Curves (Yu, Schumacher, Crane 2021)
# ---------------------------------------------------------------------------------------------

Reals = Buffer | Sequence[float] | Sequence[Sequence[float]]
Indices = Buffer | Sequence[int] | Sequence[Sequence[int]]

RepulsiveSurface = dict[str, object]
"""An implicit surface by name: ``{"type": "sphere", "center": (0, 0, 0), "radius": 1.0}``;
``"torus"`` (``center``, ``major``, ``minor``), ``"plane"`` (``point``, ``normal``),
``"double_torus"`` (``radius2``), ``"rounded_box"`` (``center``, ``half``, ``radius``),
``"metaballs"`` (``balls``: ``[(center, amplitude, width)]``, ``threshold``), ``"union"``
(``parts``: list of surface dicts), ``"smooth_union"`` (``a``, ``b``, ``k``). The optional key
``"vertices"`` restricts the constraint to those vertices (default: all)."""

RepulsiveField = dict[str, object]
"""A vector field by name: ``{"type": "constant", "vector": (x, y, z), "normalize": False}``,
``{"type": "circular"}`` (around the y axis) or ``{"type": "sphere_tangent"}``."""

class RepulsiveInputError(InputError, ValueError):
    """An argument or an option of the repulsive-curves functions is invalid (the message is
    the library's). Both an ``InputError`` and a ``ValueError``."""

REPULSIVE_BACKENDS: Final[list[str]]

def repulsive_defaults() -> dict[str, object]:
    """The defaults of all options, read from the Rust library; the stubs below equal them."""

@final
class RepulsiveSession:
    """The repulsive flow of a curve network, one step at a time (the engine of a live Blender
    operator). ``step`` and ``run`` release the GIL; ``snapshot``, ``positions``, ``status``,
    ... never wait for a running step; ``cancel`` is safe from any thread.

    Options: ``alpha``, ``beta`` (energy exponents); ``length`` (``"edges"`` keeps every edge
    length, ``"total"`` the total length, ``"free"`` neither) with ``grow_factor`` and
    ``grow_rate`` (move the kept length to a factor of itself, ``rate`` of the start length per
    step); ``max_subdivisions``; ``pinned`` (vertices that do not move), ``tangent_vertices`` /
    ``tangent_edges`` (tangents that are kept), ``surface`` (vertices stay on an implicit
    surface), ``barycenter`` (default: fixed unless something is pinned), ``planar_axis``;
    ``obstacles`` (``(vertices, triangles)`` pairs, never crossed by the collision-guarded line
    search) with ``obstacle_weight`` > 0 to add their repulsive potential; ``total_length_weight``,
    ``length_difference_weight``, ``field`` with ``field_weight`` (penalties); ``line_search``
    and its constants (``line_search_start``: ``"grow_previous"``, the first trial is at most
    ``line_search_start_factor`` times the previous accepted step, or the paper's
    ``"two_thirds"`` of the longest collision-free step); ``backend`` (``"auto"``: exact energy
    and dense solver up to ``auto_dense_max_vertices`` vertices, Barnes-Hut energy and multigrid
    above, re-chosen after a subdivision; ``"dense"``; ``"barnes_hut"``; ``"multigrid"`` = exact
    + multigrid; ``"accelerated"`` = Barnes-Hut + multigrid) with the multigrid options
    ``mg_*`` (``mg_near_field_factor`` and ``mg_smoother_factor`` for curve networks with
    near-contacts, ``mg_max_fixed_fraction``, ``mg_energy_check``, ``mg_failure_residual``) and
    ``dense_fallback_max_vertices`` (a step whose multigrid solve fails repeats it with the dense
    solver up to that size; see ``solver_fallback`` of the step report);
    ``min_edge_length_ratio`` (a step that would shrink an edge below this fraction of the
    initial mean edge length raises ``SolveError`` instead of continuing to a singular metric;
    ``None`` disables); ``threads`` and
    ``dense_threads``; the projection and the stop rule (``"gradient"`` with ``stop_threshold``
    1e-4, ``"relative"`` 1e-3, or ``"stall"``).
    ``method``: ``"projected_gradient"`` (the default; hard constraints) or the soft methods
    ``"lbfgs"`` (``lbfgs_memory`` pairs, the paper's 30) and ``"ncg"``, which solve the paper's
    soft problem: no barycenter or length constraint, a Hencky length penalty of weight
    ``hencky_weight`` (the evaluation's 10) at the start lengths. ``pinned``, ``tangent_*``,
    ``surface`` and a length schedule (``grow_factor != 1``) need the hard constraints and are
    input errors with a soft method. Not in the paper, all with the library's defaults:
    ``near_pair`` (accurate quadrature of the energy for close edge pairs, blended between
    ``near_pair_blend_start`` and ``near_pair_blend_end`` edge lengths, ``near_pair_subpanels``;
    ``False`` is the paper's energy) with ``second_neighbours`` and ``second_neighbour_rule``
    (``"fold_gated"`` with ``second_neighbour_contact_band``, ``second_neighbour_cos_shut`` and
    ``second_neighbour_cos_open``; ``"one_sided"``; ``"two_sided"``; the gated rules
    ``"gated"``, ``"one_sided_gated"``, ``"one_sided_gated_either"`` require
    ``second_neighbour_width`` in (0, 1]); ``second_neighbour_band``; ``back_projection``,
    ``second_order_correction``, ``newton_polish`` (a factor of the projection tolerance, ``None``:
    the paper's chord method), ``projection_early_abort``, ``stuck_projected_reference``,
    ``stuck_retry_chord``; ``mg_max_threads`` (cap of the multigrid's thread pool, ``None``: none).
    Raises ``RepulsiveInputError`` (a ``ValueError``) for anything wrong in the input.
    """

    def __init__(
        self,
        positions: Reals,
        edges: Indices,
        *,
        alpha: float = 3.0,
        beta: float = 6.0,
        metric: Literal["sobolev", "l2", "h1", "h2"] = "sobolev",
        quadrature: Literal["midpoint", "trapezoid"] = "midpoint",
        length: Literal["edges", "total", "free"] = "edges",
        grow_factor: float = 1.0,
        grow_rate: float = 0.01,
        max_subdivisions: int = 0,
        barycenter: bool | None = None,
        pinned: Indices | None = None,
        tangent_vertices: Indices | None = None,
        tangent_edges: Indices | None = None,
        surface: RepulsiveSurface | None = None,
        planar_axis: Literal["x", "y", "z"] | None = None,
        obstacles: Sequence[tuple[Reals, Indices]] | None = None,
        obstacle_weight: float = 0.0,
        obstacle_exponent: float | None = None,
        obstacle_theta: float | None = None,
        total_length_weight: float = 0.0,
        length_difference_weight: float = 0.0,
        field: RepulsiveField | None = None,
        field_weight: float = 1.0,
        line_search: Literal["collision", "armijo"] = "collision",
        line_search_start: Literal["two_thirds", "grow_previous"] = "grow_previous",
        line_search_start_factor: float = 2.0,
        armijo: float = 1e-4,
        shrink: float = 0.5,
        max_halvings: int = 40,
        gap: float = 0.1,
        backend: Literal["auto", "dense", "barnes_hut", "multigrid", "accelerated"] = "auto",
    auto_dense_max_vertices: int = 2500,
    dense_fallback_max_vertices: int = 5000,
    mg_coarsest_vertices: int = 512,
    mg_max_levels: int = 20,
    mg_smoothing_steps: int = 12,
    mg_initial_smoothing_steps: int = 24,
    mg_tolerance: float = 1e-3,
    mg_max_cycles: int = 20,
    mg_flexible_cg: bool = True,
    mg_near_field_factor: float = 1.0,
    mg_smoother_factor: float = 1.0,
    mg_max_fixed_fraction: float = 0.4,
    mg_energy_check: bool = True,
    mg_failure_residual: float = 1e-2,
        auto_dense_max_vertices: int = 2500,
        dense_fallback_max_vertices: int = 5000,
        mg_coarsest_vertices: int = 512,
        mg_max_levels: int = 20,
        mg_smoothing_steps: int = 12,
        mg_initial_smoothing_steps: int = 24,
        mg_tolerance: float = 1e-3,
        mg_max_cycles: int = 20,
        mg_flexible_cg: bool = True,
        mg_near_field_factor: float = 1.0,
        mg_smoother_factor: float = 1.0,
        mg_max_fixed_fraction: float = 0.4,
        mg_energy_check: bool = True,
        mg_failure_residual: float = 1e-2,
        bh_theta_x: float = 0.125,
        bh_theta_t: float = 0.25,
        bh_leaf_size: int = 8,
        threads: int | None = None,
        dense_threads: int | None = None,
    min_edge_length_ratio: float | None = 1e-3,
        min_edge_length_ratio: float | None = 1e-3,
        projection_tolerance: float = 1e-4,
        projection_absolute: bool = False,
        projection_iterations: int = 3,
        stop: Literal["gradient", "relative", "stall"] = "gradient",
        stop_threshold: float | None = None,
        stall_window: int = 10,
        gradient_norm: Literal["euclidean", "mass"] = "euclidean",
        max_steps: int = 1000,
        time_limit: float | None = None,
        method: Literal["projected_gradient", "lbfgs", "ncg"] = "projected_gradient",
        lbfgs_memory: int = 30,
        hencky_weight: float = 10.0,
        near_pair: bool = True,
        near_pair_blend_start: float = 0.25,
        near_pair_blend_end: float = 0.75,
        near_pair_subpanels: int = 4,
        second_neighbours: bool = True,
        second_neighbour_rule: Literal["fold_gated", "one_sided", "two_sided", "gated", "one_sided_gated", "one_sided_gated_either"] = "fold_gated",
        second_neighbour_width: float | None = None,
        second_neighbour_contact_band: tuple[float, float] = (0.001, 0.005),
        second_neighbour_cos_shut: float = 0.0,
        second_neighbour_cos_open: float = 0.5,
        second_neighbour_band: tuple[float, float] | None = None,
        back_projection: bool = True,
        second_order_correction: bool = False,
        newton_polish: float | None = 100.0,
        projection_early_abort: bool = True,
        stuck_projected_reference: bool = False,
        stuck_retry_chord: bool = True,
        mg_max_threads: int | None = 16,
    ) -> None: ...
    def step(self) -> dict[str, object]:
        """One step. Keys: ``iteration``, ``status`` (``"accepted"``, ``"stuck"``,
        ``"converged"``, ``"retarget_stalled"``, ``"cancelled"``), ``session_status``,
        ``n_vertices``, ``n_edges``, ``energy_before``, ``energy_after``, ``gradient_norm``,
        ``tau``, ``halvings``, ``max_violation``, ``linear``, ``retarget``, ``subdivided``,
        ``energy_backend``, ``solver_backend``, ``solver_fallback`` (``None``, or a dict with
        ``from``, ``to``, ``cause`` = ``"multigrid_diverged"`` / ``"multigrid_not_converged"``,
        ``message``, ``dense_solves``, ``build_seconds``: the multigrid solve failed and the
        dense solver answered the rest of the step), ``schedule`` (``None``, or ``progress``,
        ``span``, ``fraction``, ``factor_target``, ``factor_imposed``, ``factor_measured``,
        ``finished``: the length schedule), ``seconds`` ..."""
    def run(
        self,
        max_steps: int | None = None,
        progress: Callable[[dict[str, object]], bool | None] | None = None,
        time_limit: float | None = None,
    ) -> dict[str, object]:
        """Step until final, ``max_steps``, ``cancel()``, ``time_limit`` seconds, or
        ``progress(report)`` returns ``False``. Keys: ``end``, ``status``, ``steps_taken``,
        ``iteration``, ``energy_before``, ``energy_after``, ``seconds``, ``steps``."""
    def snapshot(self) -> tuple[array[float], array[int], int, float, str]:
        """``(positions, edges, iteration, energy, status)`` of the last completed step."""
    def cancel(self) -> None: ...
    def reset(self) -> None: ...
    def set_schedule(self, factor: float, rate: float = 0.01) -> None: ...
    @property
    def positions(self) -> array[float]:
        """``array('d')``, 3 per vertex."""
    @property
    def edges(self) -> array[int]:
        """``array('i')``, 2 per edge."""
    @property
    def energy(self) -> float: ...
    @property
    def iteration(self) -> int: ...
    @property
    def status(self) -> Literal[
        "running", "converged", "max_steps", "time_limit", "stuck", "cancelled", "error"
    ]: ...
    @property
    def error(self) -> str | None: ...
    @property
    def n_vertices(self) -> int: ...
    @property
    def n_edges(self) -> int: ...
    @property
    def subdivisions(self) -> int: ...
    @property
    def mean_edge_length0(self) -> float: ...
    @property
    def total_length(self) -> float: ...
    @property
    def vertex_source(self) -> array[int]:
        """``array('i')``: the input vertex each vertex is, -1 if created by subdivision.
        Old vertices keep their indices; the midpoint of edge ``I`` is vertex ``|V| + I`` and
        edge ``I`` becomes edges ``I`` and ``|E| + I``."""
    @property
    def edge_source(self) -> array[int]:
        """``array('i')``: the input edge each edge descends from."""
    @property
    def report(self) -> dict[str, object] | None:
        """The dict of the last step."""

def repulsive_flow(
    positions: Reals,
    edges: Indices,
    steps: int | None = None,
    progress: Callable[[dict[str, object]], bool | None] | None = None,
    *,
    alpha: float = 3.0,
    beta: float = 6.0,
    metric: Literal["sobolev", "l2", "h1", "h2"] = "sobolev",
    quadrature: Literal["midpoint", "trapezoid"] = "midpoint",
    length: Literal["edges", "total", "free"] = "edges",
    grow_factor: float = 1.0,
    grow_rate: float = 0.01,
    max_subdivisions: int = 0,
    barycenter: bool | None = None,
    pinned: Indices | None = None,
    tangent_vertices: Indices | None = None,
    tangent_edges: Indices | None = None,
    surface: RepulsiveSurface | None = None,
    planar_axis: Literal["x", "y", "z"] | None = None,
    obstacles: Sequence[tuple[Reals, Indices]] | None = None,
    obstacle_weight: float = 0.0,
    obstacle_exponent: float | None = None,
    obstacle_theta: float | None = None,
    total_length_weight: float = 0.0,
    length_difference_weight: float = 0.0,
    field: RepulsiveField | None = None,
    field_weight: float = 1.0,
    line_search: Literal["collision", "armijo"] = "collision",
    line_search_start: Literal["two_thirds", "grow_previous"] = "grow_previous",
    line_search_start_factor: float = 2.0,
    armijo: float = 1e-4,
    shrink: float = 0.5,
    max_halvings: int = 40,
    gap: float = 0.1,
    backend: Literal["auto", "dense", "barnes_hut", "multigrid", "accelerated"] = "auto",
    auto_dense_max_vertices: int = 2500,
    dense_fallback_max_vertices: int = 5000,
    mg_coarsest_vertices: int = 512,
    mg_max_levels: int = 20,
    mg_smoothing_steps: int = 12,
    mg_initial_smoothing_steps: int = 24,
    mg_tolerance: float = 1e-3,
    mg_max_cycles: int = 20,
    mg_flexible_cg: bool = True,
    mg_near_field_factor: float = 1.0,
    mg_smoother_factor: float = 1.0,
    mg_max_fixed_fraction: float = 0.4,
    mg_energy_check: bool = True,
    mg_failure_residual: float = 1e-2,
    bh_theta_x: float = 0.125,
    bh_theta_t: float = 0.25,
    bh_leaf_size: int = 8,
    threads: int | None = None,
    dense_threads: int | None = None,
    min_edge_length_ratio: float | None = 1e-3,
    projection_tolerance: float = 1e-4,
    projection_absolute: bool = False,
    projection_iterations: int = 3,
    stop: Literal["gradient", "relative", "stall"] = "gradient",
    stop_threshold: float | None = None,
    stall_window: int = 10,
    gradient_norm: Literal["euclidean", "mass"] = "euclidean",
    max_steps: int = 1000,
    time_limit: float | None = None,
    method: Literal["projected_gradient", "lbfgs", "ncg"] = "projected_gradient",
    lbfgs_memory: int = 30,
    hencky_weight: float = 10.0,
    near_pair: bool = True,
    near_pair_blend_start: float = 0.25,
    near_pair_blend_end: float = 0.75,
    near_pair_subpanels: int = 4,
    second_neighbours: bool = True,
    second_neighbour_rule: Literal["fold_gated", "one_sided", "two_sided", "gated", "one_sided_gated", "one_sided_gated_either"] = "fold_gated",
    second_neighbour_width: float | None = None,
    second_neighbour_contact_band: tuple[float, float] = (0.001, 0.005),
    second_neighbour_cos_shut: float = 0.0,
    second_neighbour_cos_open: float = 0.5,
    second_neighbour_band: tuple[float, float] | None = None,
    back_projection: bool = True,
    second_order_correction: bool = False,
    newton_polish: float | None = 100.0,
    projection_early_abort: bool = True,
    stuck_projected_reference: bool = False,
    stuck_retry_chord: bool = True,
    mg_max_threads: int | None = 16,
) -> tuple[array[float], array[int], dict[str, object]]:
    """The repulsive flow in one call: ``(positions, edges, report)``. ``steps`` overrides
    ``max_steps``. Report keys: ``status``, ``iterations``, ``energy_initial``,
    ``energy_final``, ``seconds``, ``n_subdivisions``, ``vertex_source``, ``edge_source``,
    ``steps`` (the dicts of ``RepulsiveSession.step``). A session stepped the same number of
    times gives bitwise the same curve."""

def curve_energy(
    positions: Reals,
    edges: Indices,
    alpha: float = 3.0,
    beta: float = 6.0,
    *,
    backend: Literal["exact", "barnes_hut"] = "exact",
    bh_theta_x: float = 0.125,
    bh_theta_t: float = 0.25,
    bh_leaf_size: int = 8,
    threads: int | None = None,
) -> tuple[float, array[float]]:
    """The tangent-point energy and its differential (``array('d')``, 3 per vertex)."""


# ---- IntrinsicTriangulation: integer coordinates (Gillespie, Sharp, Crane 2021) -------------

class Places(TypedDict):
    """Points of a mesh as arrays, one entry per point."""

    kind: array[int]
    """``array('i')``: 0 vertex, 1 edge, 2 face."""
    id: array[int]
    """``array('i')``: the vertex, edge or face id."""
    params: array[float]
    """``array('d')``, 3 per point: ``(0, 0, 0)`` for a vertex, ``(t, 0, 0)`` for an edge,
    barycentric coordinates for a face."""

class VertexPlaces(Places):
    positions: array[float]
    """``array('d')``, 3 per vertex of the triangulation: its place in space."""

class Counts(TypedDict):
    n_vertices: int
    n_edges: int
    n_faces: int
    n_input_vertices: int
    n_input_faces: int
    min_angle: float
    max_angle: float
    n_crossings: int

class FlipReport(TypedDict):
    status: Literal["delaunay", "flip_limit", "not_flippable"]
    n_flips: int
    n_not_delaunay: int

class RefineReport(TypedDict):
    status: Literal["done", "insertion_limit"]
    n_inserted: int
    n_face_insertions: int
    n_edge_splits: int
    n_flips: int
    n_removals: int
    n_removals_refused: int
    min_angle_valid: float
    """Radians; narrow-cone faces are left out."""
    n_below_bound: int
    certified: bool

class CommonSubdivision(TypedDict):
    positions: array[float]
    """``array('d')``, 3 per vertex of the subdivision, on the input surface. The vertices of
    the triangulation come first (``n_intrinsic_vertices`` of them, with their ids), then one
    per crossing of an intrinsic edge with an input edge."""
    on_input: Places
    on_triangulation: Places
    polygon_starts: array[int]
    """``array('i')``, where each polygon starts in ``corners``; one more entry at the end."""
    corners: array[int]
    input_face: array[int]
    """Per polygon: the input face it lies in."""
    intrinsic_face: array[int]
    """Per polygon: the face of the triangulation it lies in."""
    n_intrinsic_vertices: int

class EdgePolylines(TypedDict):
    edges: array[int]
    points: array[float]
    """``array('d')``, 3 per point; every polyline lies on the input surface."""
    offsets: array[int]
    """``array('i')``, where each polyline starts in ``points``; one more entry at the end."""
    along_input_edge: list[int]
    """Per polyline: the input edge it lies along, or -1."""
    n_crossings: list[int]
    max_relative_mismatch: float

class Certificate(TypedDict):
    n_faces: int
    n_along: int
    n_crossings: int
    n_input_edges: int
    n_inserted_on_edges: int
    max_crossings_per_input_edge: int

@final
class IntrinsicTriangulation:
    """An intrinsic triangulation of a mesh on integer coordinates (default; Gillespie, Sharp,
    Crane 2021) or signposts (Sharp, Soliman, Crane 2019). Starts as the input mesh; the input
    is never changed. Vertex ids below ``n_input_vertices`` are the input's. Raises
    ``InputError``, ``MeshError`` and ``SolveError``."""

    def __init__(
        self,
        positions: Buffer,
        triangles: Buffer,
        representation: Literal["integer", "signposts"] = "integer",
        mollify: float | None = None,
    ) -> None:
        """``mollify``: relative margin of intrinsic mollification (integer only; the paper
        uses ``1e-4``)."""
    @property
    def representation(self) -> Literal["integer", "signposts"]: ...
    def counts(self) -> Counts: ...
    def flip_to_delaunay(self) -> FlipReport: ...
    def refine(
        self,
        min_angle_degrees: float = 25.0,
        rules: Literal["chew", "ruppert"] = "chew",
        max_insertions: int | None = None,
        max_circumradius: float | None = None,
    ) -> RefineReport:
        """Delaunay refinement to an angle bound of at most 35 degrees. ``"chew"`` are the
        integer paper's rules, ``"ruppert"`` the signpost paper's."""
    def triangles(self) -> tuple[array[int], array[float]]:
        """(``array('i')`` 3 vertex ids per face, ``array('d')`` one length per edge)."""
    def vertex_places(self) -> VertexPlaces: ...
    def common_subdivision(self) -> CommonSubdivision:
        """The common subdivision of the input mesh and the triangulation (Fig. 12 of the
        integer paper): convex polygons with their input face and intrinsic face."""
    def edge_polylines(self, edges: Buffer | None = None) -> EdgePolylines:
        """The intrinsic edges (default all) as polylines on the input mesh."""
    def certify(self) -> Certificate:
        """The certificate of the integer coordinates; ``InputError`` for signposts,
        ``SolveError`` naming the item that fails."""
    def to_input(self, values: Buffer) -> array[float]:
        """Transfer a function on the triangulation's vertices to the input vertices (L2
        projection). Integer coordinates only."""
    def to_intrinsic(self, values: Buffer) -> array[float]:
        """Transfer a function on the input vertices to the triangulation's vertices."""


# ---- Geodesics by edge flips: "You Can Find Geodesic Paths in Triangle Meshes by Just Flipping
# Edges" (Sharp, Crane 2020), crate ddg-flip-geodesics ------------------------------------------

class GeodesicCertificate(TypedDict):
    """The certificate of the triangulation ``T`` (every joint straight or exempt)."""

    geodesic: bool
    tolerance: float
    n_joints: int
    n_pinned: int
    n_straight: int
    n_blocked: int
    min_angle: float
    failures: list[tuple[int, float, bool]]
    """Failing joints: (vertex of ``T``, angle in radians, flexible)."""

class InputGeodesicCertificate(TypedDict):
    """The certificate on the input mesh, computed from the input alone."""

    passed: bool
    n_crossings: int
    max_crossing_deviation: float
    n_joints: int
    n_exempt: int
    min_joint_angle: float
    max_wedge_mismatch: float
    max_length_error: float
    n_faults: int
    faults: list[str]

class BezierReport(TypedDict):
    rounds: int
    control_points: array[int]
    """Vertices of the final triangulation: the control polygon of the result."""
    n_splits: int
    n_snapped: int
    n_delaunay_retries: int
    n_nudged: int
    n_shortenings: int
    n_not_geodesic: int

class GeodesicReport(TypedDict):
    representation: Literal["signposts", "integer"]
    status: Literal[
        "geodesic", "stalled", "iteration_limit", "length_ratio", "min_angle", "sweep_limit"
    ]
    length: float
    initial_length: float
    iterations: int
    flips: int
    refused_flips: int
    stalled: int
    n_paths: int
    n_crossings: int
    max_relative_mismatch: float
    n_triangulation_vertices: int
    certificate: GeodesicCertificate
    input_certificate: InputGeodesicCertificate | None
    """None with ``certify_input=False``."""
    bezier: BezierReport | None
    """Only for ``geodesic_bezier``."""

@final
class GeodesicResult:
    @property
    def polylines(self) -> list[array[float]]:
        """Per path ``array('d')``, 3 per point on the input surface (numpy:
        ``asarray(p).reshape(-1, 3)``). A closed path repeats its first point."""
    @property
    def surface_points(self) -> list[Places]:
        """Per path the places of its points on the input mesh (vertex, edge with ``t`` from the
        lower-numbered vertex, face with barycentric coordinates)."""
    @property
    def lengths(self) -> array[float]: ...
    @property
    def closed(self) -> list[bool]: ...
    @property
    def report(self) -> GeodesicReport: ...
    @property
    def control_polygon(self) -> array[float] | None:
        """Bézier: ``array('d')``, 3 per control point of the result."""
    @property
    def triangulation_edges(self) -> EdgePolylines | None:
        """With ``triangulation_edges=True``: the edges of the final triangulation as polylines
        on the input mesh."""

def geodesic_path(
    positions: Buffer,
    triangles: Buffer,
    vertices: Buffer,
    *,
    closed: bool = False,
    pin_waypoints: bool = True,
    representation: Literal["signposts", "integer"] | None = None,
    mollify: float | None = None,
    tolerance: float = 1e-5,
    max_iterations: int | None = None,
    stop_length_ratio: float | None = None,
    certify_input: bool = True,
    triangulation_edges: bool = False,
) -> GeodesicResult:
    """A geodesic path (or loop, ``closed=True``) through the waypoints ``vertices`` by FlipOut.
    ``representation=None`` is the default of all consumers (integer coordinates; signposts
    on request). ``pin_waypoints=False``: the waypoints only seed the
    path. Raises ``InputError`` (the message names the parameter), ``MeshError`` and
    ``SolveError`` (the path crosses itself; a failed extraction on signposts says to try
    ``representation="integer"``)."""

def geodesic_network(
    positions: Buffer,
    triangles: Buffer,
    edges: Buffer,
    *,
    pinned: Buffer | None = None,
    representation: Literal["signposts", "integer"] | None = None,
    mollify: float | None = None,
    tolerance: float = 1e-5,
    max_iterations: int | None = None,
    stop_length_ratio: float | None = None,
    certify_input: bool = True,
    triangulation_edges: bool = False,
) -> GeodesicResult:
    """The network of mesh edges ``edges`` (vertex pairs, shape ``(m, 2)``) shortened without
    crossings; junctions and ``pinned`` vertices stay. One path per branch or loop."""

def geodesic_bezier(
    positions: Buffer,
    triangles: Buffer,
    control_vertices: Buffer,
    *,
    rounds: int = 4,
    representation: Literal["signposts", "integer"] | None = None,
    mollify: float | None = None,
    tolerance: float = 1e-5,
    max_iterations: int | None = None,
    stop_length_ratio: float | None = None,
    certify_input: bool = True,
    triangulation_edges: bool = False,
) -> GeodesicResult:
    """The geodesic Bézier curve with the given control points, ``rounds`` of midpoint
    subdivision (1 to 8 suit interactive use). ``SolveError`` when the shortest paths between
    control points touch or cross, or a control point is repeated."""

# Surface simplification using intrinsic error metrics (Liu, Gillespie, Chislett, Sharp, Jacobson,
# Crane, ACM TOG 42(4), 2023; crate ddg-intrinsic-simplification).

Coo: TypeAlias = tuple[array[int], array[int], array[float]]
"""A sparse matrix as triplets ``(rows, cols, values)``."""

@final
class Simplification:
    """A coarse intrinsic triangulation of a surface and the map from the input onto it. Tracked
    points are the input vertices (ids ``0 .. n_input_vertices - 1``) followed by the extra
    ``points`` given to ``simplify``."""
    @property
    def n_input_vertices(self) -> int: ...
    @property
    def n_vertices(self) -> int: ...
    @property
    def n_faces(self) -> int: ...
    @property
    def vertices(self) -> array[int]:
        """Fine id of every coarse vertex (ids at or above ``n_input_vertices`` were inserted
        by refinement)."""
    @property
    def faces(self) -> array[int]:
        """3 fine vertex ids per coarse face (a face of the Δ-complex may repeat a vertex)."""
    @property
    def faces_compact(self) -> array[int]:
        """3 coarse vertex ids (``0 .. n_vertices - 1``, the index into ``vertices``) per face."""
    @property
    def lengths(self) -> array[float]:
        """3 side lengths per face: corner ``k`` to corner ``k + 1``."""
    @property
    def point_face(self) -> array[int]:
        """The coarse face of every tracked point (a point that is a coarse vertex: its lowest
        incident face), -1 for none."""
    @property
    def point_bary(self) -> array[float]:
        """3 barycentric coordinates per tracked point, in the corners of ``point_face``."""
    @property
    def point_vertex(self) -> array[int]:
        """The compacted coarse vertex a tracked point is, -1 for a point inside a face."""
    @property
    def face_colours(self) -> array[int]:
        """Greedy colouring of the coarse faces: faces that share a vertex differ."""
    @property
    def n_colours(self) -> int: ...
    @property
    def point_colours(self) -> array[int]:
        """The colour of the coarse face of every tracked point, -1 for none."""
    @property
    def report(self) -> dict[str, object]:
        """``stop``, vertex counts, removals by kind, flips, skips and failures by reason,
        ``components``, ``keys``, ``timings``, ``seconds``, and ``anisotropy``,
        ``pre_refinement``, ``post_refinement`` (dicts or None)."""
    def prolongation(self) -> Coo:
        """``P`` (input vertices by coarse vertices), rows sum to one."""
    def prolongation_scipy(self) -> object:
        """``P`` as a ``scipy.sparse`` CSR matrix; ``ImportError`` without scipy and numpy."""
    def point_prolongation(self) -> Coo:
        """``P`` for every tracked point (extra points included)."""
    def vector_prolongation(self, n: int = 1) -> dict[str, object]:
        """``P^vec`` for ``n``-direction fields: ``rows``, ``cols``, ``re``, ``im``, ``shape``,
        ``report``. Not available with refinement or anisotropy (``InputError``)."""
    def vector_prolongation_scipy(self, n: int = 1) -> object:
        """``P^vec`` as a complex ``scipy.sparse`` CSR matrix."""
    def poisson(
        self, f: Buffer, fixed: Buffer | None = None, values: Buffer | None = None
    ) -> dict[str, object]:
        """Poisson's equation on the coarse mesh, prolongated: ``coarse``, ``fine``,
        ``report``. Without ``fixed`` the mean-zero problem on a closed surface."""
    def geodesic(
        self, source: int, method: Literal["exact", "heat"] = "exact", time_factor: float = 1.0
    ) -> dict[str, array[float]]:
        """Geodesic distance from a tracked point on the coarse mesh: ``coarse``, ``fine``."""
    def voronoi(self, sources: Buffer, threads: int = 4) -> dict[str, array[float] | array[int]]:
        """Geodesic Voronoi labels: ``fine_labels``, ``coarse_labels``, ``fine_distance``."""

def simplify(
    positions: Buffer,
    triangles: Buffer,
    *,
    target: int | None = None,
    ratio: float | None = None,
    weights: float = 0.0,
    memory: bool = True,
    fixed: Buffer | None = None,
    masses: Buffer | None = None,
    mass_weight: float = 1.0,
    points: tuple[Buffer, Buffer] | None = None,
    refine_before: bool = False,
    refine_after: bool = False,
    repair_simplicial: bool = True,
    anisotropy: float | None = None,
    anisotropy_direction: Literal["max", "min"] | Buffer | None = None,
    anisotropy_clamp: bool = False,
    representation: Literal["signposts", "integer"] = "integer",
) -> Simplification:
    """Simplifies a mesh by intrinsic vertex removal with the intrinsic error metric (Alg. 1).
    Exactly one of ``target`` (vertex count) and ``ratio`` is required. ``weights`` is the
    curvature/area weight ``w`` (0: curvature, 1: area). ``points`` is ``(faces, bary)``: extra
    tracked points. Raises ``InputError``, ``MeshError``, ``SolveError``."""

@final
class MultigridSolver:
    @property
    def info(self) -> dict[str, object]: ...
    def solve(
        self, b: Buffer, known_values: Buffer | None = None, x0: Buffer | None = None
    ) -> tuple[array[float], dict[str, object]]: ...
    def pcg(
        self,
        b: Buffer,
        known_values: Buffer | None = None,
        x0: Buffer | None = None,
        max_iterations: int = 100,
    ) -> tuple[array[float], dict[str, object]]: ...

@final
class Hierarchy:
    @property
    def n_levels(self) -> int: ...
    @property
    def level_sizes(self) -> list[int]: ...
    @property
    def report(self) -> dict[str, object]: ...
    def level_vertices(self, level: int) -> array[int]: ...
    def faces(self, level: int) -> array[int]: ...
    def laplacian(self, level: int) -> Coo: ...
    def mass(self, level: int) -> array[float]: ...
    def prolongation(self, level: int) -> Coo:
        """``P_k``: level ``k`` by level ``k + 1``."""
    def multigrid(
        self,
        matrix: tuple[Buffer, Buffer, Buffer] | None = None,
        *,
        shift: float = 0.0,
        known: Buffer | None = None,
        pre_sweeps: int = 2,
        post_sweeps: int = 2,
        tolerance: float = 1e-10,
        max_cycles: int = 50,
    ) -> MultigridSolver:
        """The multigrid solver for ``matrix`` (default ``L + shift M`` of level 0)."""

def build_hierarchy(
    positions: Buffer,
    triangles: Buffer,
    *,
    ratio: float = 0.25,
    n_min: int = 500,
    fresh_levels: bool = False,
    weights: float = 0.0,
    memory: bool = True,
    fixed: Buffer | None = None,
) -> Hierarchy: ...

@final
class QemResult:
    @property
    def positions(self) -> array[float]: ...
    @property
    def triangles(self) -> array[int]: ...
    @property
    def vertex_map(self) -> array[int]: ...
    @property
    def projection_face(self) -> array[int]: ...
    @property
    def projection_bary(self) -> array[float]: ...
    @property
    def projection_distance(self) -> array[float]: ...
    @property
    def report(self) -> dict[str, int | float | bool]: ...

def qem_simplify(
    positions: Buffer,
    triangles: Buffer,
    *,
    target: int,
    area_weighted: bool = False,
    boundary_weight: float = 1000.0,
) -> QemResult:
    """Quadric error metric simplification (Garland, Heckbert 1997), the paper's extrinsic
    baseline, with the closest-point projection of every input vertex."""
