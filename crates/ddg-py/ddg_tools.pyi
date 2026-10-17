"""Type stubs of the extension module ``ddg_tools`` (built from ``crates/ddg-py``).

Array arguments accept any C-contiguous object with the buffer protocol; see README.md for the
element types. Results are ``array.array`` objects.
"""

from array import array
from collections.abc import Buffer, Sequence
from typing import Final, Literal, final

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
        carries its own flag). ``state``: bytes from
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
) -> DirectionFieldResult:
    """The smoothest ``n``-direction field (``n`` = 1 vectors, 2 lines, 4 crosses; Knöppel,
    Crane, Pinkall, Schröder 2013) on the intrinsic Delaunay triangulation of the mesh
    (``"delaunay"``), its intrinsic Delaunay refinement to ``min_angle_degrees``
    (``"refined"``) or the mesh itself (``"input"``), with the finite element (``"fem"``) or
    finite difference (``"fd"``) connection Laplacian. Raises ``InputError``, ``MeshError``
    and ``SolveError``."""

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
) -> AdaptiveResult:
    """Solve a model problem with adaptive refinement of an intrinsic triangulation (residual
    error estimator, bulk marking, intrinsic Delaunay refinement): the harmonic Green's
    function of ``source`` (``"green"``, mean zero), one backward Euler step of length ``time``
    of the heat equation from a unit point source (``"heat"``; default time the squared mean
    edge length), or a harmonic map of a disk-shaped mesh onto the unit disk (``"disk"``).
    Returns the solution at the mesh's vertices. Raises ``InputError``, ``MeshError`` (also a
    mesh that is not a disk for ``"disk"``) and ``SolveError``."""
