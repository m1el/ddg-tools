"""Type stubs of the extension module ``ddg_tools`` (built from ``crates/ddg-py``).

Array arguments accept any C-contiguous object with the buffer protocol; see README.md for the
element types. Results are ``array.array`` objects.
"""

from array import array
from collections.abc import Buffer
from typing import Final, final

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
) -> UnwrapResult: ...
def pack_charts(
    uvs: Buffer,
    polygon_starts: Buffer,
    polygon_charts: Buffer,
    *,
    margin: float = 0.0,
) -> tuple[array[float], dict[str, int | float]]: ...
