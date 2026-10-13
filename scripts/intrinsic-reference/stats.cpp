// Statistics driver for the reference runs of "Navigating Intrinsic Triangulations" (Sharp, Soliman, Crane 2019).
//
// The authors' demo (MIT) performs the operations but prints no counts. This driver does the same calls as the
// demo's main.cpp (same reader, same SignpostIntrinsicTriangulation, same flipToDelaunay / delaunayRefine with the
// same arguments) and counts flips and insertions through geometry-central's own callback lists. It writes one
// "key value" pair per line to stdout. It does not write any triangulation: the demo does that. Both are
// deterministic, and the script checks that the face and vertex counts of the two agree.
//
// usage: stats mesh.obj flip
//        stats mesh.obj refine ANGLE_DEG [CIRCUMRADIUS|inf] [MAX_INSERTIONS]
//   MAX_INSERTIONS follows the demo: 0 = none, negative = |n| * nVertices (demo default -10).
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

using namespace geometrycentral;
using namespace geometrycentral::surface;

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: stats mesh.obj flip | refine ANGLE [CIRCUMRADIUS|inf] [MAX_INSERTIONS]\n";
    return 2;
  }
  std::string op = argv[2];
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(argv[1]);

  size_t nV0 = mesh->nVertices(), nF0 = mesh->nFaces(), nE0 = mesh->nEdges();
  SignpostIntrinsicTriangulation tri(*mesh, *geometry);

  size_t flips = 0, faceIns = 0, edgeSplits = 0;
  tri.edgeFlipCallbackList.push_back([&](Edge) { flips++; });
  tri.faceInsertionCallbackList.push_back([&](Face, Vertex) { faceIns++; });
  tri.edgeSplitCallbackList.push_back([&](Edge, Halfedge, Halfedge) { edgeSplits++; });

  auto t0 = std::chrono::steady_clock::now();
  size_t maxIns = INVALID_IND;
  if (op == "flip") {
    tri.flipToDelaunay();
  } else if (op == "refine") {
    double angle = std::atof(argv[3]);
    double circ = std::numeric_limits<double>::infinity();
    if (argc > 4 && std::string(argv[4]) != "inf") circ = std::atof(argv[4]);
    int mi = -10;
    if (argc > 5) mi = std::atoi(argv[5]);
    if (mi < 0) mi *= -(int)nV0;       // as the demo: scale by nVerts
    maxIns = (mi != 0) ? (size_t)mi : INVALID_IND;
    tri.delaunayRefine(angle, circ, maxIns);
  } else {
    std::cerr << "unknown operation\n";
    return 2;
  }
  double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  std::cout.precision(17);
  std::cout << "nV_input " << nV0 << "\nnF_input " << nF0 << "\nnE_input " << nE0 << "\n";
  std::cout << "nV " << tri.mesh.nVertices() << "\nnF " << tri.mesh.nFaces() << "\nnE " << tri.mesh.nEdges() << "\n";
  std::cout << "flips " << flips << "\nface_insertions " << faceIns << "\nedge_splits " << edgeSplits << "\n";
  std::cout << "max_insertions_allowed " << (maxIns == INVALID_IND ? -1 : (long long)maxIns) << "\n";
  std::cout << "has_boundary " << (mesh->hasBoundary() ? 1 : 0) << "\n";
  std::cout << "is_delaunay " << (tri.isDelaunay() ? 1 : 0) << "\n";
  std::cout << "min_angle_deg_gc " << tri.minAngleDegrees() << "\n";
  std::cout << "op_seconds " << sec << "\n";
  return 0;
}
