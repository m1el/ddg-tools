// bff_ref: thin command-line driver around the reference C++ implementation of
// Boundary First Flattening, for differential testing.
//
//   bff_ref MODE in.obj out.txt [cones.txt]
//
//   MODE = cones : cones.txt holds lines "index angle" (1-based OBJ vertex index,
//                  angle in radians; blank lines and '#' comments allowed).
//                  Mirrors flatten() in apps/command-line/src/CommandLine.cpp
//                  after ConePlacement, with the cones taken from the file:
//                  BFF bff(mesh); Cutter::cut(cones, mesh); bff.flattenWithCones(C, true)
//   MODE = auto  : DenseMatrix u(bff.data->bN); bff.flatten(u, true)
//   MODE = disk  : bff.flattenToDisk()
//
// No HoleFiller, no Generators, no projectUvsToPcaAxis, no UV packing.
//
// Output (all floats "%.17g"), everything in ORIGINAL 0-based OBJ vertex indices
// and ORIGINAL OBJ face order:
//
//   # bff_ref MODE
//   faces N
//   f v0 v1 v2  u0 v0  u1 v1  u2 v2
//   cut M
//   c a b
//   positions V
//   p x y z
//
// Exit status: 0 ok, 1 usage / input rejected / reference failure,
//              3 output written but its UVs are non-finite or all exactly zero
//                (BFF::normalize() overwrites NaN UVs with zero).

#include "bff/project/Bff.h"
#include "bff/mesh/MeshIO.h"
#include "bff/project/Cutter.h"

#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace bff;

static void die(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	std::fprintf(stderr, "bff_ref: error: ");
	std::vfprintf(stderr, fmt, ap);
	std::fprintf(stderr, "\n");
	va_end(ap);
	std::exit(1);
}

static double now()
{
	using namespace std::chrono;
	return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// Scans the OBJ the same way MeshIO::readOBJ tokenizes it and rejects anything
// that is not a pure triangle mesh (readOBJ would silently fan-triangulate
// polygons, which changes the face count and order).
static void scanObj(const std::string& path, int& nVLines, int& nFLines)
{
	std::ifstream in(path.c_str());
	if (!in.is_open()) die("unable to open '%s'", path.c_str());

	nVLines = 0;
	nFLines = 0;
	std::string line;
	int lineNo = 0;
	while (std::getline(in, line)) {
		lineNo++;
		std::istringstream ss(line);
		std::string token;
		ss >> token;

		if (token == "v") {
			nVLines++;

		} else if (token == "f") {
			int count = 0;
			while (ss >> token) count++;
			if (count != 3) {
				die("'%s' line %d: face has %d vertices; only pure triangle meshes are "
					"supported (MeshIO::read would triangulate polygons and change the "
					"face list)", path.c_str(), lineNo, count);
			}
			nFLines++;
		}
	}
}

struct ConeSpec {
	int orig; // 0-based original OBJ vertex index
	double angle;
	int lineNo;
};

static std::vector<ConeSpec> readCones(const std::string& path)
{
	std::ifstream in(path.c_str());
	if (!in.is_open()) die("unable to open cone file '%s'", path.c_str());

	std::vector<ConeSpec> cones;
	std::string line;
	int lineNo = 0;
	while (std::getline(in, line)) {
		lineNo++;
		size_t hash = line.find('#');
		if (hash != std::string::npos) line = line.substr(0, hash);

		std::istringstream ss(line);
		std::string a, b, extra;
		if (!(ss >> a)) continue; // blank
		if (!(ss >> b)) die("'%s' line %d: expected 'index angle'", path.c_str(), lineNo);
		if (ss >> extra) die("'%s' line %d: trailing garbage '%s'", path.c_str(), lineNo, extra.c_str());

		char *end = NULL;
		long idx = std::strtol(a.c_str(), &end, 10);
		if (*end != '\0') die("'%s' line %d: bad vertex index '%s'", path.c_str(), lineNo, a.c_str());
		double angle = std::strtod(b.c_str(), &end);
		if (*end != '\0' || !std::isfinite(angle)) {
			die("'%s' line %d: bad angle '%s'", path.c_str(), lineNo, b.c_str());
		}

		ConeSpec c;
		c.orig = (int)idx - 1;
		c.angle = angle;
		c.lineNo = lineNo;
		cones.push_back(c);
	}

	return cones;
}

int main(int argc, const char *argv[])
{
	if (argc < 4 || argc > 5) {
		std::fprintf(stderr, "USAGE: %s MODE in.obj out.txt [cones.txt]\n"
							 "  MODE = cones | auto | disk  (cones.txt required for 'cones')\n",
					 argv[0]);
		return 1;
	}

	const std::string mode = argv[1];
	const std::string inPath = argv[2];
	const std::string outPath = argv[3];
	const bool modeCones = mode == "cones";
	const bool modeAuto = mode == "auto";
	const bool modeDisk = mode == "disk";
	if (!modeCones && !modeAuto && !modeDisk) die("unknown MODE '%s' (cones|auto|disk)", mode.c_str());
	if (modeCones && argc != 5) die("MODE cones needs a cones.txt argument");
	if (!modeCones && argc != 4) die("MODE %s takes no cones.txt argument", mode.c_str());

	double t0 = now();

	// ---- pre-scan and raw polygon soup (the reference's own OBJ parser) --------------
	int nVLines = 0, nFLines = 0;
	scanObj(inPath, nVLines, nFLines);

	std::string error;
	PolygonSoup raw;
	std::vector<std::pair<int, int>> uncuttableEdges;
	if (!MeshIO::readOBJ(inPath, raw, uncuttableEdges, error)) {
		die("MeshIO::readOBJ failed on '%s': %s", inPath.c_str(), error.c_str());
	}

	const int nOrig = (int)raw.positions.size();
	const int nFaces = (int)raw.indices.size()/3;
	if (nOrig != nVLines) die("internal: readOBJ gave %d positions for %d 'v' lines", nOrig, nVLines);
	if ((int)raw.indices.size() != 3*nFLines || !uncuttableEdges.empty()) {
		die("internal: readOBJ gave %d indices for %d 'f' lines", (int)raw.indices.size(), nFLines);
	}
	if (nFaces == 0 || nOrig == 0) die("'%s' has no faces or no vertices", inPath.c_str());

	for (int I = 0; I < 3*nFaces; I += 3) {
		for (int J = 0; J < 3; J++) {
			int i = raw.indices[I + J];
			if (i < 0 || i >= nOrig) {
				die("face %d references vertex %d, out of range [1, %d]", I/3, i + 1, nOrig);
			}
		}
		if (raw.indices[I] == raw.indices[I + 1] || raw.indices[I + 1] == raw.indices[I + 2] ||
			raw.indices[I + 2] == raw.indices[I]) {
			die("face %d has a repeated vertex", I/3);
		}
	}

	// MeshIO::buildModel first calls PolygonSoup::removeIsolatedVertices: if (and only
	// if) some vertex is unreferenced, vertices are renumbered in order of first
	// appearance in the face list. Replicate that map so everything can be reported in
	// original indices (the replication is checked against positions below).
	std::vector<int> origToSoup(nOrig, -1);
	std::vector<int> soupToOrig;
	{
		std::vector<uint8_t> referenced(nOrig, 0);
		for (size_t k = 0; k < raw.indices.size(); k++) referenced[raw.indices[k]] = 1;
		int nReferenced = 0;
		for (int i = 0; i < nOrig; i++) nReferenced += referenced[i];

		if (nReferenced == nOrig) {
			soupToOrig.resize(nOrig);
			for (int i = 0; i < nOrig; i++) {
				soupToOrig[i] = i;
				origToSoup[i] = i;
			}

		} else {
			for (size_t k = 0; k < raw.indices.size(); k++) {
				int i = raw.indices[k];
				if (origToSoup[i] == -1) {
					origToSoup[i] = (int)soupToOrig.size();
					soupToOrig.push_back(i);
				}
			}
			std::fprintf(stderr, "bff_ref: warning: %d of %d OBJ vertices are not referenced by "
								 "any face; the reference renumbers the rest in order of first "
								 "appearance (output is mapped back to original indices)\n",
						 nOrig - nReferenced, nOrig);
		}
	}
	const int nUsed = (int)soupToOrig.size();

	// ---- load through the reference -------------------------------------------------
	Model model;
	if (!MeshIO::read(inPath, model, error)) {
		die("MeshIO::read failed on '%s': %s", inPath.c_str(), error.c_str());
	}
	if (model.nVertices() != nUsed) {
		// PolygonSoup::splitNonManifoldVertices appended duplicates (this is also how
		// the reference deals with non-manifold edges); those have no original index
		die("the reference built %d vertices from the %d the OBJ references, i.e. it split "
			"non-manifold vertices/edges (into %d component(s)); such meshes are not "
			"supported by this driver", model.nVertices(), nUsed, model.size());
	}
	if (model.size() != 1) {
		die("mesh has %d connected components; exactly 1 is required", model.size());
	}

	Mesh& mesh = model[0];
	if ((int)mesh.vertices.size() != nUsed) {
		die("internal: reference mesh has %d vertices, expected %d", (int)mesh.vertices.size(), nUsed);
	}
	if ((int)mesh.faces.size() != nFaces) {
		die("internal: reference mesh has %d faces, OBJ has %d", (int)mesh.faces.size(), nFaces);
	}

	// mesh vertex -> original OBJ index, via the model's own map
	std::vector<int> vertexToOrig(nUsed);
	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		int local = (int)(v - mesh.vertices.begin());
		if (v->index != local) die("internal: Vertex::index %d at position %d", v->index, local);
		int g = model.globalVertexIndex(0, v->index);
		if (g < 0 || g >= nUsed) die("internal: globalVertexIndex out of range");
		std::pair<int, int> l = model.localVertexIndex(g);
		if (l.first != 0 || l.second != v->index) die("internal: model index maps are not inverse");
		vertexToOrig[v->index] = soupToOrig[g];

		// check the map against the geometry: undo MeshIO::normalize
		Vector p = v->position*mesh.radius + mesh.cm;
		const Vector& q = raw.positions[soupToOrig[g]];
		double tol = 1e-9*(1.0 + mesh.radius + mesh.cm.norm());
		if ((p - q).norm() > tol) {
			die("internal: vertex %d does not map back to OBJ vertex %d", v->index, soupToOrig[g] + 1);
		}
	}

	// face order and corner order: compare every reference face with the OBJ face
	std::vector<HalfEdgeCIter> faceStart(nFaces);
	int nExact = 0, nRotated = 0, nFlipped = 0;
	for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		int fi = (int)(f - mesh.faces.begin());
		if (f->index != fi) die("internal: Face::index %d at position %d", f->index, fi);

		const int *o = &raw.indices[3*fi];
		HalfEdgeCIter h0 = f->halfEdge()->next(); // first halfedge created for this face
		HalfEdgeCIter h = h0;
		int n = 0;
		int t[3] = {-1, -1, -1};
		do {
			if (n < 3) t[n] = vertexToOrig[h->vertex()->index];
			n++;
			h = h->next();
		} while (h != h0);
		if (n != 3) die("internal: face %d has %d sides", fi, n);

		// start at the corner of the OBJ's first vertex
		int r = -1;
		for (int J = 0; J < 3; J++) if (t[J] == o[0]) r = J;
		if (r == -1) die("internal: face %d (%d %d %d) is not OBJ face %d (%d %d %d)",
						 fi, t[0], t[1], t[2], fi, o[0], o[1], o[2]);
		h = h0;
		for (int J = 0; J < r; J++) h = h->next();
		faceStart[fi] = h;

		int a = t[(r + 1) % 3], b = t[(r + 2) % 3];
		if (a == o[1] && b == o[2]) {
			if (r == 0) nExact++; else nRotated++;
		} else if (a == o[2] && b == o[1]) {
			nFlipped++;
		} else {
			die("internal: face %d (%d %d %d) is not OBJ face %d (%d %d %d)",
				fi, t[0], t[1], t[2], fi, o[0], o[1], o[2]);
		}
	}
	if (nFlipped > 0) {
		std::fprintf(stderr, "bff_ref: warning: the reference reversed the orientation of %d of "
							 "%d faces (inconsistently oriented OBJ); those 'f' lines list "
							 "v0 v2 v1 of the OBJ face\n", nFlipped, nFaces);
	}

	// ---- topology requirements -------------------------------------------------------
	const int nBoundaries = (int)mesh.boundaries.size();
	const int chi = mesh.eulerCharacteristic();
	if (nBoundaries > 1) {
		die("mesh has %d boundary loops; the reference needs HoleFiller for that, which this "
			"driver does not run (0 or 1 loops for 'cones', exactly 1 for 'auto'/'disk')",
			nBoundaries);
	}
	if (chi + nBoundaries != 2) {
		die("mesh is not genus 0 (V - E + F = %d with %d boundary loops); the reference needs "
			"Generators for that, which this driver does not run", chi, nBoundaries);
	}
	if (!modeCones && nBoundaries != 1) {
		die("MODE %s needs a disk (exactly one boundary loop); this mesh is closed", mode.c_str());
	}

	double t1 = now();

	// ---- flatten ---------------------------------------------------------------------
	BFF bff(mesh);
	double t2 = now(), t3 = t2;
	int nCones = 0;

	if (modeCones) {
		std::vector<ConeSpec> specs = readCones(argv[4]);
		std::vector<VertexIter> cones;
		DenseMatrix C(bff.data->iN);
		std::vector<uint8_t> isCone(nUsed, 0);

		for (size_t k = 0; k < specs.size(); k++) {
			const ConeSpec& s = specs[k];
			if (s.orig < 0 || s.orig >= nOrig) {
				die("cone file line %d: vertex index %d out of range [1, %d]", s.lineNo, s.orig + 1, nOrig);
			}
			if (origToSoup[s.orig] == -1) {
				die("cone file line %d: vertex %d is not referenced by any face", s.lineNo, s.orig + 1);
			}

			std::pair<int, int> l = model.localVertexIndex(origToSoup[s.orig]);
			VertexIter v = mesh.vertices.begin() + l.second;
			if (isCone[v->index]) die("cone file line %d: vertex %d listed twice", s.lineNo, s.orig + 1);
			if (v->onBoundary()) die("cone file line %d: vertex %d is on the boundary", s.lineNo, s.orig + 1);

			isCone[v->index] = 1;
			C(bff.data->index[v->wedge()]) = s.angle;
			cones.push_back(v);
		}
		nCones = (int)cones.size();

		Cutter::cut(cones, mesh);
		t3 = now();

		int nCut = 0;
		for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) nCut += e->onCut ? 1 : 0;
		if (nBoundaries == 0 && nCut == 0) {
			die("closed mesh and Cutter::cut produced no cut (%d cones given; at least 2 are "
				"needed to cut a closed surface)", nCones);
		}

		if (!bff.flattenWithCones(C, true)) {
			die("BFF::flattenWithCones failed (status = %s)",
				bff.status == BFF::ErrorCode::factorizationFailed ? "factorizationFailed" : "ok");
		}

	} else if (modeAuto) {
		DenseMatrix u(bff.data->bN);
		if (!bff.flatten(u, true)) {
			die("BFF::flatten failed (status = %s)",
				bff.status == BFF::ErrorCode::factorizationFailed ? "factorizationFailed" : "ok");
		}

	} else {
		if (!bff.flattenToDisk()) {
			die("BFF::flattenToDisk failed (status = %s)",
				bff.status == BFF::ErrorCode::factorizationFailed ? "factorizationFailed" : "ok");
		}
	}

	double t4 = now();

	// ---- write -----------------------------------------------------------------------
	FILE *out = std::fopen(outPath.c_str(), "w");
	if (!out) die("unable to open '%s' for writing", outPath.c_str());

	long nNonFinite = 0;
	long nZero = 0; // BFF::normalize() overwrites NaN UVs with (0, 0)
	std::fprintf(out, "# bff_ref %s\n", mode.c_str());
	std::fprintf(out, "faces %d\n", nFaces);
	for (int fi = 0; fi < nFaces; fi++) {
		// The corner at he->vertex() is he->next()->corner() (Corner::vertex() is
		// halfEdge()->prev()->vertex()); this is how Distortion.cpp and
		// MeshIO::collectModelUvs read UVs.
		HalfEdgeCIter h = faceStart[fi];
		int idx[3];
		double uv[3][2];
		for (int J = 0; J < 3; J++) {
			CornerCIter c = h->next()->corner();
			if (c->vertex()->index != h->vertex()->index) die("internal: corner/vertex mismatch in face %d", fi);
			idx[J] = vertexToOrig[h->vertex()->index];
			uv[J][0] = c->uv.x;
			uv[J][1] = c->uv.y;
			if (!std::isfinite(uv[J][0]) || !std::isfinite(uv[J][1])) nNonFinite++;
			if (uv[J][0] == 0.0 && uv[J][1] == 0.0) nZero++;
			h = h->next();
		}

		std::fprintf(out, "f %d %d %d  %.17g %.17g  %.17g %.17g  %.17g %.17g\n",
					 idx[0], idx[1], idx[2],
					 uv[0][0], uv[0][1], uv[1][0], uv[1][1], uv[2][0], uv[2][1]);
	}

	int nCut = 0;
	for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) nCut += e->onCut ? 1 : 0;
	std::fprintf(out, "cut %d\n", nCut);
	for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
		if (e->onCut) {
			HalfEdgeCIter h = e->halfEdge();
			std::fprintf(out, "c %d %d\n", vertexToOrig[h->vertex()->index],
						 vertexToOrig[h->flip()->vertex()->index]);
		}
	}

	// The wedge at which the reference starts to walk the (cut) boundary: the original
	// vertex and the face of the corner. The boundary curve is laid out from there.
	{
		CutPtrSet cb = mesh.cutBoundary();
		CutPtrIterator it = cb.begin();
		if (it != cb.end()) {
			WedgeIter w = *it;
			std::fprintf(out, "start %d %d\n", vertexToOrig[w->vertex()->index], w->face()->index);
		}
	}

	std::fprintf(out, "positions %d\n", nOrig);
	for (int i = 0; i < nOrig; i++) {
		Vector p;
		if (origToSoup[i] != -1) {
			p = mesh.vertices[model.localVertexIndex(origToSoup[i]).second].position;

		} else {
			// unreferenced vertex: not part of the reference mesh at all; apply the
			// same transform MeshIO::normalize applied to the others
			p = raw.positions[i];
			p -= mesh.cm;
			p /= mesh.radius;
		}

		std::fprintf(out, "p %.17g %.17g %.17g\n", p.x, p.y, p.z);
	}

	if (std::fclose(out) != 0) die("error writing '%s'", outPath.c_str());
	double t5 = now();

	std::fprintf(stderr, "bff_ref: %s: V=%d (OBJ %d) E=%d F=%d boundaryLoops=%d cones=%d cutEdges=%d\n",
				 mode.c_str(), nUsed, nOrig, (int)mesh.edges.size(), nFaces, nBoundaries, nCones, nCut);
	std::fprintf(stderr, "bff_ref: corner order vs OBJ: exact=%d rotated=%d flipped=%d\n",
				 nExact, nRotated, nFlipped);
	std::fprintf(stderr, "bff_ref: normalize: cm = (%.17g, %.17g, %.17g) radius = %.17g\n",
				 mesh.cm.x, mesh.cm.y, mesh.cm.z, mesh.radius);
	std::fprintf(stderr, "bff_ref: time: load %.3fs, BFF setup %.3fs, cut %.3fs, flatten %.3fs, "
						 "write %.3fs, total %.3fs\n",
				 t1 - t0, t2 - t1, t3 - t2, t4 - t3, t5 - t4, t5 - t0);

	if (nNonFinite > 0) {
		std::fprintf(stderr, "bff_ref: error: %ld corners with non-finite UVs in output\n", nNonFinite);
		return 3;
	}
	if (nZero == 3L*nFaces) {
		// the reference replaces NaN UVs by (0, 0) in BFF::normalize(), and one NaN
		// poisons the centre of mass it subtracts from every corner
		std::fprintf(stderr, "bff_ref: error: every UV is exactly (0, 0): the reference produced "
							 "NaNs and zeroed them\n");
		return 3;
	}
	if (nZero > 0) {
		std::fprintf(stderr, "bff_ref: warning: %ld corners have UV exactly (0, 0)\n", nZero);
	}

	return 0;
}
