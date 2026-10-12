// Stub for bff::BinPacking::pack (the real src/project/BinPacking.cpp needs the
// rectangle-bin-pack submodule, which is not checked out). It is referenced only
// from MeshIO::packUvs, i.e. the MeshIO::write path, which bff_ref never calls.

#include "bff/project/BinPacking.h"
#include <cstdio>
#include <cstdlib>

namespace bff {

void BinPacking::pack(const Model&, double, const std::vector<uint8_t>&,
					  std::vector<Vector>&, std::vector<Vector>&,
					  std::vector<uint8_t>&, Vector&, Vector&)
{
	std::fprintf(stderr, "bff_ref: BinPacking::pack is stubbed out and must not be called\n");
	std::abort();
}

} // namespace bff
