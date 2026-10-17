"""Write a mesh as {"positions", "triangles"} JSON for the Rust helper."""
import json, sys
from common import *
name, out = sys.argv[1], sys.argv[2]
pos, polys, _ = load_obj(model(name))
json.dump({"positions": pos.tolist(), "triangles": triangles(polys).tolist()}, open(out, "w"))
