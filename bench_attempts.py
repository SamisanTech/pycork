"""Time resolve and resolve+hull. Prints one line per phase label."""
import os, sys, time
import numpy as np
import trimesh
import pycork

path = os.environ.get("CORK_STL", r"E:\github.com\SamisanTech\slc_stl\tests\data\21.stl")
if len(sys.argv) > 1:
    path = sys.argv[1]
m = trimesh.load_mesh(path, force="mesh", process=False)
m.merge_vertices()
m.update_faces(m.nondegenerate_faces())
m.remove_unreferenced_vertices()
V = np.ascontiguousarray(m.vertices, dtype=np.float64)
F = np.ascontiguousarray(m.faces, dtype=np.uint64)
label = os.environ.get("CORK_ATTEMPT", "baseline")
mode = os.environ.get("CORK_BENCH", "both")  # resolve|hull|both

if mode in ("resolve", "both"):
    t0 = time.time()
    v, f = pycork.resolveIntersection(V, F)
    dt = time.time() - t0
    print(f"ATTEMPT {label} resolve {dt:.3f}s V={len(v)} F={len(f)}")
    sys.stdout.flush()

if mode in ("hull", "both"):
    t0 = time.time()
    hv, hf, st = pycork.outerHull(V, F, resolve=True)
    dt = time.time() - t0
    print(
        f"ATTEMPT {label} hull {dt:.3f}s V={len(hv)} F={len(hf)} "
        f"unresolved={st.get('unresolved_patches')} pruned={st.get('sliver_faces_pruned')}"
    )
    sys.stdout.flush()
