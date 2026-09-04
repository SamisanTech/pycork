"""Profiled resolveIntersection only."""
import os, time
import numpy as np
import trimesh
import pycork

path = r"E:\github.com\SamisanTech\slc_stl\tests\data\21.stl"
m = trimesh.load_mesh(path, force="mesh", process=False)
m.merge_vertices()
m.update_faces(m.nondegenerate_faces())
m.remove_unreferenced_vertices()
V = np.ascontiguousarray(m.vertices, dtype=np.float64)
F = np.ascontiguousarray(m.faces, dtype=np.uint64)
print(f"input V={len(V):,} F={len(F):,}")

# warmup + one timed run
if os.environ.get("CORK_WARMUP"):
    pycork.resolveIntersection(V, F)

t0 = time.time()
v, f = pycork.resolveIntersection(V, F)
dt = time.time() - t0
print(f"done in {dt:.3f}s  V={len(v):,} F={len(f):,}")
