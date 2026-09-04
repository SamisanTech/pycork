"""Run stock (default) pycork resolveIntersection and write STL."""
import sys
import time
from pathlib import Path

# Prefer the pristine build over the optimized site-packages copy.
sys.path.insert(0, r"E:\github.com\pycork_orig\build\lib.win-amd64-3.9")

import numpy as np
import trimesh
import pycork

STL = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data\21.stl")
OUT = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data\21_default_cork.stl")


def edge_stats(faces):
    e = np.vstack(
        [
            np.sort(faces[:, [0, 1]], 1),
            np.sort(faces[:, [1, 2]], 1),
            np.sort(faces[:, [2, 0]], 1),
        ]
    )
    _, counts = np.unique(e, axis=0, return_counts=True)
    hist = {int(k): int(v) for k, v in zip(*np.unique(counts, return_counts=True))}
    return int((counts == 1).sum()), int((counts >= 3).sum()), hist


def main():
    print(f"pycork file: {pycork.__file__}")
    print("has outerHull:", hasattr(pycork, "outerHull"))
    if hasattr(pycork, "outerHull"):
        raise SystemExit("ERROR: loaded optimized pycork, not default")

    m = trimesh.load_mesh(str(STL), force="mesh", process=False)
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    V = np.ascontiguousarray(m.vertices, dtype=np.float64)
    F = np.ascontiguousarray(m.faces, dtype=np.uint64)
    print(f"input V={len(V):,} F={len(F):,}")

    print("pycork.resolveIntersection (default cork) ...")
    t0 = time.time()
    rv, rf = pycork.resolveIntersection(V, F)
    dt = time.time() - t0
    print(f"done in {dt:.2f}s")

    out = trimesh.Trimesh(vertices=rv, faces=rf, process=False)
    o, nm, hist = edge_stats(rf.astype(np.int64))
    comps = trimesh.graph.connected_components(out.face_adjacency, nodes=np.arange(len(rf)))
    print(
        f"DEFAULT V={len(out.vertices):,} F={len(out.faces):,} open={o} nm={nm} hist={hist} "
        f"shells={len(comps)} watertight={out.is_watertight} winding_ok={out.is_winding_consistent} "
        f"volume={out.volume:.6f}"
    )
    out.export(str(OUT))
    print(f"Wrote {OUT}")
    print(f"TIME {dt:.2f}s")


if __name__ == "__main__":
    main()
