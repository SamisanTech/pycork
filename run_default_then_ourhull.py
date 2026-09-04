"""Stock cork resolveIntersection -> our optimized outerHull (resolve=False)."""
import importlib
import sys
import time
from pathlib import Path

import numpy as np
import trimesh

STL = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data\21.stl")
OUT_RESOLVE = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data\21_default_cork.stl")
OUT_HULL = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data\21_default_cork_ourhull.stl")

ORIG_LIB = r"E:\github.com\pycork_orig\build\lib.win-amd64-3.9"
FAST_PYD = Path(
    r"C:\Users\Karan\AppData\Local\Programs\Python\Python39\lib\site-packages"
    r"\pycork-0.1.3-py3.9-win-amd64.egg\pycork\pycork.cp39-win_amd64.pyd"
)


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


def load_ext(name, pyd_or_pkg_path, as_package=False):
    if as_package:
        sys.path.insert(0, pyd_or_pkg_path)
        if name in sys.modules:
            del sys.modules[name]
        for k in list(sys.modules):
            if k.startswith(name + "."):
                del sys.modules[k]
        return importlib.import_module(name)
    # load .pyd as extension under a unique module name, then wrap
    # Prefer importing the egg package after clearing orig path
    raise RuntimeError("use package path")


def main():
    m = trimesh.load_mesh(str(STL), force="mesh", process=False)
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    V = np.ascontiguousarray(m.vertices, dtype=np.float64)
    F = np.ascontiguousarray(m.faces, dtype=np.uint64)
    print(f"input V={len(V):,} F={len(F):,}")

    # 1) stock cork resolve
    print("--- stock cork resolveIntersection ---")
    sys.path.insert(0, ORIG_LIB)
    import pycork as cork_orig

    print("orig:", cork_orig.__file__, "outerHull?", hasattr(cork_orig, "outerHull"))
    t0 = time.time()
    rv, rf = cork_orig.resolveIntersection(V, F)
    dt_r = time.time() - t0
    print(f"resolve done in {dt_r:.2f}s  V={len(rv):,} F={len(rf):,}")
    trimesh.Trimesh(rv, rf, process=False).export(str(OUT_RESOLVE))
    print(f"Wrote {OUT_RESOLVE}")

    # 2) unload orig, load our optimized pycork (with outerHull)
    print("--- our outerHull(resolve=False) on stock result ---")
    del sys.modules["pycork"]
    for k in list(sys.modules):
        if k.startswith("pycork."):
            del sys.modules[k]
    if ORIG_LIB in sys.path:
        sys.path.remove(ORIG_LIB)

    # site-packages egg is already on path
    import pycork as cork_fast

    print("fast:", cork_fast.__file__, "outerHull?", hasattr(cork_fast, "outerHull"))
    if not hasattr(cork_fast, "outerHull"):
        raise SystemExit("optimized pycork without outerHull")

    rv64 = np.ascontiguousarray(rv, dtype=np.float64)
    rf64 = np.ascontiguousarray(rf, dtype=np.uint64)
    t1 = time.time()
    hv, hf, st = cork_fast.outerHull(rv64, rf64, resolve=False)
    dt_h = time.time() - t1
    print(f"outerHull done in {dt_h:.2f}s  stats={st}")

    hull = trimesh.Trimesh(hv, hf, process=False)
    o, nm, hist = edge_stats(hf.astype(np.int64))
    comps = trimesh.graph.connected_components(hull.face_adjacency, nodes=np.arange(len(hf)))
    print(
        f"HULL V={len(hull.vertices):,} F={len(hull.faces):,} open={o} nm={nm} hist={hist} "
        f"shells={len(comps)} watertight={hull.is_watertight} winding_ok={hull.is_winding_consistent} "
        f"volume={hull.volume:.6f}"
    )
    hull.export(str(OUT_HULL))
    print(f"Wrote {OUT_HULL}")
    print(f"TOTAL resolve {dt_r:.2f}s + hull {dt_h:.2f}s = {dt_r + dt_h:.2f}s")


if __name__ == "__main__":
    main()
