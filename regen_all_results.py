"""Regenerate all cork result STLs for 21.stl."""
import sys
import time
from pathlib import Path

import numpy as np
import trimesh

DATA = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data")
STL = DATA / "21.stl"
ORIG_LIB = r"E:\github.com\pycork_orig\build\lib.win-amd64-3.9"


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


def summarize(label, v, f, dt=None):
    mesh = trimesh.Trimesh(v, f, process=False)
    o, nm, hist = edge_stats(np.asarray(f, dtype=np.int64))
    comps = trimesh.graph.connected_components(
        mesh.face_adjacency, nodes=np.arange(len(f))
    )
    t = f" ({dt:.2f}s)" if dt is not None else ""
    print(
        f"{label}{t}  V={len(mesh.vertices):,} F={len(mesh.faces):,} "
        f"open={o} nm={nm} shells={len(comps)} watertight={mesh.is_watertight} "
        f"winding_ok={mesh.is_winding_consistent} volume={mesh.volume:.6f}"
    )
    return mesh


def unload_pycork():
    for k in list(sys.modules):
        if k == "pycork" or k.startswith("pycork."):
            del sys.modules[k]
    if ORIG_LIB in sys.path:
        sys.path.remove(ORIG_LIB)


def main():
    m = trimesh.load_mesh(str(STL), force="mesh", process=False)
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    V = np.ascontiguousarray(m.vertices, dtype=np.float64)
    F = np.ascontiguousarray(m.faces, dtype=np.uint64)
    print(f"input V={len(V):,} F={len(F):,}\n")

    # ---- optimized build ----
    unload_pycork()
    import pycork as cork_fast

    print(f"[fast] {cork_fast.__file__}")
    assert hasattr(cork_fast, "outerHull")

    print("1/4  optimized resolveIntersection ...")
    t0 = time.time()
    rv, rf = cork_fast.resolveIntersection(V, F)
    dt = time.time() - t0
    summarize("RESOLVE", rv, rf, dt).export(str(DATA / "21_resolve.stl"))
    print(f"  Wrote {DATA / '21_resolve.stl'}\n")

    print("2/4  optimized outerHull (resolve+hull) ...")
    t0 = time.time()
    hv, hf, st = cork_fast.outerHull(V, F, resolve=True)
    dt = time.time() - t0
    print(f"  stats={st}")
    summarize("HULL", hv, hf, dt).export(str(DATA / "21_hull.stl"))
    print(f"  Wrote {DATA / '21_hull.stl'}\n")

    # ---- stock / default cork ----
    unload_pycork()
    sys.path.insert(0, ORIG_LIB)
    import pycork as cork_orig

    print(f"[orig] {cork_orig.__file__}")
    assert not hasattr(cork_orig, "outerHull")

    print("3/4  default cork resolveIntersection ...")
    t0 = time.time()
    drv, drf = cork_orig.resolveIntersection(V, F)
    dt_r = time.time() - t0
    summarize("DEFAULT", drv, drf, dt_r).export(str(DATA / "21_default_cork.stl"))
    print(f"  Wrote {DATA / '21_default_cork.stl'}\n")

    # ---- stock resolve + our hull ----
    unload_pycork()
    import pycork as cork_fast2

    print(f"[fast] {cork_fast2.__file__}")
    assert hasattr(cork_fast2, "outerHull")

    print("4/4  our outerHull(resolve=False) on default cork result ...")
    t0 = time.time()
    ohv, ohf, ost = cork_fast2.outerHull(
        np.ascontiguousarray(drv, dtype=np.float64),
        np.ascontiguousarray(drf, dtype=np.uint64),
        resolve=False,
    )
    dt_h = time.time() - t0
    print(f"  stats={ost}")
    summarize("DEFAULT+OURHULL", ohv, ohf, dt_h).export(
        str(DATA / "21_default_cork_ourhull.stl")
    )
    print(f"  Wrote {DATA / '21_default_cork_ourhull.stl'}")
    print(f"  (default resolve {dt_r:.2f}s + hull {dt_h:.2f}s = {dt_r + dt_h:.2f}s)\n")
    print("ALL DONE")


if __name__ == "__main__":
    main()
