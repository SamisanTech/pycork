"""Stock pycork.resolveIntersection on 24_s.stl. Isolated from site-packages."""
from __future__ import annotations

import collections
import importlib.util
import time
from pathlib import Path

import numpy as np
import trimesh

PYD = Path(r"E:\github.com\pycork_orig\build\lib.win-amd64-3.9\pycork\pycork.cp39-win_amd64.pyd")
STL = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data\24_s.stl")
OUT = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data\24_s_stock_resolve.stl")


def load_stock():
    import sys

    lib = str(PYD.parent.parent)
    sys.path.insert(0, lib)
    import pycork

    print("imported", pycork.__file__, flush=True)
    if "pycork_orig" not in pycork.__file__.replace("\\", "/"):
        raise SystemExit(f"got site-packages pycork: {pycork.__file__}")
    return pycork


def estat(faces):
    f = np.asarray(faces, np.int64)
    e = np.vstack(
        [np.sort(f[:, [0, 1]], 1), np.sort(f[:, [1, 2]], 1), np.sort(f[:, [2, 0]], 1)]
    )
    dt = np.dtype([("a", np.int64), ("b", np.int64)])
    u, c = np.unique(np.ascontiguousarray(e).view(dt).ravel(), return_counts=True)
    return int((c == 1).sum()), int((c >= 3).sum())


def main():
    print("stock pyd", PYD, "exists", PYD.exists(), flush=True)
    stock = load_stock()
    mesh = trimesh.load_mesh(str(STL), force="mesh", process=False)
    mesh.merge_vertices()
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    v = np.ascontiguousarray(mesh.vertices, np.float64)
    f = np.ascontiguousarray(mesh.faces, np.uint64)
    o, n = estat(f)
    print(f"INPUT V={len(v):,} F={len(f):,} open={o} nm={n}", flush=True)
    t0 = time.perf_counter()
    vo, fo = stock.resolveIntersection(v, f)
    dt = time.perf_counter() - t0
    vo = np.ascontiguousarray(vo, np.float64)
    fo = np.ascontiguousarray(fo, np.uint64)
    o, n = estat(fo)
    m = trimesh.Trimesh(vo, fo, process=False)
    comps = trimesh.graph.connected_components(m.face_adjacency, nodes=np.arange(len(fo)))
    print(
        f"STOCK resolve {dt:.3f}s V={len(vo):,} F={len(fo):,} open={o} nm={n} "
        f"shells={len(comps)} wt={m.is_watertight} vol={m.volume:.6f}",
        flush=True,
    )
    rows = []
    for c in comps:
        sub = m.submesh([c], append=True, repair=False)
        rows.append((len(c), float(sub.area), bool(sub.is_watertight)))
    rows.sort(key=lambda r: r[1], reverse=True)
    maxA = rows[0][1] if rows else 1.0
    for i, (nf, area, wt) in enumerate(rows[:15]):
        print(f"  shell F={nf:>9,}  area={area:12.4f}  frac={area/maxA:.6f}  wt={wt}")
    trimesh.Trimesh(vo, fo, process=False).export(str(OUT))
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
