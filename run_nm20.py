"""Default pycork.repair on NM 20.stl → pycork_repair/20_repair.stl"""
from __future__ import annotations

import collections
import time
from pathlib import Path

import numpy as np
import trimesh
import pycork

STL = Path(r"E:\github.com\SamisanTech\manifold\test\Repair Test files\Non Manifold\20.stl")
OUT = STL.parent / "pycork_repair" / "20_repair.stl"


def estat(faces):
    f = np.asarray(faces, np.int64)
    e = np.vstack(
        [np.sort(f[:, [0, 1]], 1), np.sort(f[:, [1, 2]], 1), np.sort(f[:, [2, 0]], 1)]
    )
    dt = np.dtype([("a", np.int64), ("b", np.int64)])
    u, c = np.unique(np.ascontiguousarray(e).view(dt).ravel(), return_counts=True)
    return int((c == 1).sum()), int((c >= 3).sum()), dict(
        sorted(collections.Counter(c.tolist()).items())
    )


def main():
    print("pycork", pycork.__file__, flush=True)
    v, f = pycork.readSTL(str(STL))
    v = np.ascontiguousarray(v, np.float64)
    f = np.ascontiguousarray(f, np.uint64)
    o, n, h = estat(f)
    print(f"IN  V={len(v):,} F={len(f):,} open={o} nm={n} hist={h}", flush=True)
    t0 = time.perf_counter()
    vo, fo, st = pycork.repair(v, f)
    dt = time.perf_counter() - t0
    print(f"repair {dt:.3f}s  {st}", flush=True)
    o, n, h = estat(fo)
    m = trimesh.Trimesh(vo, fo, process=False)
    comps = trimesh.graph.connected_components(m.face_adjacency, nodes=np.arange(len(fo)))
    print(
        f"OUT V={len(vo):,} F={len(fo):,} open={o} nm={n} hist={h} "
        f"shells={len(comps)} wt={m.is_watertight} vol={m.volume:.6f}",
        flush=True,
    )
    OUT.parent.mkdir(exist_ok=True)
    pycork.writeSTL(
        str(OUT),
        np.ascontiguousarray(vo, np.float64),
        np.ascontiguousarray(fo, np.uint64),
    )
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
