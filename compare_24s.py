"""Compare 24_s.stl: input shells, our resolve, our hull leftover=0, our repair."""
from __future__ import annotations

import collections
import time
from pathlib import Path

import numpy as np
import trimesh
import pycork

STL = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data\24_s.stl")
REPAIR = Path(r"E:\github.com\SamisanTech\slc_stl\tests\data\24_s_repair.stl")


def estat(faces):
    f = np.asarray(faces, np.int64)
    e = np.vstack(
        [np.sort(f[:, [0, 1]], 1), np.sort(f[:, [1, 2]], 1), np.sort(f[:, [2, 0]], 1)]
    )
    dt = np.dtype([("a", np.int64), ("b", np.int64)])
    u, c = np.unique(np.ascontiguousarray(e).view(dt).ravel(), return_counts=True)
    return int((c == 1).sum()), int((c >= 3).sum())


def shells(v, f):
    m = trimesh.Trimesh(v, f, process=False)
    comps = trimesh.graph.connected_components(m.face_adjacency, nodes=np.arange(len(f)))
    rows = []
    for i, c in enumerate(comps):
        sub = m.submesh([c], append=True, repair=False)
        rows.append(
            {
                "i": i,
                "F": int(len(c)),
                "area": float(sub.area),
                "vol": float(sub.volume) if sub.is_volume else None,
                "wt": bool(sub.is_watertight),
            }
        )
    rows.sort(key=lambda r: r["area"], reverse=True)
    return rows


def show(name, v, f, dt=None):
    o, n = estat(f)
    sh = shells(v, f)
    m = trimesh.Trimesh(v, f, process=False)
    print(
        f"\n=== {name} ===  V={len(v):,} F={len(f):,} open={o} nm={n} "
        f"shells={len(sh)} wt={m.is_watertight} vol={m.volume:.6f}"
        + (f"  {dt:.3f}s" if dt is not None else "")
    )
    maxA = sh[0]["area"] if sh else 1.0
    for r in sh[:12]:
        print(
            f"  shell F={r['F']:>9,}  area={r['area']:12.4f}  "
            f"frac={r['area']/maxA:.6f}  wt={r['wt']}  vol={r['vol']}"
        )
    if len(sh) > 12:
        print(f"  ... {len(sh) - 12} more")
    return sh


def main():
    print("ours", pycork.__file__)
    v, f = pycork.readSTL(str(STL))
    v = np.ascontiguousarray(v, np.float64)
    f = np.ascontiguousarray(f, np.uint64)
    show("INPUT readSTL", v, f)

    if REPAIR.exists():
        vr, fr = pycork.readSTL(str(REPAIR))
        show("OUR repair STL (re-weld)", vr, fr)

    print("\nOUR resolveIntersection ...", flush=True)
    t0 = time.perf_counter()
    vo, fo = pycork.resolveIntersection(v, f)
    dt = time.perf_counter() - t0
    show("OUR resolve", vo, fo, dt)

    print("\nOUR outerHull leftover=0 (stock hull prune only) ...", flush=True)
    t0 = time.perf_counter()
    vh, fh, st = pycork.outerHull(vo, fo, raysPerPatch=5, resolve=False)
    dt = time.perf_counter() - t0
    print("hull stats", dict(st))
    show("OUR hull leftover=0", vh, fh, dt)


if __name__ == "__main__":
    main()
