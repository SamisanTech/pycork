"""AABB neighbor degree of NM 20 shells after drop exact dups."""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import trimesh
import pycork

STL = Path(r"E:\github.com\SamisanTech\manifold\test\Repair Test files\Non Manifold\20.stl")


def face_key(f):
    a = np.sort(f, axis=1)
    dt = np.dtype([("a", np.int64), ("b", np.int64), ("c", np.int64)])
    return np.ascontiguousarray(a).view(dt).ravel()


def main():
    v, f = pycork.readSTL(str(STL))
    keys = face_key(f.astype(np.int64))
    keep = np.zeros(len(f), dtype=bool)
    seen = set()
    for i, k in enumerate(keys.tolist()):
        if k not in seen:
            seen.add(k)
            keep[i] = True
    f = f[keep]
    m = trimesh.Trimesh(v, f, process=False)
    comps = trimesh.graph.connected_components(m.face_adjacency, nodes=np.arange(len(f)))
    print(f"unique-face shells={len(comps)}", flush=True)
    boxes = []
    for c in comps:
        pts = v[f[c].ravel().astype(np.int64)]
        boxes.append((pts.min(0), pts.max(0), len(c)))
    n = len(boxes)
    deg = np.zeros(n, dtype=np.int32)
    t0 = time.perf_counter()
    for i in range(n):
        a0, a1, _ = boxes[i]
        for j in range(i + 1, n):
            b0, b1, _ = boxes[j]
            if np.all(a1 >= b0) and np.all(b1 >= a0):
                deg[i] += 1
                deg[j] += 1
    print(f"AABB pairs {time.perf_counter()-t0:.3f}s", flush=True)
    print(
        f"degree max={deg.max()} >=10={(deg>=10).sum()} "
        f">=20={(deg>=20).sum()} >=30={(deg>=30).sum()} "
        f">=50={(deg>=50).sum()}",
        flush=True,
    )
    order = np.argsort(-deg)[:12]
    for i in order:
        print(
            f"  shell F={boxes[i][2]:>6,}  aabb_neighbors={deg[i]}  "
            f"box={boxes[i][0].tolist()} .. {boxes[i][1].tolist()}",
            flush=True,
        )
    # how many shells sit in the single densest pile (neighbors of the max-deg shell + itself)
    hot = int(np.argmax(deg))
    pile = {hot}
    a0, a1, _ = boxes[hot]
    for j in range(n):
        b0, b1, _ = boxes[j]
        if np.all(a1 >= b0) and np.all(b1 >= a0):
            pile.add(j)
    print(f"densest pile size={len(pile)} faces={sum(boxes[i][2] for i in pile)}", flush=True)


if __name__ == "__main__":
    main()
