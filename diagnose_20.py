"""NM 20.stl: topology + why stock resolve hangs + noise-then-repair try."""
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


def face_key(f):
    a = np.sort(f, axis=1)
    dt = np.dtype([("a", np.int64), ("b", np.int64), ("c", np.int64)])
    return np.ascontiguousarray(a).view(dt).ravel()


def main():
    print("pycork", pycork.__file__, flush=True)
    t0 = time.perf_counter()
    v, f = pycork.readSTL(str(STL))
    v = np.ascontiguousarray(v, np.float64)
    f = np.ascontiguousarray(f, np.uint64)
    print(f"readSTL {time.perf_counter()-t0:.3f}s  V={len(v):,} F={len(f):,}", flush=True)

    o, n, hist = estat(f)
    print(f"edges open={o} nm={n} hist={hist}", flush=True)

    vf = np.bincount(f.ravel().astype(np.int64), minlength=len(v))
    print(
        f"vert face-valence max={vf.max()} >=10={(vf>=10).sum()} "
        f">=20={(vf>=20).sum()} >=30={(vf>=30).sum()} >=52={(vf>=52).sum()}",
        flush=True,
    )
    tops = np.argsort(-vf)[:8]
    print("  top verts", [(int(i), int(vf[i]), v[i].tolist()) for i in tops], flush=True)

    keys = face_key(f.astype(np.int64))
    uk, kc = np.unique(keys, return_counts=True)
    print(f"exact face dups: extra={int((kc-1).sum())}  groups={(kc>1).sum()}", flush=True)

    m = trimesh.Trimesh(v, f, process=False)
    comps = trimesh.graph.connected_components(m.face_adjacency, nodes=np.arange(len(f)))
    sizes = np.array([len(c) for c in comps])
    print(
        f"shells={len(comps)}  size hist "
        f"1={(sizes==1).sum()} <=5={(sizes<=5).sum()} "
        f">5={(sizes>5).sum()} max={sizes.max()}",
        flush=True,
    )
    order = np.argsort(-sizes)
    for i in order[:10]:
        c = comps[i]
        sub = m.submesh([c], append=True, repair=False)
        print(
            f"  shell F={len(c):>7,} area={sub.area:10.3f} wt={sub.is_watertight} "
            f"vol={sub.volume if sub.is_volume else None}",
            flush=True,
        )

    # drop exact dups + comps <=5  (same as drop_noise minFaces=5)
    keep_uniq = np.zeros(len(f), dtype=bool)
    seen = set()
    for i, k in enumerate(keys.tolist()):
        if k not in seen:
            seen.add(k)
            keep_uniq[i] = True
    uniq_i = np.nonzero(keep_uniq)[0]
    m2 = trimesh.Trimesh(v, f[uniq_i], process=False)
    comps2 = trimesh.graph.connected_components(
        m2.face_adjacency, nodes=np.arange(len(uniq_i))
    )
    keep2 = np.zeros(len(uniq_i), dtype=bool)
    for c in comps2:
        if len(c) > 5:
            keep2[c] = True
    fc = f[uniq_i][keep2]
    # compact unused verts
    used = np.unique(fc)
    remap = np.full(len(v), -1, np.int64)
    remap[used] = np.arange(len(used))
    vc = v[used]
    fc = remap[fc.astype(np.int64)].astype(np.uint64)
    o2, n2, h2 = estat(fc)
    vf2 = np.bincount(fc.ravel().astype(np.int64), minlength=len(vc))
    print(
        f"\nAFTER drop dups+<=5: V={len(vc):,} F={len(fc):,} "
        f"dropped_faces={len(f)-len(fc):,} open={o2} nm={n2} hist={h2}",
        flush=True,
    )
    print(
        f"  valence max={vf2.max()} >=10={(vf2>=10).sum()} "
        f">=30={(vf2>=30).sum()} >=52={(vf2>=52).sum()}",
        flush=True,
    )
    m3 = trimesh.Trimesh(vc, fc, process=False)
    comps3 = trimesh.graph.connected_components(m3.face_adjacency, nodes=np.arange(len(fc)))
    print(f"  shells={len(comps3)}", flush=True)

    # Default repair hangs in resolve (noise skipped when hull=True).
    # Think path: drop_noise then resolve (hull=False), then outerHull.
    print("\n=== repair hull=False (drop_noise then resolve) ===", flush=True)
    t0 = time.perf_counter()
    vr, fr, st = pycork.repair(v, f, hull=False)
    print(f"noise+resolve {time.perf_counter()-t0:.3f}s  {st}", flush=True)
    o3, n3, h3 = estat(fr)
    print(f"  V={len(vr):,} F={len(fr):,} open={o3} nm={n3} hist={h3}", flush=True)

    print("\n=== hull+leftover on resolved (resolve=False) ===", flush=True)
    t0 = time.perf_counter()
    vh, fh, hs = pycork.repair(vr, fr, resolve=False, hull=True)
    print(f"hull {time.perf_counter()-t0:.3f}s  {hs}", flush=True)
    o4, n4, h4 = estat(fh)
    mh = trimesh.Trimesh(vh, fh, process=False)
    ch = trimesh.graph.connected_components(mh.face_adjacency, nodes=np.arange(len(fh)))
    print(
        f"  V={len(vh):,} F={len(fh):,} open={o4} nm={n4} hist={h4} "
        f"shells={len(ch)} wt={mh.is_watertight} vol={mh.volume:.6f}",
        flush=True,
    )
    OUT.parent.mkdir(exist_ok=True)
    pycork.writeSTL(str(OUT), np.ascontiguousarray(vh, np.float64),
                    np.ascontiguousarray(fh, np.uint64))
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
