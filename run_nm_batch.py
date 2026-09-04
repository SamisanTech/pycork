"""Batch pycork.repair on the Mira Non Manifold test set. One _repair.stl each."""
from __future__ import annotations

import collections
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import trimesh
import pycork


def edge_stats(faces):
    f = np.asarray(faces, dtype=np.int64)
    if len(f) == 0:
        return 0, 0, {}
    edges = np.vstack(
        [
            np.sort(f[:, [0, 1]], axis=1),
            np.sort(f[:, [1, 2]], axis=1),
            np.sort(f[:, [2, 0]], axis=1),
        ]
    )
    dtype = np.dtype([("a", np.int64), ("b", np.int64)])
    structured = np.ascontiguousarray(edges).view(dtype).ravel()
    uniq, counts = np.unique(structured, return_counts=True)
    open_n = int((counts == 1).sum())
    nm_n = int((counts >= 3).sum())
    hist = dict(sorted(collections.Counter(counts.tolist()).items()))
    return open_n, nm_n, hist


def run_one(path: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / (path.stem + "_repair.stl")
    rec = {
        "file": path.name,
        "in_bytes": path.stat().st_size,
        "ok": False,
    }
    t0 = time.perf_counter()
    verts, faces = pycork.readSTL(str(path))
    rec["io_in"] = round(time.perf_counter() - t0, 3)
    verts = np.ascontiguousarray(verts, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    o0, n0, h0 = edge_stats(faces)
    rec.update(
        in_V=int(len(verts)),
        in_F=int(len(faces)),
        in_open=o0,
        in_nm=n0,
        in_hist=h0,
    )

    t0 = time.perf_counter()
    vout, fout, rstats = pycork.repair(verts, faces)
    rec["repair_s"] = round(time.perf_counter() - t0, 3)
    rec["flags"] = {k: rstats[k] for k in rstats if k in (
        "resolve", "hull", "noise", "unify", "cluster", "puzzle"
    )}

    vout = np.ascontiguousarray(vout, dtype=np.float64)
    fout = np.ascontiguousarray(fout, dtype=np.uint64)
    mesh = trimesh.Trimesh(vertices=vout, faces=fout, process=False)
    o1, n1, h1 = edge_stats(fout)
    if len(fout):
        comps = trimesh.graph.connected_components(
            mesh.face_adjacency, nodes=np.arange(len(fout))
        )
        n_shells = len(comps)
    else:
        n_shells = 0
    rec.update(
        out_V=int(len(vout)),
        out_F=int(len(fout)),
        out_open=o1,
        out_nm=n1,
        out_hist=h1,
        shells=n_shells,
        watertight=bool(mesh.is_watertight) if len(fout) else False,
        winding_ok=bool(mesh.is_winding_consistent) if len(fout) else False,
    )
    try:
        rec["volume"] = float(round(float(mesh.volume), 6))
    except Exception:
        rec["volume"] = None
    try:
        rec["solid"] = bool(pycork.isSolid(vout, fout)) if len(fout) else False
    except Exception:
        rec["solid"] = None

    t0 = time.perf_counter()
    pycork.writeSTL(str(out_path), vout, fout)
    rec["io_out"] = round(time.perf_counter() - t0, 3)
    rec["out_bytes"] = out_path.stat().st_size
    rec["out_path"] = str(out_path)
    rec["ok"] = True
    return rec


def main():
    if len(sys.argv) < 3:
        raise SystemExit("usage: run_nm_batch.py <stl> <out_dir>")
    path = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    if int(np.__version__.split(".")[0]) >= 2:
        raise SystemExit("pycork needs numpy<2")
    try:
        rec = run_one(path, out_dir)
    except Exception as e:
        rec = {
            "file": path.name,
            "ok": False,
            "error": f"{type(e).__name__}: {e}",
            "trace": traceback.format_exc(),
        }
    print(json.dumps(rec, separators=(",", ":")))


if __name__ == "__main__":
    main()
