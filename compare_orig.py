"""Run pycork.resolveIntersection with whichever pycork is on sys.path and dump
a full set of result metrics to JSON so two builds can be compared.

usage: python compare_orig.py --tag orig --runs 3 [stl]
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
import time
from pathlib import Path

import numpy as np
import trimesh
import pycork


def edge_stats(faces):
    f = np.asarray(faces, dtype=np.int64)
    edges = np.vstack([np.sort(f[:, [0, 1]], axis=1),
                       np.sort(f[:, [1, 2]], axis=1),
                       np.sort(f[:, [2, 0]], axis=1)])
    dtype = np.dtype([("a", np.int64), ("b", np.int64)])
    structured = np.ascontiguousarray(edges).view(dtype).ravel()
    uniq, counts = np.unique(structured, return_counts=True)
    return {
        "unique_edges": int(len(uniq)),
        "open_edges": int((counts == 1).sum()),
        "nonmanifold_edges": int((counts >= 3).sum()),
        "edge_hist": {str(k): int(v) for k, v in sorted(collections.Counter(counts.tolist()).items())},
    }


def metrics(v, f):
    v = np.asarray(v, dtype=np.float64)
    f = np.asarray(f, dtype=np.int64)
    m = trimesh.Trimesh(vertices=v, faces=f, process=False)
    # connected components by shared edges -> "shells"
    comps = trimesh.graph.connected_components(m.face_adjacency, nodes=np.arange(len(f)))
    comp_sizes = sorted((len(c) for c in comps), reverse=True)
    # duplicate / degenerate faces
    fs = np.sort(f, axis=1)
    dup = len(fs) - len(np.unique(fs, axis=0))
    degen = int((~m.nondegenerate_faces()).sum())
    # unique vertex positions (cork may emit coincident verts)
    uniq_pos = len(np.unique(np.round(v, 9), axis=0))
    out = {
        "V": int(len(v)),
        "F": int(len(f)),
        "unique_positions": int(uniq_pos),
        "volume": float(m.volume),
        "area": float(m.area),
        "bbox_min": [float(x) for x in m.bounds[0]],
        "bbox_max": [float(x) for x in m.bounds[1]],
        "centroid": [float(x) for x in m.centroid],
        "center_mass": [float(x) for x in m.center_mass] if m.is_volume or True else None,
        "shells": len(comp_sizes),
        "shell_sizes_top10": comp_sizes[:10],
        "watertight": bool(m.is_watertight),
        "winding_consistent": bool(m.is_winding_consistent),
        "euler_number": int(m.euler_number),
        "duplicate_faces": int(dup),
        "degenerate_faces": degen,
        "pycork_isSolid": bool(pycork.isSolid(np.ascontiguousarray(v), np.ascontiguousarray(f, dtype=np.uint64))),
    }
    out.update(edge_stats(f))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stl", nargs="?", default=r"E:\github.com\SamisanTech\slc_stl\tests\data\21.stl")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--outdir", default=r"E:\github.com\pycork\compare_out")
    args = ap.parse_args()

    mesh = trimesh.load_mesh(args.stl, force="mesh", process=False)
    mesh.merge_vertices()
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    verts = np.ascontiguousarray(mesh.vertices, dtype=np.float64)
    faces = np.ascontiguousarray(mesh.faces, dtype=np.uint64)

    res = {"tag": args.tag, "module": pycork.__file__, "input": metrics(verts, faces.astype(np.int64)), "runs": []}
    outdir = Path(args.outdir); outdir.mkdir(exist_ok=True)
    for r in range(args.runs):
        t0 = time.perf_counter()
        vo, fo = pycork.resolveIntersection(verts, faces)
        dt = time.perf_counter() - t0
        m = metrics(vo, fo)
        m["seconds"] = dt
        res["runs"].append(m)
        print(f"[{args.tag}] run {r}: {dt:.2f}s V={m['V']} F={m['F']} shells={m['shells']} "
              f"vol={m['volume']:.6f} open={m['open_edges']} nm={m['nonmanifold_edges']}", flush=True)
        if r == 0:
            trimesh.Trimesh(vo, fo, process=False).export(str(outdir / f"{args.tag}_run0.stl"))
    (outdir / f"{args.tag}.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
