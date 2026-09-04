"""Simple pycork.resolveIntersection + Open3D preview of result only."""
from __future__ import annotations

import argparse
import collections
import time
from pathlib import Path

import numpy as np
import open3d as o3d
import trimesh
import pycork


def edge_stats(faces):
    f = np.asarray(faces, dtype=np.int64)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "stl",
        nargs="?",
        default=r"E:\github.com\SamisanTech\slc_stl\tests\data\21.stl",
    )
    ap.add_argument("--no-preview", action="store_true")
    args = ap.parse_args()

    if int(np.__version__.split(".")[0]) >= 2:
        raise SystemExit("pycork needs numpy<2")

    path = Path(args.stl)
    print(f"Loading {path}")
    verts, faces = pycork.readSTL(str(path))
    verts = np.ascontiguousarray(verts, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    o0, n0, h0 = edge_stats(faces)
    print(f"INPUT  V={len(verts):,} F={len(faces):,} open={o0} nm={n0} hist={h0}")

    # cork has no separate outerHull API — resolveIntersections only remeshes at SI
    print("pycork.resolveIntersection ...")
    t0 = time.time()
    vout, fout = pycork.resolveIntersection(verts, faces)
    dt = time.time() - t0
    print(f"done in {dt:.2f}s")

    out = trimesh.Trimesh(vertices=vout, faces=fout, process=False)
    o1, n1, h1 = edge_stats(fout)
    print(
        f"OUT    V={len(out.vertices):,} F={len(out.faces):,} open={o1} nm={n1} hist={h1} "
        f"solid={pycork.isSolid(np.ascontiguousarray(out.vertices, dtype=np.float64), np.ascontiguousarray(out.faces, dtype=np.uint64))}"
    )

    out_path = path.with_name(path.stem + "_resolve.stl")
    pycork.writeSTL(str(out_path), np.ascontiguousarray(vout, dtype=np.float64),
                    np.ascontiguousarray(fout, dtype=np.uint64))
    print(f"Wrote {out_path}")

    # combined: resolve + outer hull (winding-number classification)
    print("pycork.outerHull (resolve + hull) ...")
    t0 = time.time()
    hv, hf, hstats = pycork.outerHull(verts, faces)
    dth = time.time() - t0
    print(f"done in {dth:.2f}s  stats={hstats}")
    hull = trimesh.Trimesh(vertices=hv, faces=hf, process=False)
    o2, n2, h2 = edge_stats(hf)
    comps = trimesh.graph.connected_components(hull.face_adjacency, nodes=np.arange(len(hf)))
    print(
        f"HULL   V={len(hull.vertices):,} F={len(hull.faces):,} open={o2} nm={n2} hist={h2} "
        f"shells={len(comps)} watertight={hull.is_watertight} winding_ok={hull.is_winding_consistent} "
        f"volume={hull.volume:.6f} (resolved volume={out.volume:.6f}) "
        f"solid={pycork.isSolid(np.ascontiguousarray(hull.vertices, dtype=np.float64), np.ascontiguousarray(hull.faces, dtype=np.uint64))}"
    )
    hull_path = path.with_name(path.stem + "_hull.stl")
    pycork.writeSTL(str(hull_path), np.ascontiguousarray(hv, dtype=np.float64),
                    np.ascontiguousarray(hf, dtype=np.uint64))
    print(f"Wrote {hull_path}")

    if args.no_preview:
        return

    title = f"pycork outerHull (resolve {dt:.2f}s, resolve+hull {dth:.2f}s)"
    print(f"Open3D preview: {title}")
    o3 = o3d.geometry.TriangleMesh()
    o3.vertices = o3d.utility.Vector3dVector(np.asarray(hull.vertices, dtype=np.float64))
    o3.triangles = o3d.utility.Vector3iVector(np.asarray(hull.faces, dtype=np.int32))
    o3.compute_vertex_normals()
    o3d.visualization.draw_geometries([o3], window_name=title)


if __name__ == "__main__":
    main()
