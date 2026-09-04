"""Run pycork.repair (one integrated pass) and write a single STL."""
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
    out_path = path.with_name(path.stem + "_repair.stl")
    print(f"Loading {path}")
    print(f"pycork {pycork.__file__}")

    t0 = time.perf_counter()
    verts, faces = pycork.readSTL(str(path))
    io_in = time.perf_counter() - t0
    verts = np.ascontiguousarray(verts, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    o0, n0, h0 = edge_stats(faces)
    print(f"INPUT  V={len(verts):,} F={len(faces):,} open={o0} nm={n0} hist={h0}")
    print(f"IO in  {io_in:.3f}s")

    print("pycork.repair (resolve + our hull) ...")
    t0 = time.perf_counter()
    vout, fout, rstats = pycork.repair(verts, faces)
    dt = time.perf_counter() - t0
    print(f"repair {dt:.3f}s  {rstats}")

    mesh = trimesh.Trimesh(vertices=vout, faces=fout, process=False)
    o1, n1, h1 = edge_stats(fout)
    comps = trimesh.graph.connected_components(mesh.face_adjacency, nodes=np.arange(len(fout)))
    print(
        f"OUT    V={len(mesh.vertices):,} F={len(mesh.faces):,} open={o1} nm={n1} hist={h1} "
        f"shells={len(comps)} watertight={mesh.is_watertight} winding_ok={mesh.is_winding_consistent} "
        f"volume={mesh.volume:.6f} "
        f"solid={pycork.isSolid(np.ascontiguousarray(mesh.vertices, dtype=np.float64), np.ascontiguousarray(mesh.faces, dtype=np.uint64))}"
    )

    t0 = time.perf_counter()
    pycork.writeSTL(
        str(out_path),
        np.ascontiguousarray(vout, dtype=np.float64),
        np.ascontiguousarray(fout, dtype=np.uint64),
    )
    io_out = time.perf_counter() - t0
    print(f"Wrote {out_path}")
    print(f"IO out {io_out:.3f}s  IO+repair {io_in + dt:.3f}s")

    if args.no_preview:
        return

    title = f"pycork repair {dt:.2f}s"
    print(f"Open3D preview: {title}")
    o3 = o3d.geometry.TriangleMesh()
    o3.vertices = o3d.utility.Vector3dVector(np.asarray(mesh.vertices, dtype=np.float64))
    o3.triangles = o3d.utility.Vector3iVector(np.asarray(mesh.faces, dtype=np.int32))
    o3.compute_vertex_normals()
    o3d.visualization.draw_geometries([o3], window_name=title)


if __name__ == "__main__":
    main()
