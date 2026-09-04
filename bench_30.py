"""30.stl: our default path, then each experimental flag."""
import collections
import os
import sys
import time

import numpy as np
import trimesh
import pycork

STL = sys.argv[1] if len(sys.argv) > 1 else r"E:\github.com\SamisanTech\slc_stl\tests\data\30.stl"


def edge_stats(faces):
    f = np.asarray(faces, dtype=np.int64)
    e = np.vstack([np.sort(f[:, [0, 1]], 1), np.sort(f[:, [1, 2]], 1), np.sort(f[:, [2, 0]], 1)])
    _, counts = np.unique(e, axis=0, return_counts=True)
    hist = {int(k): int(v) for k, v in zip(*np.unique(counts, return_counts=True))}
    return int((counts == 1).sum()), int((counts >= 3).sum()), hist


def summarize(label, v, f, dt):
    m = trimesh.Trimesh(v, f, process=False)
    o, nm, hist = edge_stats(f)
    comps = trimesh.graph.connected_components(m.face_adjacency, nodes=np.arange(len(f)))
    print(
        f"{label}  {dt:.3f}s  V={len(m.vertices):,} F={len(m.faces):,} "
        f"open={o} nm={nm} shells={len(comps)} watertight={m.is_watertight} "
        f"winding={m.is_winding_consistent} vol={m.volume:.6f} hist={hist}"
    )


def load():
    mesh = trimesh.load_mesh(STL, force="mesh", process=False)
    mesh.merge_vertices()
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    V = np.ascontiguousarray(mesh.vertices, dtype=np.float64)
    F = np.ascontiguousarray(mesh.faces, dtype=np.uint64)
    return V, F


def clear_flags():
    for k in ("CORK_LBVH", "CORK_LBVH_TRI", "CORK_HULL_FLOOD", "CORK_AUTOPOLICY"):
        os.environ.pop(k, None)


def run(label, resolve=True, hull=True):
    V, F = load()
    if resolve:
        t0 = time.time()
        rv, rf = pycork.resolveIntersection(V, F)
        summarize(f"{label} RESOLVE", rv, rf, time.time() - t0)
    if hull:
        t0 = time.time()
        hv, hf, st = pycork.outerHull(V, F, resolve=True)
        print(f"  hull stats {st}")
        summarize(f"{label} HULL   ", hv, hf, time.time() - t0)


def main():
    V, F = load()
    o, nm, hist = edge_stats(F)
    print(f"FILE {STL}")
    print(f"INPUT V={len(V):,} F={len(F):,} open={o} nm={nm} hist={hist}")
    print(f"pycork {pycork.__file__}\n")

    print("======== OUR CODE (no flags) ========")
    clear_flags()
    run("OURS")

    print("\n======== CORK_LBVH=1 ========")
    clear_flags()
    os.environ["CORK_LBVH"] = "1"
    run("LBVH")

    print("\n======== CORK_LBVH=1 CORK_LBVH_TRI=1 ========")
    clear_flags()
    os.environ["CORK_LBVH"] = "1"
    os.environ["CORK_LBVH_TRI"] = "1"
    run("LBVH_TRI", hull=False)

    print("\n======== CORK_HULL_FLOOD=1 ========")
    clear_flags()
    os.environ["CORK_HULL_FLOOD"] = "1"
    run("FLOOD")

    print("\n======== CORK_AUTOPOLICY=1 ========")
    clear_flags()
    os.environ["CORK_AUTOPOLICY"] = "1"
    run("AUTOPOL")


if __name__ == "__main__":
    main()
