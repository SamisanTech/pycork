"""pycork — triangle-soup repair and CSG.

Callable APIs (optional kwargs shown):
  readSTL(path) -> (verts, tris)
  writeSTL(path, verts, tris)
  resolveIntersection(verts, tris) -> (verts, tris)
  outerHull(verts, tris, raysPerPatch=5, resolve=True, hull_backend="ours", exact=False)
      -> (verts, tris, stats)
  repair(verts, tris, resolve=True, hull=True, unify=False, cluster=False,
         puzzle=False, clean=False, perturb=False, noise=True, collapse=False,
         si_subset=False, raysPerPatch=5, minFaces=5, hardDegree=30,
         perturbIntensity=1e-3, collapseRel=0.02, hull_backend="ours",
         exactHull=False)
      -> (verts, tris, stats)
      exactHull: Manifold-style +Z winding in our patch hull (no Manifold lib).
      hull_backend: "ours" (winding, default) or "manifold" / "exact"
      (Exact Manifold Decompose + BatchBoolean Add — not convex Hull).
      Call set_manifold(mod) first; on failure falls back to ours.
  set_manifold(mod) / has_manifold()
  union / difference / intersection
  isSolid(verts, tris)

verts: (N,3) any float dtype — cast to float64 at the CorkMesh boundary
tris:  (M,3) any integer dtype
Kernels stay double (exact predicates). File STL is native float32.
"""
from .pycork import *
