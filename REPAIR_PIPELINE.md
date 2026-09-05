# Repair pipeline comparison

Three pipelines, same job: take a triangle soup and emit a watertight outer solid.

- **Mira live** — `src/repair_thread.py` (Auto-repair when `NEW_REPAIR` is off)
- **Listed spec** — the cleaned Mira stages (merge → TR → split → puzzle → stitch/SI → 30+ union → noise → Exact → fill)
- **pycork default** — `pycork.repair()` (this repo)

Mira never sends the raw soup through cork. It splits first and only corks subsets. That is why it survives NM piles and why it is slow (Python + trimesh + pyvista + Open3D + decimator + two ohebvh passes).

pycork welds while reading, then one full resolve + our winding hull. On slc `21.stl` that is **~1.05 s repair / ~1.12 s IO+repair**.

Do **not** port Mira ohebvh / ohebvh2 / ohebvh2_box. Our winding hull is the hull.

---

## slc 21.stl gate (must hold)

File: `E:\github.com\SamisanTech\slc_stl\tests\data\21.stl`  
After weld: **V=631,213 F=1,261,664**, open=0, nm=2.

Default `pycork.repair()`:

| step | time |
|---|---:|
| `readSTL` (mmap + exact-float weld) | **0.07 s** |
| eigen → CorkMesh | 0.04 s |
| **resolve** | **0.58 s** |
| **hull** (winding) | **0.41 s** |
| **repair()** | **1.05 s** |
| writeSTL | 0.02 s |
| **IO + repair** | **1.12 s** |

Out: **V=464,990 F=930,606**, open=0, nm=0, **1 shell**, watertight, vol=**484.637424**.

Trimesh on the same file is **0.75 s IO + 1.41 s merge** before any repair. We weld while reading.

---

## Step map

| Listed spec | Mira live (`repair_thread`) | pycork default `repair()` |
|---|---|---|
| Merge vertices, dedup ← replace trimesh | trimesh `merge_vertices` + dup/degen/unref | **done.** `readSTL` exact-float weld, **0.07 s**. No trimesh. |
| TR (0.004 µm, easier puzzle) | pyvista `decimate` on slivers `extent < 0.004`, factor cap 0.8 | **off.** Would change slc 21. `collapse` exists, default off. |
| Split shells → open/closed | `split(only_watertight=False)`, invert if vol < 0, COM-perturb closed | **written, off.** `cluster=True` does face-adj split. Default is one soup. |
| Find connected components (scrap) | AABB `connection_matrix` → clusters; tiny comps dropped later | **written, off.** Same AABB idea. Default never clusters. |
| Puzzle on open (conformal) | hole-vert match, EC 0/1 only, concat + split, HoleFix if ≤20 hole verts | **written, off.** `puzzle=True` (only with `cluster`). No HoleFix. |
| Stitch / SI on leftover open (non-conformal) | optional Delaunay `stitchkaroinko`; then `repair_open_shells`: perturb → **SI-subset cork** → `ohebvh2_box` → fill | **not ported.** No Delaunay stitch. `si_subset` exists, **must stay off** (broke 21: vol 484→416). |
| 30+ multiply-SI shells, batch union (scrap) | AABB degree > 30 → extra TR + pairwise `pycork.union` | **written, off.** `unify=True`. Tried; slower, worse volume. |
| Remove noise | `remove_noise` before cork (“if kept, cork crashes”) | **partial.** Stacked-dup drop only on valence piles (NM 20). slc 21 skips it. |
| SI + ohebvh2×2 + collapse + ohebox | `resolveIntersectionbvh2` (SI faces only + cork) → **ohebvh2 twice** → `graph_collapse` → **box-ray ohe** | **different kernel.** Full `resolveIntersections` + **our winding hull**. No ohebvh. Collapse off. |
| Remove noise (small-area shells) | drop shells `< 0.001 × max area` after box hull | **done, inside hull.** leftover `0.001 × max` when `noise=True`. |
| Fill holes | `mira_repair.post_repair_fix` + HoleFix/TMesh | **ours, only if hull left open.** 3-cycles always; centroid-star if hole verts ≤ 20. Closed slc 21 skips fill. |

---

## What actually runs

```
pycork:  weld 0.07s → resolve 0.58s → winding hull 0.41s → write 0.02s
Mira:    trimesh weld → TR → split → puzzle → stitch/open-SI → hard-union
         → noise → SI-subset cork → ohebvh2 ×2 → collapse → box hull → noise → HoleFix
```

Default `pycork.repair()` flags:

```
clean=False, perturb=False, cluster=False, puzzle=False,
unify=False, collapse=False, si_subset=False
resolve=True, hull=True, noise=True
```

---

## Do not turn on by default

Already measured landmines:

- `si_subset` — AABB-overlap only; destroyed slc 21 (vol 484→416)
- whole-mesh `perturb` / `clean` — opened slc 21
- always-`cluster` — 3.5 s, V/F wobble
- star-cap on holes > 20 verts — NM `6.stl` vol 231→3.7
- drop zero-area faces that share edges with the body — opened slc 21
- ohebvh / ohebvh2 / ohebvh2_box — not our hull

---

## Mira integration

| piece | name |
|---|---|
| submodule | `external/pycork` (this repo) |
| Python module | `pycork` |
| encrypted DLL | `DLLs/unify.dll` |
| load | `main.py` `load_from_bytes("pycork", …/unify.dll)` → `init_repair_mod(…, pycork, …)` |
| C++ repair entry | `pycork.repair()` via `repair_thread.pycork_repair()` / `mode="pycork"` |

Build path (same as other native modules):

1. `scripts/setup.bat` — `cork.dll` for `mira_module` (`BUILD_PYTHON=OFF`), then `setup.py` → `server_files`
2. `scripts/encrypt_pyds.py` — `pycork*.pyd` → `DLLs/unify.dll`
3. `repair_thread` receives the module through `init_repair_mod`, same as holefix / decimator / mira_repair
