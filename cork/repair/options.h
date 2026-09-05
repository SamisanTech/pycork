#pragma once

namespace cork {
namespace repair {

// Optional stages. Defaults match today's 21.stl resolve-only kernel.
// Extra Mira stages stay off unless the caller asks.
struct Options {
    bool resolve = true;
    bool hull = true;         // Mira Exact path always extracts outer hull
    bool exactHull = false;   // false = our rays; true = Manifold +Z winding
    bool unify = false;       // pairwise cork union of hard closed shells
    bool cluster = false;     // split-first AABB groups (opt-in; slc 21 gate)
    bool puzzle = false;      // hole-vert puzzle (on in cluster_repair)
    bool si_subset = false;   // intersecting-face cork (opt-in; slc 21 gate)
    bool clean = false;       // degen drop — opens slc 21.stl if always on
    bool perturb = false;     // Mira jitter — opens slc 21.stl if always on
    bool noise = true;        // pre-SI dups+tiny; leftover shells in hull
    bool collapse = false;    // short-edge weld (mira_module.collapse_edges)
    bool fill = true;         // boundary fill after hull; off when Manifold hull follows
    bool deferHull = false;   // hull=true noise path, but skip our outerHull (Manifold follows)
    int raysPerPatch = 5;
    int minFaces = 5;         // noise: Mira post_repair_fix / remove_noise
    int hardDegree = 30;      // Mira SIcount: AABB-neighbor degree for "hard"
    double perturbIntensity = 1e-3;
    double collapseRel = 0.02; // vs mean edge length (repair_thread.graph_collapse)
    double noiseAreaFrac = 0.001; // after hull: drop shells < this * max area
};

struct Stats {
    int shells = 0;
    int clusters = 0;
    int closed = 0;
    int open = 0;
    int hard = 0;
    int unified = 0;
    int puzzled = 0;
    int noiseDropped = 0;
};

}  // namespace repair
}  // namespace cork
