#pragma once

#include <exception>

#include <cork/cork.h>
#include <cork/repair/cluster.h>
#include <cork/repair/holes.h>
#include <cork/repair/ops.h>
#include <cork/repair/options.h>
#include <cork/util/profile.h>

namespace cork {
namespace repair {

class Pipeline {
public:
    static void run_piece(CorkMesh &mesh, const Options &opt) {
        if (opt.clean) drop_degen_unref(mesh);
        // Mira: drop tiny/dup scraps before cork or it can hang / exit(1).
        // Hull path: only drop exact dups when a busy vertex has an edge
        // used 8+ times (NM 20 valence-52).  slc 21 has busy verts but
        // manifold spokes — skip the FaceKey walk.
        if (opt.noise && !opt.hull) drop_noise(mesh, opt.minFaces, 0.0);
        else if (opt.noise && opt.hull && mesh.hasStackedDuplicateFaces()) {
            drop_noise(mesh, opt.minFaces, 0.0, /*dropTiny=*/false);
            mesh.preferLbvhIsct = true;
        }
        if (opt.perturb) perturb_along_normals(mesh, opt.perturbIntensity);
        if (opt.resolve) {
            CORK_PROF("repair.resolve");
            auto go = [&]() {
                if (opt.si_subset) resolve_si_subset(mesh);
                else mesh.resolveIntersections();
            };
            try {
                go();
            } catch (const std::exception &) {
                // Mira: "if the noise is kept, cork crashes"
                drop_noise(mesh, opt.minFaces);
                go();
            }
        }
        if (opt.collapse) collapse_short_edges(mesh, opt.collapseRel);
        bool needFill = true;
        if (opt.hull && !opt.deferHull) {
            CORK_PROF("repair.hull");
            cork_hull::HullStats hs;
            extract_hull(mesh, opt, &hs);
            needFill = !hs.closed;
        } else if (opt.hull && opt.deferHull) {
            needFill = false;
        } else if (opt.noise) {
            drop_small_shells(mesh, opt.noiseAreaFrac);
        }
        // close_leftover (weld+puzzle+fill) is opt.puzzle — leftover split
        // on a 100k+ soup is seconds and did not close these rims.
        if (needFill && opt.fill) {
            if (opt.puzzle) close_leftover(mesh, opt, nullptr);
            else fill_boundary_loops(mesh);
        }
    }

    static void run(CorkMesh &mesh, const Options &opt, Stats *st = nullptr) {
        if (opt.cluster || opt.unify) {
            cluster_repair(mesh, opt, st, run_piece);
            return;
        }
        run_piece(mesh, opt);
    }
};

}  // namespace repair
}  // namespace cork
