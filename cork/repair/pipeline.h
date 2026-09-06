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
        // Compact stacked pile (NM 26: 2.9M in a 22mm AABB). Full resolve
        // shreds to ~76M faces (~289s) and leaves a 23k-open belt. Hull
        // the dup-dropped soup instead (LBVH). NM 20 is 328k — resolves.
        // slc 21 never preferLbvhIsct.
        const bool compactPile =
            mesh.preferLbvhIsct && mesh.numTris() > 1000000;
        if (opt.resolve && !compactPile) {
            CORK_PROF("repair.resolve");
            auto go = [&]() {
                if (opt.si_subset) resolve_si_subset(mesh);
                else mesh.resolveIntersections();
            };
            try {
                go();
            } catch (const std::exception &) {
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
            // 21 is one closed shell — the FaceKey/UF walk is a ~0.6s no-op.
            // Pile path (20/26) still drops leftover dust.
            if (opt.noise && mesh.preferLbvhIsct)
                drop_small_shells(mesh, opt.noiseAreaFrac);
            needFill = !hs.closed;
        } else if (opt.hull && opt.deferHull) {
            needFill = false;
        } else if (opt.noise) {
            drop_small_shells(mesh, opt.noiseAreaFrac);
        }
        if (needFill && opt.fill) {
            if (opt.puzzle) close_leftover(mesh, opt, nullptr);
            else fill_boundary_loops(mesh, compactPile);
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
