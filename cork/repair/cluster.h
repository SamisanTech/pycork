#pragma once

// Face-adjacency shells, AABB clusters, hole-vertex puzzle, pairwise union.
// Mirrors mira_module.non_intersecting_clusters + repair_thread easy/hard.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <vector>

#include <cork/repair/holes.h>
#include <cork/repair/ops.h>
#include <cork/repair/options.h>
#include <cork/util/profile.h>
#include <cork/util/unionFind.h>

namespace cork {
namespace repair {

struct Shell {
    Raw raw;
    BBox3d box;
    Vec3d com;
    double volume = 0;
    double area = 0;
    int openEdges = 0;
    int degree = 0;
    int euler = 0;
    std::vector<uint64_t> holeKeys;
};

inline Raw concat_raw(Raw a, const Raw &b) {
    const int off = (int)a.vertices.size();
    a.vertices.insert(a.vertices.end(), b.vertices.begin(), b.vertices.end());
    a.triangles.reserve(a.triangles.size() + b.triangles.size());
    for (CorkTriangle t : b.triangles) {
        t.a += off;
        t.b += off;
        t.c += off;
        a.triangles.push_back(t);
    }
    return a;
}

inline uint64_t quant_xyz(const Vec3d &p) {
    auto q = [](double x) -> uint32_t {
        return (uint32_t)(int32_t)std::lround(x * 1000.0);
    };
    uint64_t h = q(p.x) * 0x9e3779b1u;
    h ^= q(p.y) + 0x9e3779b1u + (h << 6) + (h >> 2);
    h ^= q(p.z) + 0x9e3779b1u + (h << 6) + (h >> 2);
    return h;
}

inline void measure(Shell &s) {
    const size_t nT = s.raw.triangles.size();
    s.box = BBox3d();
    s.com = Vec3d(0, 0, 0);
    s.volume = 0;
    s.area = 0;
    s.openEdges = 0;
    std::unordered_map<uint64_t, int> ec;
    ec.reserve(nT * 2);
    for (size_t t = 0; t < nT; ++t) {
        const auto &tr = s.raw.triangles[t];
        const Vec3d &a = s.raw.vertices[tr.a].pos;
        const Vec3d &b = s.raw.vertices[tr.b].pos;
        const Vec3d &c = s.raw.vertices[tr.c].pos;
        BBox3d tb(a, a);
        tb = convex(tb, BBox3d(b, b));
        tb = convex(tb, BBox3d(c, c));
        s.box = convex(s.box, tb);
        const Vec3d n = cross(b - a, c - a);
        s.area += 0.5 * len(n);
        s.volume += dot(a, cross(b, c)) / 6.0;
        s.com += (a + b + c) * (1.0 / 3.0);
        auto add = [&](int i, int j) { ec[edge_key(i, j)]++; };
        add(tr.a, tr.b);
        add(tr.b, tr.c);
        add(tr.c, tr.a);
    }
    if (nT) s.com /= (double)nT;
    s.holeKeys.clear();
    for (const auto &kv : ec) {
        if (kv.second == 1) {
            ++s.openEdges;
            const int a = (int)(kv.first >> 32);
            const int b = (int)(kv.first & 0xffffffffu);
            if ((size_t)a < s.raw.vertices.size())
                s.holeKeys.push_back(quant_xyz(s.raw.vertices[a].pos));
            if ((size_t)b < s.raw.vertices.size())
                s.holeKeys.push_back(quant_xyz(s.raw.vertices[b].pos));
        }
    }
    std::sort(s.holeKeys.begin(), s.holeKeys.end());
    s.holeKeys.erase(std::unique(s.holeKeys.begin(), s.holeKeys.end()), s.holeKeys.end());
    const int nE = (int)ec.size();
    s.euler = (int)s.raw.vertices.size() - nE + (int)nT;
}

inline void invert_shell(Shell &s) {
    for (auto &tr : s.raw.triangles) std::swap(tr.b, tr.c);
    s.volume = -s.volume;
}

inline std::vector<Shell> split_shells(const Raw &raw, int minFaces) {
    const size_t nT = raw.triangles.size();
    std::vector<Shell> out;
    if (nT == 0) return out;
    std::unordered_map<uint64_t, int> first;
    first.reserve(nT * 2);
    UnionFind uf((uint)nT);
    for (size_t t = 0; t < nT; ++t) {
        const auto &tr = raw.triangles[t];
        auto add = [&](int a, int b) {
            const uint64_t k = edge_key(a, b);
            auto it = first.find(k);
            if (it == first.end())
                first.emplace(k, (int)t);
            else
                uf.unionIds((uint)it->second, (uint)t);
        };
        add(tr.a, tr.b);
        add(tr.b, tr.c);
        add(tr.c, tr.a);
    }
    std::vector<int> root(nT), count(nT, 0);
    for (size_t t = 0; t < nT; ++t) {
        root[t] = (int)uf.find((uint)t);
        count[root[t]]++;
    }
    std::unordered_map<int, std::vector<int>> faces;
    for (size_t t = 0; t < nT; ++t) {
        if (count[root[t]] < minFaces) continue;
        faces[root[t]].push_back((int)t);
    }
    out.reserve(faces.size());
    for (auto &kv : faces) {
        std::vector<char> keep(nT, 0);
        for (int t : kv.second) keep[t] = 1;
        Shell s;
        s.raw = compact(raw, keep);
        if (s.raw.triangles.empty()) continue;
        measure(s);
        if (s.volume < 0.0 && s.openEdges == 0) invert_shell(s);
        out.push_back(std::move(s));
    }
    return out;
}

inline void aabb_degrees(std::vector<Shell> &shells) {
    const int n = (int)shells.size();
    for (int i = 0; i < n; ++i) shells[i].degree = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (hasIsct(shells[i].box, shells[j].box)) {
                ++shells[i].degree;
                ++shells[j].degree;
            }
}

inline std::vector<std::vector<int>> aabb_clusters(const std::vector<Shell> &shells) {
    const int n = (int)shells.size();
    UnionFind uf((uint)std::max(n, 1));
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (hasIsct(shells[i].box, shells[j].box))
                uf.unionIds((uint)i, (uint)j);
    std::unordered_map<int, std::vector<int>> g;
    for (int i = 0; i < n; ++i) g[(int)uf.find((uint)i)].push_back(i);
    std::vector<std::vector<int>> out;
    out.reserve(g.size());
    for (auto &kv : g) out.push_back(std::move(kv.second));
    return out;
}

inline bool share_hole(const Shell &a, const Shell &b) {
    size_t i = 0, j = 0;
    while (i < a.holeKeys.size() && j < b.holeKeys.size()) {
        if (a.holeKeys[i] == b.holeKeys[j]) return true;
        if (a.holeKeys[i] < b.holeKeys[j])
            ++i;
        else
            ++j;
    }
    return false;
}

// repair_thread puzzle: attach EC 0/1 opens that share hole vertices.
inline int puzzle_opens(std::vector<Shell> &shells) {
    const int n = (int)shells.size();
    if (n < 2) return 0;
    UnionFind uf((uint)n);
    int merges = 0;
    for (int i = 0; i < n; ++i) {
        if (shells[i].openEdges == 0) continue;
        for (int j = i + 1; j < n; ++j) {
            if (shells[j].openEdges == 0) continue;
            const int ecj = shells[j].euler;
            if (ecj != 0 && ecj != 1) continue;
            if (!share_hole(shells[i], shells[j])) continue;
            uf.unionIds((uint)i, (uint)j);
            ++merges;
        }
    }
    if (!merges) return 0;
    std::unordered_map<int, std::vector<int>> g;
    for (int i = 0; i < n; ++i) g[(int)uf.find((uint)i)].push_back(i);
    std::vector<Shell> next;
    next.reserve(g.size());
    int puzzled = 0;
    for (auto &kv : g) {
        auto &idx = kv.second;
        Shell acc = std::move(shells[idx[0]]);
        for (size_t k = 1; k < idx.size(); ++k) {
            acc.raw = concat_raw(std::move(acc.raw), shells[idx[k]].raw);
            ++puzzled;
        }
        measure(acc);
        if (acc.volume < 0.0 && acc.openEdges == 0) invert_shell(acc);
        next.push_back(std::move(acc));
    }
    shells.swap(next);
    return puzzled;
}

inline void com_perturb(Shell &s, const Vec3d &clusterCom, double intensity) {
    Vec3d dir = clusterCom - s.com;
    const double L = len(dir);
    if (L <= 0.0) return;
    dir /= L;
    uint64_t x = (uint64_t)std::rand() ^ 0xA11CEULL;
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x ^= x >> 31;
    const double u = (double)(x >> 11) * (1.0 / 9007199254740992.0);
    const Vec3d d = dir * (intensity * u);
    for (auto &v : s.raw.vertices) v.pos += d;
    measure(s);
}

// Stacked-sheet soups only: SI the AABB-overlapping groups, pass the
// disjoint shells through a single hull.  Full-soup cork was 1.3s of
// 30-sheet hits plus 129 easy parts that never needed SI.
inline void pile_split_repair(CorkMesh &mesh, const Options &opt) {
    CORK_PROF("repair.pile_split");
    std::vector<Shell> shells = split_shells(mesh.raw(), std::max(1, opt.minFaces));
    if (shells.size() <= 1) {
        if (opt.resolve) mesh.resolveIntersections();
        if (opt.hull) {
            cork_hull::HullStats hs;
            mesh.outerHull(opt.raysPerPatch, &hs, opt.noise ? opt.noiseAreaFrac : 0.0);
            if (!hs.closed) fill_boundary_loops(mesh);
        }
        return;
    }
    auto groups = aabb_clusters(shells);
    Raw easy;
    std::vector<std::vector<int>> hard;
    for (auto &g : groups) {
        if (g.size() <= 1) {
            if (easy.triangles.empty())
                easy = std::move(shells[g[0]].raw);
            else
                easy = concat_raw(std::move(easy), shells[g[0]].raw);
        } else {
            hard.push_back(std::move(g));
        }
    }
    cork_prof::note("  pile shells", (double)shells.size());
    cork_prof::note("  pile hard groups", (double)hard.size());

    Raw acc;
    bool have = false;
    auto take = [&](CorkMesh &&m) {
        Raw r = m.raw();
        if (!have) {
            acc = std::move(r);
            have = true;
        } else {
            acc = concat_raw(std::move(acc), r);
        }
    };
    auto finish = [&](CorkMesh &m, bool resolve) {
        if (resolve && opt.resolve) {
            m.preferLbvhIsct = true;
            m.resolveIntersections();
        }
        if (opt.hull) {
            cork_hull::HullStats hs;
            m.outerHull(opt.raysPerPatch, &hs, opt.noise ? opt.noiseAreaFrac : 0.0);
            if (!hs.closed) fill_boundary_loops(m);
        }
    };
    for (auto &g : hard) {
        Raw r;
        for (int i : g) {
            if (r.triangles.empty())
                r = shells[i].raw;
            else
                r = concat_raw(std::move(r), shells[i].raw);
        }
        if (r.triangles.empty()) continue;
        CorkMesh part(std::move(r));
        finish(part, true);
        take(std::move(part));
    }
    if (!easy.triangles.empty()) {
        CorkMesh part(std::move(easy));
        finish(part, false);
        take(std::move(part));
    }
    if (have) adopt(mesh, std::move(acc));
}

inline Raw concat_indices(const std::vector<Shell> &shells, const std::vector<int> &idx) {
    Raw acc;
    for (int i : idx) {
        if (i < 0 || i >= (int)shells.size()) continue;
        if (acc.triangles.empty())
            acc = shells[i].raw;
        else
            acc = concat_raw(std::move(acc), shells[i].raw);
    }
    return acc;
}

inline void pairwise_union(std::vector<CorkMesh> &solids) {
    while (solids.size() > 1) {
        std::vector<CorkMesh> next;
        next.reserve((solids.size() + 1) / 2);
        size_t i = 0;
        if (solids.size() % 2) {
            CorkMesh a = std::move(solids[solids.size() - 2]);
            a.boolUnion(solids[solids.size() - 1]);
            solids[solids.size() - 2] = std::move(a);
            solids.pop_back();
        }
        for (; i + 1 < solids.size(); i += 2) {
            solids[i].boolUnion(solids[i + 1]);
            next.push_back(std::move(solids[i]));
        }
        solids.swap(next);
    }
}

inline void cluster_repair(CorkMesh &mesh, const Options &opt, Stats *st,
                           void (*run_piece)(CorkMesh &, const Options &)) {
    CORK_PROF("repair.cluster");
    std::vector<Shell> shells = split_shells(mesh.raw(), opt.minFaces);
    if (st) st->shells = (int)shells.size();
    {
        const int p = puzzle_opens(shells);
        if (st) st->puzzled = p;
    }
    aabb_degrees(shells);
    auto groups = aabb_clusters(shells);
    if (st) st->clusters = (int)groups.size();

    CorkMesh acc;
    bool have = false;
    for (auto &g : groups) {
        Vec3d com(0, 0, 0);
        double vol = 0;
        int nClosed = 0, nOpen = 0;
        for (int i : g) {
            if (shells[i].openEdges == 0) {
                com += shells[i].com * std::abs(shells[i].volume);
                vol += std::abs(shells[i].volume);
                ++nClosed;
            } else
                ++nOpen;
        }
        if (st) {
            st->closed += nClosed;
            st->open += nOpen;
        }
        if (vol > 0.0) com /= vol;

        std::vector<int> easy, hard, opens;
        for (int i : g) {
            if (opt.perturb && shells[i].openEdges == 0)
                com_perturb(shells[i], com, opt.perturbIntensity);
            if (shells[i].openEdges > 0) {
                opens.push_back(i);
                continue;
            }
            bool isHard = false;
            if (opt.unify && shells[i].degree > opt.hardDegree) {
                int hi = 0;
                for (int j : g)
                    if (j != i && hasIsct(shells[i].box, shells[j].box) &&
                        shells[j].degree > opt.hardDegree)
                        ++hi;
                isHard = hi > opt.hardDegree;
            }
            if (isHard)
                hard.push_back(i);
            else
                easy.push_back(i);
        }
        if (st) st->hard += (int)hard.size();

        Raw pieceRaw;
        if (opt.unify && hard.size() >= 2) {
            std::vector<CorkMesh> solids;
            solids.reserve(hard.size());
            for (int i : hard) solids.emplace_back(shells[i].raw);
            pairwise_union(solids);
            if (!solids.empty()) {
                pieceRaw = solids[0].raw();
                if (st) st->unified += (int)hard.size() - 1;
            }
        } else {
            easy.insert(easy.end(), hard.begin(), hard.end());
        }

        auto take = [&](CorkMesh &&part) {
            if (!have) {
                acc = std::move(part);
                have = true;
            } else {
                acc.disjointUnion(part);
            }
        };
        if (!pieceRaw.triangles.empty()) {
            CorkMesh hardPart(std::move(pieceRaw));
            run_piece(hardPart, opt);
            take(std::move(hardPart));
        }
        if (!easy.empty()) {
            Raw er = concat_indices(shells, easy);
            if (!er.triangles.empty()) {
                CorkMesh easyPart(std::move(er));
                run_piece(easyPart, opt);
                take(std::move(easyPart));
            }
        }
        if (!opens.empty()) {
            Raw oraw = concat_indices(shells, opens);
            if (!oraw.triangles.empty()) {
                CorkMesh openPart(std::move(oraw));
                run_piece(openPart, opt);
                take(std::move(openPart));
            }
        }
    }
    if (have) mesh = std::move(acc);
}

// Mira post_repair_fix on leftover boundary only (hull already closed → skip).
// weld coincident verts, drop exact dups, puzzle hole-sharing opens, fill.
inline int close_leftover(CorkMesh &mesh, const Options &opt, Stats *st) {
    CORK_PROF("repair.close_leftover");
    weld_coincident(mesh);
    drop_noise(mesh, opt.minFaces, 0.0, /*dropTiny=*/false);
    std::vector<Shell> shells = split_shells(mesh.raw(), opt.minFaces);
    if (st) st->shells = (int)shells.size();
    if (shells.size() >= 2) {
        const int p = puzzle_opens(shells);
        if (st) st->puzzled = p;
        Raw acc;
        for (auto &s : shells) {
            if (s.raw.triangles.empty()) continue;
            if (acc.triangles.empty())
                acc = std::move(s.raw);
            else
                acc = concat_raw(std::move(acc), s.raw);
        }
        if (!acc.triangles.empty()) adopt(mesh, std::move(acc));
    }
    return fill_boundary_loops(mesh);
}

}  // namespace repair
}  // namespace cork
