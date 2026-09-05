// +-------------------------------------------------------------------------
// | mesh.hull.tpp
// |
// | Outer hull of a single, already intersection-resolved mesh.
// |
// | Cork's BoolProblem classifies faces of *two* labelled operands by parity
// | and cannot handle a self-overlapping mesh (regions with winding number 2
// | would be treated as "outside").  This module instead uses generalized
// | winding numbers:
// |
// |   1. faces are grouped into patches connected through manifold edges
// |      (edges with exactly two incident faces); intersection curves and
// |      other non-manifold edges separate patches
// |   2. seed faces per patch (more for larger patches).  Rays both ways
// |      through the SI uniform grid on a normal soup, or our LBVH when
// |      the arrangement is a sheet pile (NM edges / faces > 0.15).
// |      Crossing count is w+ / w-; require w- == w+ + 1.
// |   3. a face is kept if its normal side is exterior (w+ == 0), kept and
// |      flipped if its back side is exterior (w- == 0), otherwise deleted.
// |      If the seeds of a patch disagree (this happens where cork leaves a
// |      gap in an intersection curve so both sides of the curve fell into
// |      one patch) the patch is split by a multi-source BFS from the seeds
// |      and every face takes the class of its nearest seed.
// |
// | Only vertices referenced by surviving faces are kept.
// +-------------------------------------------------------------------------
#pragma once

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <mutex>
#include <iostream>
#include <unordered_map>

namespace cork_hull {

struct HullStats {
    size_t patches = 0, split_patches = 0, seeds = 0, kept = 0, flipped = 0,
           deleted = 0, rays = 0, unresolved = 0, pruned = 0;
    bool closed = true; // no open/NM boundary left after prune
};

// splitmix64 -- deterministic per (patch, sample) random stream
inline uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
inline double unit01(uint64_t x) { return (double)(x >> 11) * (1.0 / 9007199254740992.0); }

// Moller-Trumbore.  On hit returns +1 if the ray leaves the solid through
// this face (d . n > 0) and -1 if it enters; 0 for a miss.
inline int rayTriSign(const Vec3d &p, const Vec3d &d,
                      const Vec3d &a, const Vec3d &b, const Vec3d &c, double &t,
                      double tMin)
{
    Vec3d e1 = b - a, e2 = c - a;
    Vec3d pv = cross(d, e2);
    double det = dot(e1, pv);           // = -(d . n)
    if (std::fabs(det) < 1e-300) return 0;
    double inv = 1.0 / det;
    Vec3d tv = p - a;
    double u = dot(tv, pv) * inv;
    if (u < 0.0 || u > 1.0) return 0;
    Vec3d qv = cross(tv, e1);
    double v = dot(d, qv) * inv;
    if (v < 0.0 || u + v > 1.0) return 0;
    t = dot(e2, qv) * inv;
    // faces exactly coincident with the source face (duplicated sheets) sit
    // at t ~ 0 and are not crossings of the ray leaving the source point
    if (t <= tMin) return 0;
    return (det < 0.0) ? +1 : -1;
}

enum : int8_t { HULL_UNSET = -1, HULL_KEEP = 0, HULL_FLIP = 1, HULL_DELETE = 2 };

} // namespace cork_hull

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::hasStackedDuplicateFaces() const
{
    const size_t nT = tris.size();
    const size_t nV = verts.size();
    if (nT == 0 || nV == 0) return false;
    std::vector<int> deg(nV, 0);
    for (size_t t = 0; t < nT; ++t) {
        const Tri &tr = tris[t];
        if (tr.a == tr.b || tr.b == tr.c || tr.c == tr.a) continue;
        deg[tr.a]++; deg[tr.b]++; deg[tr.c]++;
    }
    // Stacked duplicate sheets make extreme fans (NM 20: 234, NM 26: 160).
    // slc 21's busiest vertex is a feature fan (~96) with manifold spokes
    // and emax=4 — stop before any edge walk.
    int vmax = 0;
    for (int d : deg) if (d > vmax) vmax = d;
    if (vmax < 100) return false;

    // Any undirected edge used 8+ times. NM 26 has 111 such edges, none
    // of them incident to the valence>=100 verts, so a hot-vert-only
    // scan missed the pile and Triangle later exit(1)'d.
    auto ek = [](int a, int b) -> uint64_t {
        uint32_t u = (uint32_t)std::min(a, b), w = (uint32_t)std::max(a, b);
        return ((uint64_t)u << 32) | w;
    };
    std::unordered_map<uint64_t, int> ec;
    ec.reserve(nT * 2);
    for (size_t t = 0; t < nT; ++t) {
        const Tri &tr = tris[t];
        if (tr.a == tr.b || tr.b == tr.c || tr.c == tr.a) continue;
        const int vs[3] = { (int)tr.a, (int)tr.b, (int)tr.c };
        for (int k = 0; k < 3; ++k)
            if (++ec[ek(vs[k], vs[(k + 1) % 3])] >= 8)
                return true;
    }
    return false;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::outerHull(int raysPerPatch, cork_hull::HullStats *stats,
                                       double leftoverAreaFrac, bool exact)
{
    using namespace cork_si;
    using namespace cork_hull;
    CORK_PROF("outerHull total");
    const size_t nt = tris.size();
    if (nt == 0) return;
    const int maxTries = std::max(3, raysPerPatch);
    const bool dbg = std::getenv("CORK_HULL_DEBUG") != nullptr;
    const bool flood = std::getenv("CORK_HULL_FLOOD") != nullptr;

    // ------------------------------------------------------------------
    // 1. patches + face adjacency across manifold edges
    // ------------------------------------------------------------------
    std::vector<uint32_t> patchOf(nt);
    RawArray<uint32_t> nb(3 * nt);          // neighbour across edge k (or ~0u)
    RawArray<uint32_t> nmEdgeOf(3 * nt);    // non-manifold edge id of slot k (or ~0u)
    std::vector<uint32_t> nmStart(1, 0), nmFace;   // CSR: nm edge -> incident faces
    uint32_t np = 0;
    {
        CORK_PROF("  hull: patches");
        struct EK { uint64_t key; uint32_t tri; uint32_t k; };
        RawArray<EK> ek(3 * nt);
        cork_par::for_each_idx(nt, 8192, [&](size_t i) {
            const Tri &t = tris[i];
            for (uint k = 0; k < 3; ++k) {
                uint a = t.v[(k + 1) % 3], b = t.v[(k + 2) % 3];
                uint lo = std::min(a, b), hi = std::max(a, b);
                ek[3 * i + k].key = ((uint64_t)lo << 32) | (uint64_t)hi;
                ek[3 * i + k].tri = (uint32_t)i;
                ek[3 * i + k].k   = k;
                nb[3 * i + k] = ~0u;
                nmEdgeOf[3 * i + k] = ~0u;
            }
        });
        cork_par::sort(ek.begin(), ek.end(), [](const EK &x, const EK &y) {
            return x.key < y.key || (x.key == y.key && x.tri < y.tri);
        });
        std::vector<uint32_t> parent(nt);
        for (size_t i = 0; i < nt; ++i) parent[i] = (uint32_t)i;
        auto find = [&](uint32_t x) {
            while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
            return x;
        };
        const size_t n3 = ek.size();
        for (size_t i = 0; i < n3;) {
            size_t j = i + 1;
            while (j < n3 && ek[j].key == ek[i].key) ++j;
            if (j - i == 2) {                       // manifold edge -> same patch
                uint32_t fa = ek[i].tri, fb = ek[i + 1].tri;
                nb[3 * fa + ek[i].k]     = fb;
                nb[3 * fb + ek[i + 1].k] = fa;
                uint32_t a = find(fa), b = find(fb);
                if (a != b) parent[std::max(a, b)] = std::min(a, b);
            } else {                                // open or non-manifold edge
                uint32_t eid = (uint32_t)(nmStart.size() - 1);
                for (size_t k = i; k < j; ++k) {
                    nmEdgeOf[3 * ek[k].tri + ek[k].k] = eid;
                    nmFace.push_back(ek[k].tri);
                }
                nmStart.push_back((uint32_t)nmFace.size());
            }
            i = j;
        }
        std::vector<uint32_t> pid(nt, ~0u);
        for (size_t i = 0; i < nt; ++i) {
            uint32_t r = find((uint32_t)i);
            if (pid[r] == ~0u) pid[r] = np++;
            patchOf[i] = pid[r];
        }
    }
    // patch -> faces (CSR)
    std::vector<uint32_t> pstart(np + 1, 0), plist(nt);
    {
        for (size_t i = 0; i < nt; ++i) pstart[patchOf[i] + 1]++;
        for (size_t p = 0; p < np; ++p) pstart[p + 1] += pstart[p];
        std::vector<uint32_t> fill(pstart.begin(), pstart.end() - 1);
        for (size_t i = 0; i < nt; ++i) plist[fill[patchOf[i]]++] = (uint32_t)i;
    }

    // ------------------------------------------------------------------
    // 2. Spatial index for winding rays.
    //    Normal soup (slc 21): uniform-grid DDA — cache-friendly on a
    //    near-2-manifold.  Sheet pile (NM 20): LBVH — grid cells are
    //    packed 30-deep and DDA hangs.  Switch on NM-edge density
    //    after patches (free): pile if nNm/nFaces > 0.15 or nNm > 200k.
    // ------------------------------------------------------------------
    const size_t nNm = nmStart.size() > 0 ? nmStart.size() - 1 : 0;
    const bool pile = nNm > 200000 || (nt > 0 && nNm * 20 > nt * 3);
    RawArray<BBox3d> tb(nt);
    CellGrid grid;
    RawArray<uint64_t> ent;
    RawArray<uint32_t> cellStart;
    cork_lbvh::LBVH tree;
    double cellH = 1.0;
    double tMin = 0.0;
    {
        CORK_PROF(pile ? "  hull: lbvh" : "  hull: grid");
        struct Acc { BBox3d box; double sum = 0.0; };
        cork_par::Local<Acc> acc;
        cork_par::for_range(nt, 8192, [&](size_t b0, size_t b1) {
            Acc &a = acc.local();
            for (size_t i = b0; i < b1; ++i) {
                BBox3d bb;
                for (uint k = 0; k < 3; ++k) {
                    const Vec3d &p = verts[tris[i].v[k]].pos;
                    bb.minp = min(bb.minp, p); bb.maxp = max(bb.maxp, p);
                }
                tb[i] = bb;
                a.box = convex(a.box, bb);
                Vec3d d = bb.maxp - bb.minp;
                a.sum += std::max(d.x, std::max(d.y, d.z));
            }
        });
        BBox3d world; double meanExt = 0.0;
        acc.combine_each([&](const Acc &a) { world = convex(world, a.box); meanExt += a.sum; });
        meanExt /= (double)nt;
        Vec3d ext = world.maxp - world.minp;
        double maxExt = std::max(ext.x, std::max(ext.y, ext.z));
        tMin = 1e-13 * std::max(maxExt, 1.0);
        cork_prof::note("  hull: #patches", (double)np);
        cork_prof::note("  hull: #nm edges", (double)nNm);
        if (pile) {
            tree.build(tb.data(), nt, world);
            cork_prof::note("  hull: #lbvh leaves", (double)tree.n);
        } else {
            double h = 3.0 * meanExt;
            if (!(h > 0.0)) h = (maxExt > 0.0) ? maxExt : 1.0;
            for (;;) {
                long long nx = std::max(1LL, (long long)std::ceil(ext.x / h) + 1);
                long long ny = std::max(1LL, (long long)std::ceil(ext.y / h) + 1);
                long long nz = std::max(1LL, (long long)std::ceil(ext.z / h) + 1);
                if (nx * ny * nz <= (1LL << 24)) {
                    grid.nx = (int)nx; grid.ny = (int)ny; grid.nz = (int)nz;
                    break;
                }
                h *= 1.25;
            }
            cellH = h;
            double pad = 1e-9 * std::max(maxExt, 1.0);
            grid.org = world.minp - Vec3d(pad, pad, pad);
            grid.inv = Vec3d(1.0 / h, 1.0 / h, 1.0 / h);
            const size_t ncells = (size_t)grid.nx * grid.ny * grid.nz;
            buildEntries(grid, tb, ent);
            const size_t ne = ent.size();
            cellStart.alloc(ncells + 1);
            cork_par::for_each_idx(ncells + 1, 65536, [&](size_t c) { cellStart[c] = ~0u; });
            cork_par::for_each_idx(ne, 16384, [&](size_t i) {
                if (i == 0 || entryKey(ent[i]) != entryKey(ent[i - 1]))
                    cellStart[entryKey(ent[i])] = (uint32_t)i;
            });
            cellStart[ncells] = (uint32_t)ne;
            for (size_t c = ncells; c-- > 0;)
                if (cellStart[c] == ~0u) cellStart[c] = cellStart[c + 1];
            cork_prof::note("  hull: #grid entries", (double)ne);
        }
    }

    // ------------------------------------------------------------------
    // 3. winding numbers by ray casting
    // ------------------------------------------------------------------
    auto finishHits = [&](std::vector<std::pair<uint32_t,int>> &hits, int &w) -> bool {
        std::sort(hits.begin(), hits.end());
        w = 0;
        for (size_t i = 0; i < hits.size();) {
            size_t j = i + 1;
            while (j < hits.size() && hits[j].first == hits[i].first) ++j;
            for (size_t k = i + 1; k < j; ++k) if (hits[k].second != hits[i].second) return false;
            w += hits[i].second;
            i = j;
        }
        return true;
    };
    auto castWinding = [&](const Vec3d &p, const Vec3d &d, uint32_t skip,
                           std::vector<std::pair<uint32_t,int>> &hits, int &w) -> bool
    {
        hits.clear();
        if (pile) {
            const Vec3d inv(
                (std::fabs(d.x) > 1e-300) ? 1.0 / d.x : (d.x >= 0.0 ? 1e300 : -1e300),
                (std::fabs(d.y) > 1e-300) ? 1.0 / d.y : (d.y >= 0.0 ? 1e300 : -1e300),
                (std::fabs(d.z) > 1e-300) ? 1.0 / d.z : (d.z >= 0.0 ? 1e300 : -1e300));
            tree.ray(p, inv, [&](uint32_t tid) {
                if (tid == skip) return;
                const Tri &t = tris[tid];
                double tt;
                int s = rayTriSign(p, d, verts[t.a].pos, verts[t.b].pos, verts[t.c].pos, tt, tMin);
                if (s != 0) hits.push_back({tid, s});
            });
        } else {
            int ix, iy, iz;
            grid.cellOf(p, ix, iy, iz);
            int    step[3];
            double tMax[3], tDelta[3];
            const double dv[3] = { d.x, d.y, d.z };
            const int    ic[3] = { ix, iy, iz };
            const int    nn[3] = { grid.nx, grid.ny, grid.nz };
            const double pv[3] = { p.x, p.y, p.z };
            const double ov[3] = { grid.org.x, grid.org.y, grid.org.z };
            for (int k = 0; k < 3; ++k) {
                if (dv[k] > 0)      { step[k] = 1;  tMax[k] = (ov[k] + (ic[k] + 1) * cellH - pv[k]) / dv[k]; tDelta[k] = cellH / dv[k]; }
                else if (dv[k] < 0) { step[k] = -1; tMax[k] = (ov[k] + ic[k] * cellH - pv[k]) / dv[k];       tDelta[k] = -cellH / dv[k]; }
                else                { step[k] = 0;  tMax[k] = INFINITY; tDelta[k] = INFINITY; }
            }
            int cur[3] = { ix, iy, iz };
            for (;;) {
                uint32_t key = grid.key(cur[0], cur[1], cur[2]);
                for (uint32_t e = cellStart[key]; e < cellStart[key + 1]; ++e) {
                    uint32_t tid = entryIdx(ent[e]);
                    if (tid == skip) continue;
                    const Tri &t = tris[tid];
                    double tt;
                    int s = rayTriSign(p, d, verts[t.a].pos, verts[t.b].pos, verts[t.c].pos, tt, tMin);
                    if (s != 0) hits.push_back({tid, s});
                }
                int ax = (tMax[0] < tMax[1]) ? ((tMax[0] < tMax[2]) ? 0 : 2) : ((tMax[1] < tMax[2]) ? 1 : 2);
                cur[ax] += step[ax];
                if (cur[ax] < 0 || cur[ax] >= nn[ax]) break;
                tMax[ax] += tDelta[ax];
            }
        }
        return finishHits(hits, w);
    };

    // classify one face by its own winding numbers; HULL_UNSET if no
    // consistent ray pair could be found
    auto classifyFace = [&](uint32_t f, uint64_t salt,
                            std::vector<std::pair<uint32_t,int>> &hits, size_t &rays) -> int8_t
    {
        const Tri &t = tris[f];
        const Vec3d &a = verts[t.a].pos, &b = verts[t.b].pos, &c = verts[t.c].pos;
        Vec3d n = cross(b - a, c - a);
        double nl = len(n);
        if (!(nl > 0.0)) return HULL_UNSET;
        n = n / nl;
        Vec3d ctr = (a + b + c) / 3.0;
        for (int s = 0; s < maxTries; ++s) {
            if (exact && s >= 3) break;
            Vec3d d;
            if (exact && s < 3) {
                // Manifold PointWinding is +Z signed crossings.  +X/+Y are
                // the same axis-aligned fallback when +Z grazes the face.
                d = Vec3d(s == 2 ? 1.0 : 0.0, s == 1 ? 1.0 : 0.0, s == 0 ? 1.0 : 0.0);
            } else if (pile && s < 3) {
                d = Vec3d(s == 0 ? 1.0 : 0.0, s == 1 ? 1.0 : 0.0, s == 2 ? 1.0 : 0.0);
            } else {
                uint64_t h0 = mix64(salt ^ ((uint64_t)f << 8) ^ (uint64_t)s);
                uint64_t h1 = mix64(h0), h2 = mix64(h1);
                d = Vec3d(unit01(h0) * 2.0 - 1.0, unit01(h1) * 2.0 - 1.0, unit01(h2) * 2.0 - 1.0);
                double dl = len(d);
                if (!(dl > 1e-3)) continue;
                d = d / dl;
            }
            double cs = dot(d, n);
            if (cs < 0.0) { d = -d; cs = -cs; }
            if (cs < 0.15) continue;                 // avoid grazing the face itself
            int wp, wm;
            rays += 2;
            if (!castWinding(ctr,  d, f, hits, wp)) continue;
            if (!castWinding(ctr, -d, f, hits, wm)) continue;
            if (wm != wp + 1) continue;               // grazed an edge somewhere
            if (wp == 0) return HULL_KEEP;
            if (wm == 0) return HULL_FLIP;
            return HULL_DELETE;
        }
        return HULL_UNSET;
    };

    std::vector<int8_t> faceCode(nt, HULL_UNSET);
    size_t nSplit = 0, nSeeds = 0, nRays = 0, nUnresolved = 0;

    // seeds: every face of a small patch, otherwise ~1 per 50 faces spread
    // by stride sampling (a wrong-side region of a few hundred faces then
    // gets several seeds).  Stored as one global CSR so the expensive ray
    // casts are parallel over *seeds*, not over patches.
    // Seeds are drawn only from the *larger* faces of a patch: slivers along
    // intersection curves sit numerically on top of the other sheet and give
    // unreliable crossing counts.
    std::vector<uint32_t> sstart(np + 1, 0), seedFace;
    {
        RawArray<double> area(nt);
        cork_par::for_each_idx(nt, 8192, [&](size_t i) {
            const Vec3d &a = verts[tris[i].a].pos, &b = verts[tris[i].b].pos, &c = verts[tris[i].c].pos;
            area[i] = len(cross(b - a, c - a));
        });
        // sort each patch's face list by area (descending), in place
        cork_par::for_each_idx(np, 16, [&](size_t p) {
            std::sort(plist.begin() + pstart[p], plist.begin() + pstart[p + 1],
                      [&](uint32_t x, uint32_t y) { return area[x] > area[y]; });
        });
        std::vector<uint32_t> ncand(np);
        for (size_t p = 0; p < np; ++p) {
            const uint32_t *faces = &plist[pstart[p]];
            size_t size = pstart[p + 1] - pstart[p];
            double amax = area[faces[0]];
            // candidates: faces at least 10% of the largest face's area,
            // but at least the top quarter of the patch
            size_t nc = 0;
            while (nc < size && area[faces[nc]] >= 0.1 * amax) ++nc;
            nc = std::max<size_t>(nc, std::max<size_t>(1, size / 4));
            ncand[p] = (uint32_t)nc;
            // Pile: ray the real sheets first (size >= 16).  Slivers follow
            // by topology, then a 1-seed fallback if still UNSET.
            size_t S;
            if (flood) S = 1;
            else if (pile && size < 16) S = 0;
            else if (pile) S = 1;
            else S = std::max<size_t>(5, nc / 50);
            if (S > 65536) S = 65536;
            if (S > nc) S = nc;
            sstart[p + 1] = sstart[p] + (uint32_t)S;
        }
        seedFace.resize(sstart[np]);
        cork_par::for_each_idx(np, 64, [&](size_t p) {
            const uint32_t *faces = &plist[pstart[p]];
            size_t nc = ncand[p];
            size_t S = sstart[p + 1] - sstart[p];
            if (S == 0) return;
            uint32_t *out = &seedFace[sstart[p]];
            if (S >= nc) { for (size_t i = 0; i < nc; ++i) out[i] = faces[i]; }
            else {
                size_t stride = nc / S;
                size_t off = (size_t)(mix64(p) % stride);
                for (size_t i = 0; i < S; ++i) out[i] = faces[off + i * stride];
            }
        });
        nSeeds = seedFace.size();
    }
    std::vector<int8_t> seedCode(seedFace.size(), HULL_UNSET);
    {
        CORK_PROF("  hull: rays");
        cork_par::Local<size_t> tls_rays([] { return (size_t)0; });
        cork_par::Local<std::vector<std::pair<uint32_t,int>>> tls_hits;
        cork_par::for_each_idx(seedFace.size(), 16, [&](size_t i) {
            seedCode[i] = classifyFace(seedFace[i], 0x9E3779B1ULL, tls_hits.local(), tls_rays.local());
        });
        tls_rays.combine_each([&](size_t r) { nRays += r; });
    }
    {
        CORK_PROF("  hull: classify patches");
        struct Cnt { size_t split = 0, unresolved = 0; };
        cork_par::Local<Cnt> tls_cnt;
        cork_par::Local<std::vector<uint32_t>> tls_queue;
        cork_par::for_each_idx(np, 1, [&](size_t p) {
            auto &cnt = tls_cnt.local();
            const uint32_t *faces = &plist[pstart[p]];
            const size_t size = pstart[p + 1] - pstart[p];
            const uint32_t *seeds = &seedFace[sstart[p]];
            const int8_t   *codes = &seedCode[sstart[p]];
            const size_t nS = sstart[p + 1] - sstart[p];

            int hist[3] = {0, 0, 0};
            size_t resolved = 0;
            for (size_t i = 0; i < nS; ++i)
                if (codes[i] != HULL_UNSET) { hist[codes[i]]++; ++resolved; }
            if (resolved == 0) {
                cnt.unresolved++;
                for (size_t i = 0; i < size; ++i) faceCode[faces[i]] = HULL_DELETE;
                return;
            }
            int best = 0;
            for (int k = 1; k < 3; ++k) if (hist[k] > hist[best]) best = k;
            const size_t minority = resolved - (size_t)hist[best];
            // Unanimous, or a small patch / lone dissenter: majority vote.
            // Only a large patch with several dissenting seeds is treated as
            // two regions glued through a gap in an intersection curve.
            if (minority == 0 || size < 200 || minority < 2) {
                for (size_t i = 0; i < size; ++i) faceCode[faces[i]] = (int8_t)best;
                return;
            }
            // seeds disagree: multi-source BFS, nearest seed wins
            cnt.split++;
            if (dbg) {
                static std::mutex mu; std::lock_guard<std::mutex> lk(mu);
                std::cerr << "[hull] patch " << p << " (" << size << " faces) split: keep/flip/del seeds = "
                          << hist[0] << "/" << hist[1] << "/" << hist[2] << std::endl;
            }
            auto &queue = tls_queue.local();
            queue.clear();
            for (size_t i = 0; i < nS; ++i) {
                if (codes[i] == HULL_UNSET) continue;
                if (faceCode[seeds[i]] != HULL_UNSET) continue;   // duplicate seed
                faceCode[seeds[i]] = codes[i];
                queue.push_back(seeds[i]);
            }
            for (size_t qi = 0; qi < queue.size(); ++qi) {
                uint32_t f = queue[qi];
                int8_t code = faceCode[f];
                for (int k = 0; k < 3; ++k) {
                    uint32_t g = nb[3 * f + k];
                    if (g == ~0u || faceCode[g] != HULL_UNSET) continue;
                    faceCode[g] = code;
                    queue.push_back(g);
                }
            }
            for (size_t i = 0; i < size; ++i)
                if (faceCode[faces[i]] == HULL_UNSET) faceCode[faces[i]] = HULL_DELETE;
        });
        tls_cnt.combine_each([&](const Cnt &c) { nSplit += c.split; nUnresolved += c.unresolved; });
        cork_prof::note("  hull: #seeds", (double)nSeeds);
        cork_prof::note("  hull: #rays", (double)nRays);
        cork_prof::note("  hull: #patches split", (double)nSplit);
    }

    // ------------------------------------------------------------------
    // 3a. topological vote for small patches.  Rays from tiny sliver
    //     patches (a few faces between nearly coincident sheets) are the
    //     least reliable, but their fate is dictated by their neighbours:
    //     along every intersection edge exactly two faces must survive.
    //     For each small patch count, per incident non-manifold edge, the
    //     kept faces of *other* patches: one -> this patch is needed,
    //     two or more (or none) -> it must go.
    // ------------------------------------------------------------------
    size_t nTopoFixed = 0;
    {
        CORK_PROF("  hull: topo vote (small patches)");
        const size_t SMALL = flood ? ~(size_t)0 : 100;
        for (int iter = 0; iter < 3; ++iter) {
            size_t changed = 0;
            for (size_t p = 0; p < np; ++p) {
                const uint32_t *faces = &plist[pstart[p]];
                const size_t size = pstart[p + 1] - pstart[p];
                if (size >= SMALL) continue;
                int keepV = 0, delV = 0;
                for (size_t i = 0; i < size; ++i) {
                    uint32_t f = faces[i];
                    for (int k = 0; k < 3; ++k) {
                        uint32_t eid = nmEdgeOf[3 * f + k];
                        if (eid == ~0u) continue;
                        int others = 0, total = 0;
                        for (uint32_t q = nmStart[eid]; q < nmStart[eid + 1]; ++q) {
                            uint32_t g = nmFace[q];
                            if (patchOf[g] == p) continue;
                            ++total;
                            if (faceCode[g] != HULL_DELETE) ++others;
                        }
                        if (total == 0) continue;            // open edge of the input
                        if (others == 1) ++keepV; else ++delV;
                    }
                }
                if (keepV + delV == 0) continue;
                int8_t cur = faceCode[faces[0]];
                int8_t want;
                if (keepV > delV) want = (cur == HULL_FLIP) ? HULL_FLIP : HULL_KEEP;
                else              want = HULL_DELETE;
                if (want != cur) {
                    for (size_t i = 0; i < size; ++i) faceCode[faces[i]] = want;
                    ++changed;
                }
            }
            nTopoFixed += changed;
            if (changed == 0) break;
        }
        cork_prof::note("  hull: small patches re-decided by topology", (double)nTopoFixed);
    }

    // ------------------------------------------------------------------
    // 3b. prune leftover scraps on the fly (no second mesh walk).
    //     - slivers: small kept comps (<100 faces) with a non-manifold /
    //       open boundary (between coincident sheets)
    //     - Mira leftover shells: any kept component whose area is
    //       < 0.001 * largest component, including closed junk
    // ------------------------------------------------------------------
    size_t nPruned = 0;
    bool hullClosed = true;
    {
        CORK_PROF("  hull: prune slivers");
        const size_t SMALL = 100;
        for (int iter = 0; iter < 4; ++iter) {
            struct EK { uint64_t key; uint32_t tri; };
            std::vector<uint32_t> kept;
            kept.reserve(nt);
            for (size_t i = 0; i < nt; ++i) if (faceCode[i] != HULL_DELETE) kept.push_back((uint32_t)i);
            const size_t nk = kept.size();
            if (nk == 0) break;
            RawArray<EK> ek(3 * nk);
            cork_par::for_each_idx(nk, 8192, [&](size_t j) {
                const Tri &t = tris[kept[j]];
                for (uint k = 0; k < 3; ++k) {
                    uint a = t.v[(k + 1) % 3], b = t.v[(k + 2) % 3];
                    uint lo = std::min(a, b), hi = std::max(a, b);
                    ek[3 * j + k].key = ((uint64_t)lo << 32) | (uint64_t)hi;
                    ek[3 * j + k].tri = (uint32_t)j;
                }
            });
            cork_par::sort(ek.begin(), ek.end(), [](const EK &x, const EK &y) {
                return x.key < y.key || (x.key == y.key && x.tri < y.tri);
            });
            std::vector<uint32_t> parent(nk);
            std::vector<uint8_t> hasBoundary(nk, 0);
            for (size_t j = 0; j < nk; ++j) parent[j] = (uint32_t)j;
            auto find = [&](uint32_t x) {
                while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
                return x;
            };
            const size_t n3 = ek.size();
            for (size_t i = 0; i < n3;) {
                size_t j = i + 1;
                while (j < n3 && ek[j].key == ek[i].key) ++j;
                if (j - i == 2) {
                    uint32_t a = find(ek[i].tri), b = find(ek[i + 1].tri);
                    if (a != b) parent[std::max(a, b)] = std::min(a, b);
                } else {
                    for (size_t k = i; k < j; ++k) hasBoundary[ek[k].tri] = 1;
                }
                i = j;
            }
            std::vector<uint32_t> compSize(nk, 0);
            std::vector<uint8_t>  compBoundary(nk, 0);
            std::vector<double>   compArea(nk, 0);
            std::vector<uint32_t> root(nk);
            for (size_t j = 0; j < nk; ++j) {
                root[j] = find((uint32_t)j);
                compSize[root[j]]++;
                compBoundary[root[j]] |= hasBoundary[j];
                const Tri &t = tris[kept[j]];
                compArea[root[j]] += 0.5 * len(cross(verts[t.a].pos - verts[t.b].pos,
                                                     verts[t.c].pos - verts[t.b].pos));
            }
            double maxA = 0;
            for (size_t j = 0; j < nk; ++j)
                if (compArea[j] > maxA) maxA = compArea[j];
            const double areaCut = leftoverAreaFrac > 0 ? leftoverAreaFrac * maxA : -1;
            size_t removed = 0;
            for (size_t j = 0; j < nk; ++j) {
                uint32_t r = root[j];
                const bool sliver = compBoundary[r] && compSize[r] < SMALL;
                const bool leftover = leftoverAreaFrac > 0 && compArea[r] < areaCut;
                if (sliver || leftover) { faceCode[kept[j]] = HULL_DELETE; ++removed; }
            }
            nPruned += removed;
            hullClosed = true;
            for (size_t j = 0; j < nk; ++j) {
                if (faceCode[kept[j]] == HULL_DELETE) continue;
                if (compBoundary[root[j]]) { hullClosed = false; break; }
            }
            if (removed == 0) break;
        }
        cork_prof::note("  hull: sliver faces pruned", (double)nPruned);
    }

    // Cork isClosed() is directed: each edge (a→b)+1, (b→a)−1.
    // Undirected valence-2 can still fail isSolid if neighbors agree.
    // One BFS on the last kept graph — same cost class as one prune pass.
    if (hullClosed) {
        CORK_PROF("  hull: orient (cork isClosed)");
        std::vector<uint32_t> kept;
        kept.reserve(nt);
        for (size_t i = 0; i < nt; ++i)
            if (faceCode[i] != HULL_DELETE) kept.push_back((uint32_t)i);
        const size_t nk = kept.size();
        if (nk >= 2) {
            struct EK { uint64_t key; uint32_t tri; };
            RawArray<EK> ek(3 * nk);
            cork_par::for_each_idx(nk, 8192, [&](size_t j) {
                const Tri &t = tris[kept[j]];
                for (uint k = 0; k < 3; ++k) {
                    uint a = t.v[(k + 1) % 3], b = t.v[(k + 2) % 3];
                    uint lo = std::min(a, b), hi = std::max(a, b);
                    ek[3 * j + k].key = ((uint64_t)lo << 32) | (uint64_t)hi;
                    ek[3 * j + k].tri = (uint32_t)j;
                }
            });
            cork_par::sort(ek.begin(), ek.end(), [](const EK &x, const EK &y) {
                return x.key < y.key || (x.key == y.key && x.tri < y.tri);
            });
            struct Nbr { uint32_t tri; uint32_t lo, hi; };
            std::vector<Nbr> adj(nk * 3);
            std::vector<uint8_t> nadj(nk, 0);
            const size_t n3 = ek.size();
            for (size_t i = 0; i < n3;) {
                size_t j = i + 1;
                while (j < n3 && ek[j].key == ek[i].key) ++j;
                if (j - i == 2) {
                    const uint32_t a = ek[i].tri, b = ek[i + 1].tri;
                    const uint32_t lo = (uint32_t)(ek[i].key >> 32);
                    const uint32_t hi = (uint32_t)ek[i].key;
                    if (nadj[a] < 3) adj[a * 3 + nadj[a]++] = {b, lo, hi};
                    if (nadj[b] < 3) adj[b * 3 + nadj[b]++] = {a, lo, hi};
                }
                i = j;
            }
            auto uses_dir = [&](uint32_t tid, bool flip, uint32_t u, uint32_t v) {
                uint a = tris[tid].a, b = tris[tid].b, c = tris[tid].c;
                if (flip) std::swap(a, b);
                return (a == u && b == v) || (b == u && c == v) || (c == u && a == v);
            };
            std::vector<char> vis(nk, 0);
            std::vector<uint32_t> stack;
            stack.reserve(nk);
            for (size_t s = 0; s < nk; ++s) {
                if (vis[s]) continue;
                vis[s] = 1;
                stack.push_back((uint32_t)s);
                while (!stack.empty()) {
                    const uint32_t j = stack.back();
                    stack.pop_back();
                    const uint32_t tid = kept[j];
                    const bool flipJ = faceCode[tid] == HULL_FLIP;
                    for (uint e = 0; e < nadj[j]; ++e) {
                        const Nbr &nb = adj[j * 3 + e];
                        if (vis[nb.tri]) continue;
                        const uint32_t uid = kept[nb.tri];
                        const bool flipK = faceCode[uid] == HULL_FLIP;
                        if (uses_dir(tid, flipJ, nb.lo, nb.hi) ==
                            uses_dir(uid, flipK, nb.lo, nb.hi))
                            faceCode[uid] = flipK ? HULL_KEEP : HULL_FLIP;
                        vis[nb.tri] = 1;
                        stack.push_back(nb.tri);
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // 4. rebuild verts / tris
    // ------------------------------------------------------------------
    {
        CORK_PROF("  hull: rebuild");
        const size_t nv = verts.size();
        std::vector<uint8_t> keepT(nt);
        cork_par::for_each_idx(nt, 8192, [&](size_t i) { keepT[i] = faceCode[i] != HULL_DELETE; });
        std::vector<uint32_t> tmap(nt + 1);
        tmap[0] = 0;
        for (size_t i = 0; i < nt; ++i) tmap[i + 1] = tmap[i] + keepT[i];
        const size_t nt2 = tmap[nt];

        std::vector<uint8_t> usedV(nv, 0);
        cork_par::for_each_idx(nt, 8192, [&](size_t i) {
            if (!keepT[i]) return;
            usedV[tris[i].a] = 1; usedV[tris[i].b] = 1; usedV[tris[i].c] = 1;
        });
        std::vector<uint32_t> vmap(nv + 1);
        vmap[0] = 0;
        for (size_t i = 0; i < nv; ++i) vmap[i + 1] = vmap[i] + usedV[i];
        const size_t nv2 = vmap[nv];

        std::vector<VertData> nverts(nv2);
        cork_par::for_each_idx(nv, 8192, [&](size_t i) { if (usedV[i]) nverts[vmap[i]] = verts[i]; });
        std::vector<Tri> ntris(nt2);
        size_t kept = 0, flipped = 0;
        cork_par::Local<std::pair<size_t,size_t>> tls_cnt([] { return std::pair<size_t,size_t>(0, 0); });
        cork_par::for_each_idx(nt, 8192, [&](size_t i) {
            if (!keepT[i]) return;
            Tri &dst = ntris[tmap[i]];
            dst = tris[i];
            dst.a = vmap[tris[i].a]; dst.b = vmap[tris[i].b]; dst.c = vmap[tris[i].c];
            auto &cnt = tls_cnt.local();
            if (faceCode[i] == HULL_FLIP) { std::swap(dst.a, dst.b); cnt.second++; }
            else cnt.first++;
        });
        tls_cnt.combine_each([&](const std::pair<size_t,size_t> &c) { kept += c.first; flipped += c.second; });
        verts.swap(nverts);
        tris.swap(ntris);
        if (stats) {
            stats->patches = np; stats->split_patches = nSplit; stats->seeds = nSeeds;
            stats->kept = kept; stats->flipped = flipped; stats->deleted = nt - nt2;
            stats->rays = nRays; stats->unresolved = nUnresolved; stats->pruned = nPruned;
            stats->closed = hullClosed;
        }
        cork_prof::note("  hull: faces kept", (double)kept);
        cork_prof::note("  hull: faces flipped", (double)flipped);
        cork_prof::note("  hull: faces deleted", (double)(nt - nt2));
    }
}
