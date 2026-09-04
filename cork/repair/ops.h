#pragma once

// Mira native kernels (external/mira_module, mira_repair wrapper), not the
// Python glue. All work on CorkMesh via raw() so mesh.isct.tpp stays put.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <vector>

#include <cork/accel/uniform_grid.h>
#include <cork/cork.h>
#include <cork/math/bbox.h>
#include <cork/repair/options.h>
#include <cork/util/profile.h>
#include <cork/util/unionFind.h>

namespace cork {
namespace repair {

using Raw = RawMesh<CorkVertex, CorkTriangle>;

inline void adopt(CorkMesh &mesh, Raw &&raw) {
    CorkMesh next(std::move(raw));
    mesh = std::move(next);
}

inline Raw compact(const Raw &in, const std::vector<char> &keep_face) {
    const size_t nV = in.vertices.size();
    const size_t nT = in.triangles.size();
    std::vector<int> remap(nV, -1);
    Raw out;
    out.vertices.reserve(nV);
    out.triangles.reserve(nT);
    for (size_t t = 0; t < nT; ++t) {
        if (t < keep_face.size() && !keep_face[t]) continue;
        CorkTriangle nt = in.triangles[t];
        const int ids[3] = {nt.a, nt.b, nt.c};
        int mapped[3];
        for (int k = 0; k < 3; ++k) {
            const int v = ids[k];
            if (v < 0 || (size_t)v >= nV) {
                mapped[k] = 0;
                continue;
            }
            if (remap[v] < 0) {
                remap[v] = (int)out.vertices.size();
                out.vertices.push_back(in.vertices[v]);
            }
            mapped[k] = remap[v];
        }
        if (mapped[0] == mapped[1] || mapped[1] == mapped[2] || mapped[2] == mapped[0])
            continue;
        nt.a = mapped[0];
        nt.b = mapped[1];
        nt.c = mapped[2];
        out.triangles.push_back(nt);
    }
    return out;
}

inline void drop_degen_unref(CorkMesh &mesh) {
    CORK_PROF("repair.clean");
    Raw raw = mesh.raw();
    const size_t nT = raw.triangles.size();
    std::vector<char> keep(nT, 0);
    size_t dropped = 0;
    for (size_t t = 0; t < nT; ++t) {
        const auto &tr = raw.triangles[t];
        if (tr.a == tr.b || tr.b == tr.c || tr.c == tr.a) { ++dropped; continue; }
        const Vec3d &a = raw.vertices[tr.a].pos;
        const Vec3d &b = raw.vertices[tr.b].pos;
        const Vec3d &c = raw.vertices[tr.c].pos;
        const Vec3d e1 = b - a;
        const Vec3d e2 = c - a;
        const Vec3d e3 = a - c;
        const double area = 0.5 * len(cross(e1, e2));
        const double emax = std::max({len(e1), len(b - c), len(e3)});
        if (emax <= 0.0 || area / emax <= 1e-7) { ++dropped; continue; }
        keep[t] = 1;
    }
    if (!dropped) return;
    adopt(mesh, compact(raw, keep));
}

inline void perturb_along_normals(CorkMesh &mesh, double intensity) {
    CORK_PROF("repair.perturb");
    if (intensity == 0.0) return;
    Raw raw = mesh.raw();
    const size_t nV = raw.vertices.size();
    const size_t nT = raw.triangles.size();
    std::vector<Vec3d> nrm(nV, Vec3d(0, 0, 0));
    for (size_t t = 0; t < nT; ++t) {
        const auto &tr = raw.triangles[t];
        const Vec3d n = cross(raw.vertices[tr.b].pos - raw.vertices[tr.a].pos,
                              raw.vertices[tr.c].pos - raw.vertices[tr.a].pos);
        nrm[tr.a] += n;
        nrm[tr.b] += n;
        nrm[tr.c] += n;
    }
    const uint64_t seed = (uint64_t)std::rand() ^ 0xC0FFEEULL;
    for (size_t i = 0; i < nV; ++i) {
        double L = len(nrm[i]);
        if (L <= 0.0) continue;
        uint64_t x = seed + 0x9e3779b97f4a7c15ULL * (i + 1);
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        x ^= x >> 31;
        const double u = (double)(x >> 11) * (1.0 / 9007199254740992.0);
        raw.vertices[i].pos += (intensity * u) * (nrm[i] / L);
    }
    adopt(mesh, std::move(raw));
}

inline uint64_t edge_key(int a, int b) {
    const uint32_t u = (uint32_t)std::min(a, b);
    const uint32_t v = (uint32_t)std::max(a, b);
    return ((uint64_t)u << 32) | v;
}

struct FaceKey {
    uint32_t a, b, c;
    bool operator==(const FaceKey &o) const { return a == o.a && b == o.b && c == o.c; }
};
struct FaceKeyHash {
    size_t operator()(const FaceKey &k) const {
        size_t h = k.a;
        h ^= (size_t)k.b + 0x9e3779b9u + (h << 6) + (h >> 2);
        h ^= (size_t)k.c + 0x9e3779b9u + (h << 6) + (h >> 2);
        return h;
    }
};
// Mira post_repair_fix / merge_vertices: identical positions only.
// Bit-identity — no distance cutoff (those punched features).
inline int weld_coincident(CorkMesh &mesh) {
    CORK_PROF("repair.weld");
    Raw raw = mesh.raw();
    const size_t nV = raw.vertices.size();
    const size_t nT = raw.triangles.size();
    if (nV == 0 || nT == 0) return 0;

    struct PosKey {
        uint64_t x, y, z;
        bool operator==(const PosKey &o) const {
            return x == o.x && y == o.y && z == o.z;
        }
    };
    struct PosKeyHash {
        size_t operator()(const PosKey &k) const {
            size_t h = (size_t)k.x;
            h ^= (size_t)k.y + 0x9e3779b9u + (h << 6) + (h >> 2);
            h ^= (size_t)k.z + 0x9e3779b9u + (h << 6) + (h >> 2);
            return h;
        }
    };
    auto key_of = [](const Vec3d &p) -> PosKey {
        PosKey k;
        std::memcpy(&k.x, &p.x, sizeof(double));
        std::memcpy(&k.y, &p.y, sizeof(double));
        std::memcpy(&k.z, &p.z, sizeof(double));
        return k;
    };

    UnionFind uf((uint)nV);
    std::unordered_map<PosKey, int, PosKeyHash> first;
    first.reserve(nV);
    int merges = 0;
    for (size_t i = 0; i < nV; ++i) {
        const PosKey k = key_of(raw.vertices[i].pos);
        auto it = first.find(k);
        if (it == first.end())
            first.emplace(k, (int)i);
        else {
            uf.unionIds((uint)it->second, (uint)i);
            ++merges;
        }
    }
    if (!merges) return 0;

    for (auto &tr : raw.triangles) {
        tr.a = (int)uf.find((uint)tr.a);
        tr.b = (int)uf.find((uint)tr.b);
        tr.c = (int)uf.find((uint)tr.c);
    }
    std::vector<char> keep(nT, 1);
    adopt(mesh, compact(raw, keep));
    return merges;
}

inline FaceKey face_key(int a, int b, int c) {
    uint32_t i = (uint32_t)a, j = (uint32_t)b, k = (uint32_t)c;
    if (i > j) std::swap(i, j);
    if (j > k) std::swap(j, k);
    if (i > j) std::swap(i, j);
    return {i, j, k};
}

// Mira remove_noise BEFORE cork ("if the noise is kept, cork crashes"):
// 1) drop exact duplicate faces
// 2) drop face-adjacency comps with nFaces <= minFaces
// Planar-vs-first-face-of-mesh is NOT used here — that punched slc 21.stl
// (a large planar shell || n0). Post-hull leftover area covers planar scraps.
inline int drop_noise(CorkMesh &mesh, int minFaces, double areaFrac = 0.0,
                      bool dropTiny = true) {
    CORK_PROF("repair.noise");
    Raw raw = mesh.raw();
    const size_t nT0 = raw.triangles.size();
    if (nT0 == 0) return 0;

    std::unordered_map<FaceKey, char, FaceKeyHash> seenF;
    seenF.reserve(nT0);
    std::vector<char> uniq(nT0, 0);
    size_t nT = 0;
    for (size_t t = 0; t < nT0; ++t) {
        const auto &tr = raw.triangles[t];
        if (tr.a == tr.b || tr.b == tr.c || tr.c == tr.a) continue;
        if (seenF.emplace(face_key(tr.a, tr.b, tr.c), 1).second) {
            uniq[t] = 1;
            ++nT;
        }
    }
    if (nT == 0) {
        adopt(mesh, Raw());
        return (int)nT0;
    }
    // Hull leftover already drops tiny scraps. Skip the component walk
    // when we only needed exact dups (Mira's cork-crash dups, no UF).
    if (!dropTiny) {
        if (nT == nT0) return 0;
        std::vector<char> keep(nT0, 0);
        for (size_t t = 0; t < nT0; ++t) keep[t] = uniq[t];
        adopt(mesh, compact(raw, keep));
        return (int)(nT0 - nT);
    }

    std::vector<int> id(nT0, -1);
    int next = 0;
    for (size_t t = 0; t < nT0; ++t)
        if (uniq[t]) id[t] = next++;

    UnionFind uf((uint)nT);
    std::unordered_map<uint64_t, int> first;
    first.reserve(nT * 2);
    auto add = [&](int a, int b, int t) {
        const uint64_t k = edge_key(a, b);
        auto it = first.find(k);
        if (it == first.end())
            first.emplace(k, t);
        else
            uf.unionIds((uint)it->second, (uint)t);
    };
    std::vector<double> area(nT, 0);
    for (size_t t = 0; t < nT0; ++t) {
        if (!uniq[t]) continue;
        const int ti = id[t];
        const auto &tr = raw.triangles[t];
        add(tr.a, tr.b, ti);
        add(tr.b, tr.c, ti);
        add(tr.c, tr.a, ti);
        area[ti] = 0.5 * len(cross(raw.vertices[tr.b].pos - raw.vertices[tr.a].pos,
                                   raw.vertices[tr.c].pos - raw.vertices[tr.a].pos));
    }

    std::vector<int> root(nT), count(nT, 0);
    std::vector<double> carea(nT, 0);
    for (int t = 0; t < (int)nT; ++t) {
        const int r = (int)uf.find((uint)t);
        root[t] = r;
        count[r]++;
        carea[r] += area[t];
    }
    double maxA = 0;
    for (int t = 0; t < (int)nT; ++t)
        if (carea[t] > maxA) maxA = carea[t];
    const double areaCut = (areaFrac > 0.0) ? areaFrac * maxA : 0.0;

    std::vector<char> keepRoot(nT, 0);
    for (int t = 0; t < (int)nT; ++t) {
        if (count[t] <= minFaces) continue;
        if (areaCut > 0.0 && carea[t] < areaCut) continue;
        keepRoot[t] = 1;
    }
    std::vector<char> keep(nT0, 0);
    int dropped = 0;
    for (size_t t = 0; t < nT0; ++t) {
        if (!uniq[t]) {
            ++dropped;
            continue;
        }
        if (keepRoot[root[id[t]]])
            keep[t] = 1;
        else
            ++dropped;
    }
    if (!dropped) return 0;
    adopt(mesh, compact(raw, keep));
    return dropped;
}

// Mira post-hull: keep face-adjacency shells with area > frac * max area.
// Does not delete faces from the largest body (avoids punching holes).
inline int drop_small_shells(CorkMesh &mesh, double areaFrac) {
    CORK_PROF("repair.noise.shells");
    if (areaFrac <= 0.0) return 0;
    Raw raw = mesh.raw();
    const size_t nT = raw.triangles.size();
    if (nT == 0) return 0;
    UnionFind uf((uint)nT);
    std::unordered_map<uint64_t, int> first;
    first.reserve(nT * 2);
    std::vector<double> area(nT, 0);
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
        area[t] = 0.5 * len(cross(raw.vertices[tr.b].pos - raw.vertices[tr.a].pos,
                                  raw.vertices[tr.c].pos - raw.vertices[tr.a].pos));
    }
    std::vector<int> root(nT);
    std::vector<double> carea(nT, 0);
    for (size_t t = 0; t < nT; ++t) {
        root[t] = (int)uf.find((uint)t);
        carea[root[t]] += area[t];
    }
    double maxA = 0;
    for (size_t t = 0; t < nT; ++t)
        if (carea[t] > maxA) maxA = carea[t];
    const double cut = areaFrac * maxA;
    std::vector<char> keep(nT, 0);
    int dropped = 0;
    for (size_t t = 0; t < nT; ++t) {
        if (carea[root[t]] > cut)
            keep[t] = 1;
        else
            ++dropped;
    }
    if (!dropped) return 0;
    adopt(mesh, compact(raw, keep));
    return dropped;
}

// mira_module.collapse_edges + repair_thread.graph_collapse
inline void collapse_short_edges(CorkMesh &mesh, double rel) {
    CORK_PROF("repair.collapse");
    Raw raw = mesh.raw();
    const size_t nV = raw.vertices.size();
    const size_t nT = raw.triangles.size();
    if (nV == 0 || nT == 0) return;
    double sum = 0.0;
    size_t ne = 0;
    for (size_t t = 0; t < nT; ++t) {
        const auto &tr = raw.triangles[t];
        sum += len(raw.vertices[tr.b].pos - raw.vertices[tr.a].pos);
        sum += len(raw.vertices[tr.c].pos - raw.vertices[tr.b].pos);
        sum += len(raw.vertices[tr.a].pos - raw.vertices[tr.c].pos);
        ne += 3;
    }
    const double tol = rel * (sum / (double)std::max<size_t>(ne, 1));
    UnionFind uf((uint)nV);
    for (size_t t = 0; t < nT; ++t) {
        const int id[3] = {raw.triangles[t].a, raw.triangles[t].b, raw.triangles[t].c};
        for (int k = 0; k < 3; ++k) {
            const int a = id[k], b = id[(k + 1) % 3];
            if (len(raw.vertices[b].pos - raw.vertices[a].pos) < tol)
                uf.unionIds((uint)a, (uint)b);
        }
    }
    std::vector<Vec3d> acc(nV, Vec3d(0, 0, 0));
    std::vector<int> cnt(nV, 0);
    for (size_t i = 0; i < nV; ++i) {
        const int r = (int)uf.find((uint)i);
        acc[r] += raw.vertices[i].pos;
        cnt[r]++;
    }
    for (size_t i = 0; i < nV; ++i) {
        const int r = (int)uf.find((uint)i);
        raw.vertices[i].pos = acc[r] / (double)cnt[r];
    }
    adopt(mesh, std::move(raw));
    drop_degen_unref(mesh);
}

inline bool share_vertex(const CorkTriangle &a, const CorkTriangle &b) {
    return a.a == b.a || a.a == b.b || a.a == b.c || a.b == b.a || a.b == b.b ||
           a.b == b.c || a.c == b.a || a.c == b.b || a.c == b.c;
}

// Mira get_intersecting_faces: edge-ray vs BVH, skip faces that share the
// edge endpoints. We mark non-adjacent AABB overlaps (same skip rule, our grid).
inline void resolve_si_subset(CorkMesh &mesh) {
    CORK_PROF("repair.si_subset");
    Raw raw = mesh.raw();
    const size_t nT = raw.triangles.size();
    if (nT < 2) {
        mesh.resolveIntersections();
        return;
    }
    std::vector<BBox3d> boxes(nT);
    BBox3d world;
    for (size_t t = 0; t < nT; ++t) {
        const auto &tr = raw.triangles[t];
        const Vec3d &pa = raw.vertices[tr.a].pos;
        const Vec3d &pb = raw.vertices[tr.b].pos;
        const Vec3d &pc = raw.vertices[tr.c].pos;
        BBox3d b(pa, pa);
        b = convex(b, BBox3d(pb, pb));
        b = convex(b, BBox3d(pc, pc));
        boxes[t] = b;
        world = convex(world, b);
    }
    UniformGrid grid;
    grid.init(world, (int)nT);
    for (size_t t = 0; t < nT; ++t) grid.insert((int)t, boxes[t]);

    std::vector<char> mark(nT, 0);
    StampSet seen;
    seen.ensure(nT);
    for (size_t t = 0; t < nT; ++t) {
        int x0, y0, z0, x1, y1, z1;
        grid.cell_range(boxes[t], x0, y0, z0, x1, y1, z1);
        seen.clear();
        for (int z = z0; z <= z1; ++z)
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) {
                    const auto &bin = grid.bins[(size_t)grid.flat(x, y, z)];
                    for (int u : bin) {
                        if ((size_t)u <= t || !seen.insert(u)) continue;
                        if (share_vertex(raw.triangles[t], raw.triangles[u])) continue;
                        if (!hasIsct(boxes[t], boxes[u])) continue;
                        mark[t] = mark[u] = 1;
                    }
                }
    }
    size_t nMark = 0;
    for (char m : mark) nMark += (size_t)m;
    if (nMark == 0) return;
    if (nMark == nT) {
        mesh.resolveIntersections();
        return;
    }
    Raw si = compact(raw, mark);
    std::vector<char> restKeep(nT);
    for (size_t t = 0; t < nT; ++t) restKeep[t] = (char)!mark[t];
    Raw rest = compact(raw, restKeep);
    CorkMesh siMesh(std::move(si));
    siMesh.resolveIntersections();
    CorkMesh restMesh(std::move(rest));
    mesh = std::move(siMesh);
    mesh.disjointUnion(restMesh);
}

}  // namespace repair
}  // namespace cork
