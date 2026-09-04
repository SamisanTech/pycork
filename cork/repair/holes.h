#pragma once

// Boundary-loop fill. Purpose of Mira fillholes_nc / fill_holes_only
// without TMesh: walk open-edge cycles, patch 3-cycles, centroid-star
// the rest. Surface/CorkMesh stay in the caller's dtype (CorkMesh = double).

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include <cork/cork.h>
#include <cork/repair/ops.h>
#include <cork/util/profile.h>

namespace cork {
namespace repair {

inline int fill_boundary_loops(CorkMesh &mesh) {
    CORK_PROF("repair.fill_holes");
    Raw raw = mesh.raw();
    const size_t nT0 = raw.triangles.size();
    if (nT0 == 0) return 0;

    std::unordered_map<uint64_t, int> ec;
    ec.reserve(nT0 * 2);
    for (size_t t = 0; t < nT0; ++t) {
        const auto &tr = raw.triangles[t];
        ec[edge_key(tr.a, tr.b)]++;
        ec[edge_key(tr.b, tr.c)]++;
        ec[edge_key(tr.c, tr.a)]++;
    }

    // Directed hole half-edges: opposite the face edge so the walk
    // goes around the hole (mesh interior stays on the face side).
    std::vector<std::pair<int, int>> he;
    he.reserve(64);
    std::unordered_map<int, std::vector<int>> out;
    auto push_he = [&](int from, int to) {
        const int id = (int)he.size();
        he.push_back({from, to});
        out[from].push_back(id);
    };
    for (size_t t = 0; t < nT0; ++t) {
        const auto &tr = raw.triangles[t];
        const int v[3] = {tr.a, tr.b, tr.c};
        for (int k = 0; k < 3; ++k) {
            const int a = v[k], b = v[(k + 1) % 3];
            if (ec[edge_key(a, b)] == 1) push_he(b, a);
        }
    }
    if (he.empty()) return 0;

    std::vector<char> used(he.size(), 0);
    std::unordered_map<FaceKey, char, FaceKeyHash> existing;
    existing.reserve(nT0);
    for (size_t t = 0; t < nT0; ++t) {
        const auto &tr = raw.triangles[t];
        existing.emplace(face_key(tr.a, tr.b, tr.c), 1);
    }

    int filled = 0;
    for (size_t s = 0; s < he.size(); ++s) {
        if (used[s]) continue;
        std::vector<int> loop;
        int cur = (int)s;
        while (!used[cur]) {
            used[cur] = 1;
            loop.push_back(he[cur].first);
            const int nxt = he[cur].second;
            int follow = -1;
            for (int id : out[nxt]) {
                if (!used[id]) {
                    follow = id;
                    break;
                }
            }
            if (follow < 0) break;
            cur = follow;
        }
        if (loop.size() < 3) continue;
        // unique consecutive verts
        std::vector<int> ring;
        ring.reserve(loop.size());
        for (int v : loop) {
            if (ring.empty() || ring.back() != v) ring.push_back(v);
        }
        if (ring.size() >= 2 && ring.front() == ring.back()) ring.pop_back();
        if (ring.size() < 3) continue;

        auto add_tri = [&](int a, int b, int c) {
            if (a == b || b == c || c == a) return;
            if (!existing.emplace(face_key(a, b, c), 1).second) return;
            CorkTriangle tr{};
            tr.a = a;
            tr.b = b;
            tr.c = c;
            raw.triangles.push_back(tr);
            ++filled;
        };

        if (ring.size() == 3) {
            add_tri(ring[0], ring[1], ring[2]);
            continue;
        }
        // Mira puzzle: fillholes_nc only when hole verts <= 20.
        // Star-capping a leftover hull rim (hundreds of verts) wrecks volume.
        if (ring.size() > 20) continue;

        // centroid-star (TMesh StarTriangulateHole purpose)
        CorkVertex cv;
        cv.pos = Vec3d(0, 0, 0);
        for (int v : ring) cv.pos += raw.vertices[v].pos;
        cv.pos /= (double)ring.size();
        const int cid = (int)raw.vertices.size();
        raw.vertices.push_back(cv);
        for (size_t i = 0; i < ring.size(); ++i)
            add_tri(cid, ring[i], ring[(i + 1) % ring.size()]);
    }

    if (!filled) return 0;
    adopt(mesh, std::move(raw));
    return filled;
}

}  // namespace repair
}  // namespace cork
