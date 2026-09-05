// Linear BVH (Karras 2012 / Manifold Collider): Morton-sort leaves, build a
// binary radix tree in parallel, query by DFS.  Used as a drop-in broadphase
// alternative to the uniform-grid merge-join.
#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <vector>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#include <cork/math/bbox.h>
#include <cork/math/vec.h>
#include <cork/util/parallel.h>

namespace cork_lbvh {

struct BoxF {
    float mnx, mny, mnz, mxx, mxy, mxz;
};

inline float fdown(double v) {
    float f = (float)v;
    return ((double)f > v) ? std::nextafter(f, -INFINITY) : f;
}
inline float fup(double v) {
    float f = (float)v;
    return ((double)f < v) ? std::nextafter(f, INFINITY) : f;
}
inline BoxF toBoxF(const BBox3d &b) {
    BoxF f;
    f.mnx = fdown(b.minp.x); f.mny = fdown(b.minp.y); f.mnz = fdown(b.minp.z);
    f.mxx = fup(b.maxp.x);   f.mxy = fup(b.maxp.y);   f.mxz = fup(b.maxp.z);
    return f;
}
inline bool overlapF(const BoxF &a, const BoxF &b) {
    return a.mnx <= b.mxx && b.mnx <= a.mxx &&
           a.mny <= b.mxy && b.mny <= a.mxy &&
           a.mnz <= b.mxz && b.mnz <= a.mxz;
}

// Slab test. inv = 1/dir with a huge value on a zero component.
inline bool rayHitBoxF(const Vec3d &p, const Vec3d &inv, const BoxF &b) {
    double tmin = -1e300, tmax = 1e300;
    const double mn[3] = { b.mnx, b.mny, b.mnz };
    const double mx[3] = { b.mxx, b.mxy, b.mxz };
    const double pv[3] = { p.x, p.y, p.z };
    const double iv[3] = { inv.x, inv.y, inv.z };
    for (int k = 0; k < 3; ++k) {
        if (std::fabs(iv[k]) > 1e200) {
            if (pv[k] < mn[k] || pv[k] > mx[k]) return false;
            continue;
        }
        double t0 = (mn[k] - pv[k]) * iv[k];
        double t1 = (mx[k] - pv[k]) * iv[k];
        if (t0 > t1) std::swap(t0, t1);
        if (t0 > tmin) tmin = t0;
        if (t1 < tmax) tmax = t1;
        if (tmax < tmin) return false;
    }
    return tmax >= 0.0;
}

inline uint32_t spread3(uint32_t v) {
    v = 0xFF0000FFu & (v * 0x00010001u);
    v = 0x0F00F00Fu & (v * 0x00000101u);
    v = 0xC30C30C3u & (v * 0x00000011u);
    v = 0x49249249u & (v * 0x00000005u);
    return v;
}

inline uint32_t clz32(uint32_t x) {
#if defined(_MSC_VER)
    unsigned long i;
    if (!_BitScanReverse(&i, x)) return 32;
    return 31 - (uint32_t)i;
#else
    return x ? (uint32_t)__builtin_clz(x) : 32;
#endif
}

inline uint32_t mortonOf(const Vec3d &p, const Vec3d &org, const Vec3d &ext) {
    auto q = [](double x, double o, double e) -> uint32_t {
        double t = (e > 0.0) ? (x - o) / e : 0.0;
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        uint32_t u = (uint32_t)(t * 1024.0);
        return u > 1023u ? 1023u : u;
    };
    uint32_t x = spread3(q(p.x, org.x, ext.x));
    uint32_t y = spread3(q(p.y, org.y, ext.y));
    uint32_t z = spread3(q(p.z, org.z, ext.z));
    return x * 4 + y * 2 + z;
}

inline bool isLeaf(int node) { return (node & 1) == 0; }
inline int leaf2node(int leaf) { return leaf * 2; }
inline int node2leaf(int node) { return node / 2; }
inline int internal2node(int in) { return in * 2 + 1; }
inline int node2internal(int node) { return (node - 1) / 2; }

struct LBVH {
    std::vector<BoxF> nodeBox;          // 2n-1
    std::vector<int> child1, child2;    // n-1 internals
    std::vector<int> parent;            // 2n-1
    std::vector<uint32_t> leafPrim;     // n, sorted
    int n = 0;

    void build(const BBox3d *boxes, size_t nboxes, const BBox3d &world) {
        n = (int)nboxes;
        if (n <= 0) return;
        leafPrim.resize((size_t)n);
        Vec3d ext = world.maxp - world.minp;
        struct Item { uint32_t code; uint32_t idx; };
        std::vector<Item> items((size_t)n);
        cork_par::for_each_idx((size_t)n, 4096, [&](size_t i) {
            Vec3d c = (boxes[i].minp + boxes[i].maxp) * 0.5;
            items[i].code = mortonOf(c, world.minp, ext);
            items[i].idx = (uint32_t)i;
        });
        cork_par::sort(items.begin(), items.end(), [](const Item &a, const Item &b) {
            return a.code < b.code || (a.code == b.code && a.idx < b.idx);
        });
        cork_par::for_each_idx((size_t)n, 8192, [&](size_t i) { leafPrim[i] = items[i].idx; });

        const int nint = std::max(0, n - 1);
        const int nnode = 2 * n - 1;
        nodeBox.resize((size_t)nnode);
        child1.assign((size_t)nint, 0);
        child2.assign((size_t)nint, 0);
        parent.assign((size_t)nnode, 0);

        if (n == 1) {
            nodeBox[0] = toBoxF(boxes[leafPrim[0]]);
            return;
        }

        std::vector<uint32_t> morton((size_t)n);
        cork_par::for_each_idx((size_t)n, 8192, [&](size_t i) { morton[i] = items[i].code; });

        auto prefix = [&](int i, int j) -> int {
            if (j < 0 || j >= n) return -1;
            if (morton[i] == morton[j])
                return 32 + (int)clz32((uint32_t)i ^ (uint32_t)j);
            return (int)clz32(morton[i] ^ morton[j]);
        };
        auto rangeEnd = [&](int i) -> int {
            int dir = prefix(i, i + 1) - prefix(i, i - 1);
            dir = (dir > 0) - (dir < 0);
            int common = prefix(i, i - dir);
            int maxL = 128;
            while (prefix(i, i + dir * maxL) > common) maxL *= 4;
            int length = 0;
            for (int step = maxL / 2; step > 0; step /= 2)
                if (prefix(i, i + dir * (length + step)) > common) length += step;
            return i + dir * length;
        };
        auto findSplit = [&](int first, int last) -> int {
            int common = prefix(first, last);
            int split = first;
            int step = last - first;
            do {
                step = (step + 1) >> 1;
                int ns = split + step;
                if (ns < last && prefix(first, ns) > common) split = ns;
            } while (step > 1);
            return split;
        };

        cork_par::for_each_idx((size_t)nint, 1024, [&](size_t internal) {
            int first = (int)internal;
            int last = rangeEnd(first);
            if (first > last) std::swap(first, last);
            int split = findSplit(first, last);
            int c1 = (split == first) ? leaf2node(split) : internal2node(split);
            ++split;
            int c2 = (split == last) ? leaf2node(split) : internal2node(split);
            child1[internal] = c1;
            child2[internal] = c2;
            int node = internal2node((int)internal);
            parent[c1] = node;
            parent[c2] = node;
        });

        cork_par::for_each_idx((size_t)n, 4096, [&](size_t leaf) {
            nodeBox[leaf2node((int)leaf)] = toBoxF(boxes[leafPrim[leaf]]);
        });
        std::vector<std::atomic<int>> ctr((size_t)nint);
        cork_par::for_each_idx((size_t)nint, 4096, [&](size_t i) { ctr[i].store(0, std::memory_order_relaxed); });
        cork_par::for_each_idx((size_t)n, 1024, [&](size_t leaf) {
            int node = leaf2node((int)leaf);
            for (;;) {
                node = parent[node];
                int in = node2internal(node);
                if (in < 0 || in >= nint) break;
                int old = ctr[in].fetch_add(1, std::memory_order_acq_rel);
                if (old == 0) return;   // first child here; the other will union
                const BoxF &a = nodeBox[child1[in]];
                const BoxF &b = nodeBox[child2[in]];
                BoxF u;
                u.mnx = std::min(a.mnx, b.mnx); u.mny = std::min(a.mny, b.mny); u.mnz = std::min(a.mnz, b.mnz);
                u.mxx = std::max(a.mxx, b.mxx); u.mxy = std::max(a.mxy, b.mxy); u.mxz = std::max(a.mxz, b.mxz);
                nodeBox[node] = u;
                if (node == 1) break;
            }
        });
    }

    // Infinite ray from p with inv = 1/direction (huge if a component is 0).
    // Calls fn(leafPrim[leaf]) for every leaf box the ray overlaps.
    template<class Fn>
    void ray(const Vec3d &p, const Vec3d &inv, Fn &&fn) const {
        if (n <= 0) return;
        if (n == 1) {
            if (rayHitBoxF(p, inv, nodeBox[0])) fn(leafPrim[0]);
            return;
        }
        int fixed[1024];
        std::vector<int> extra;
        int top = -1;
        auto push = [&](int n) {
            if (top + 1 < 1024) fixed[++top] = n;
            else extra.push_back(n);
        };
        auto pop = [&]() -> int {
            if (!extra.empty()) { int n = extra.back(); extra.pop_back(); return n; }
            return fixed[top--];
        };
        int node = 1;
        while (1) {
            int in = node2internal(node);
            int c1 = child1[in], c2 = child2[in];
            bool t1 = rayHitBoxF(p, inv, nodeBox[c1]);
            bool t2 = rayHitBoxF(p, inv, nodeBox[c2]);
            int go1 = t1 && !isLeaf(c1);
            int go2 = t2 && !isLeaf(c2);
            if (t1 && isLeaf(c1)) fn(leafPrim[node2leaf(c1)]);
            if (t2 && isLeaf(c2)) fn(leafPrim[node2leaf(c2)]);
            if (!go1 && !go2) {
                if (extra.empty() && top < 0) break;
                node = pop();
            } else {
                node = go1 ? c1 : c2;
                if (go1 && go2) push(c2);
            }
        }
    }

    // Call fn(leafPrim[leaf]) for every leaf whose box overlaps q.
    template<class Fn>
    void query(const BoxF &q, Fn &&fn) const {
        if (n <= 0) return;
        if (n == 1) {
            if (overlapF(q, nodeBox[0])) fn(leafPrim[0]);
            return;
        }
        int stack[64];
        int top = -1;
        int node = 1;
        while (1) {
            int in = node2internal(node);
            int c1 = child1[in], c2 = child2[in];
            bool t1 = overlapF(q, nodeBox[c1]);
            bool t2 = overlapF(q, nodeBox[c2]);
            int go1 = t1 && !isLeaf(c1);
            int go2 = t2 && !isLeaf(c2);
            if (t1 && isLeaf(c1)) fn(leafPrim[node2leaf(c1)]);
            if (t2 && isLeaf(c2)) fn(leafPrim[node2leaf(c2)]);
            if (!go1 && !go2) {
                if (top < 0) break;
                node = stack[top--];
            } else {
                node = go1 ? c1 : c2;
                if (go1 && go2) stack[++top] = c2;
            }
        }
    }
};

} // namespace cork_lbvh
