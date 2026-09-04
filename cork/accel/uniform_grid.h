// Uniform grid + stamp-set uniqueness (no hash clears per query)
#ifndef CORK_UNIFORM_GRID_H
#define CORK_UNIFORM_GRID_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <cork/math/bbox.h>
#include <cork/math/vec.h>

struct UniformGrid {
    BBox3d world;
    int nx = 1, ny = 1, nz = 1;
    Vec3d cell;
    Vec3d inv_cell;
    std::vector<std::vector<int>> bins;

    void init(const BBox3d &w, int nprims)
    {
        world = w;
        Vec3d d = world.maxp - world.minp;
        double extent = std::max({std::abs(d.x), std::abs(d.y), std::abs(d.z), 1.0});
        double pad = 1e-6 * extent;
        world.minp = world.minp - Vec3d(pad, pad, pad);
        world.maxp = world.maxp + Vec3d(pad, pad, pad);
        d = world.maxp - world.minp;

        int target = (int)std::cbrt(std::max(1.0, nprims / 6.0));
        target = std::max(24, std::min(target, 192));
        nx = ny = nz = target;
        cell = Vec3d(d.x / nx, d.y / ny, d.z / nz);
        inv_cell = Vec3d(1.0 / cell.x, 1.0 / cell.y, 1.0 / cell.z);
        bins.assign((size_t)nx * (size_t)ny * (size_t)nz, {});
    }

    inline int flat(int x, int y, int z) const {
        return (z * ny + y) * nx + x;
    }

    inline void clamp_cell(int &x, int &y, int &z) const {
        x = std::max(0, std::min(nx - 1, x));
        y = std::max(0, std::min(ny - 1, y));
        z = std::max(0, std::min(nz - 1, z));
    }

    inline void cell_range(const BBox3d &b, int &x0, int &y0, int &z0,
                           int &x1, int &y1, int &z1) const
    {
        x0 = (int)std::floor((b.minp.x - world.minp.x) * inv_cell.x);
        y0 = (int)std::floor((b.minp.y - world.minp.y) * inv_cell.y);
        z0 = (int)std::floor((b.minp.z - world.minp.z) * inv_cell.z);
        x1 = (int)std::floor((b.maxp.x - world.minp.x) * inv_cell.x);
        y1 = (int)std::floor((b.maxp.y - world.minp.y) * inv_cell.y);
        z1 = (int)std::floor((b.maxp.z - world.minp.z) * inv_cell.z);
        clamp_cell(x0, y0, z0);
        clamp_cell(x1, y1, z1);
    }

    void insert(int prim, const BBox3d &b)
    {
        int x0, y0, z0, x1, y1, z1;
        cell_range(b, x0, y0, z0, x1, y1, z1);
        for (int z = z0; z <= z1; ++z)
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x)
                    bins[(size_t)flat(x, y, z)].push_back(prim);
    }
};

struct BBoxSoA {
    std::vector<float> minx, miny, minz, maxx, maxy, maxz;
    void resize(size_t n) {
        minx.resize(n); miny.resize(n); minz.resize(n);
        maxx.resize(n); maxy.resize(n); maxz.resize(n);
    }
    void set(size_t i, const BBox3d &b) {
        minx[i] = (float)b.minp.x; miny[i] = (float)b.minp.y; minz[i] = (float)b.minp.z;
        maxx[i] = (float)b.maxp.x; maxy[i] = (float)b.maxp.y; maxz[i] = (float)b.maxp.z;
    }
};

// Thread-local stamp uniqueness (epoch++) — O(1) clear
struct StampSet {
    std::vector<uint32_t> stamp;
    uint32_t epoch = 1;
    void ensure(size_t n) {
        if (stamp.size() < n) stamp.assign(n, 0);
    }
    void clear() {
        if (++epoch == 0) {
            std::fill(stamp.begin(), stamp.end(), 0);
            epoch = 1;
        }
    }
    bool insert(int id) {
        if (stamp[(size_t)id] == epoch) return false;
        stamp[(size_t)id] = epoch;
        return true;
    }
};

#endif
