// AVX2 / SSE helpers for AABB tests (cork)
#ifndef CORK_BBOX_AVX_H
#define CORK_BBOX_AVX_H

#include <cork/math/bbox.h>

#if defined(_MSC_VER)
#include <intrin.h>
#endif
#include <immintrin.h>

namespace cork_simd {

// True if boxes A and B overlap. Uses SSE2 (always available on x64).
inline bool aabb_overlap(const BBox3d &a, const BBox3d &b)
{
    // a.min <= b.max && b.min <= a.max  (componentwise)
    __m128d amin_xy = _mm_set_pd(a.minp.y, a.minp.x);
    __m128d amax_xy = _mm_set_pd(a.maxp.y, a.maxp.x);
    __m128d bmin_xy = _mm_set_pd(b.minp.y, b.minp.x);
    __m128d bmax_xy = _mm_set_pd(b.maxp.y, b.maxp.x);

    __m128d c0 = _mm_cmple_pd(amin_xy, bmax_xy);
    __m128d c1 = _mm_cmple_pd(bmin_xy, amax_xy);
    __m128d ok_xy = _mm_and_pd(c0, c1);
    // movemask: bit0=x, bit1=y; need both
    if ((_mm_movemask_pd(ok_xy) & 3) != 3)
        return false;
    return a.minp.z <= b.maxp.z && b.minp.z <= a.maxp.z;
}

// Test one query box against N axis-aligned boxes packed as SoA floats.
// out_mask[i] = 1 if overlap. Returns count of overlaps.
// N should be multiple of 8 for AVX2 path when available.
inline int aabb_overlap_soa_f32(
    float qminx, float qminy, float qminz,
    float qmaxx, float qmaxy, float qmaxz,
    const float *minx, const float *miny, const float *minz,
    const float *maxx, const float *maxy, const float *maxz,
    int n, int *out_idx, int out_cap)
{
    int w = 0;
#if defined(__AVX2__) || defined(_MSC_VER)
    const __m256 qminx8 = _mm256_set1_ps(qminx);
    const __m256 qminy8 = _mm256_set1_ps(qminy);
    const __m256 qminz8 = _mm256_set1_ps(qminz);
    const __m256 qmaxx8 = _mm256_set1_ps(qmaxx);
    const __m256 qmaxy8 = _mm256_set1_ps(qmaxy);
    const __m256 qmaxz8 = _mm256_set1_ps(qmaxz);

    int i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 mn_x = _mm256_loadu_ps(minx + i);
        __m256 mn_y = _mm256_loadu_ps(miny + i);
        __m256 mn_z = _mm256_loadu_ps(minz + i);
        __m256 mx_x = _mm256_loadu_ps(maxx + i);
        __m256 mx_y = _mm256_loadu_ps(maxy + i);
        __m256 mx_z = _mm256_loadu_ps(maxz + i);

        __m256 ok = _mm256_cmp_ps(qminx8, mx_x, _CMP_LE_OQ);
        ok = _mm256_and_ps(ok, _mm256_cmp_ps(mn_x, qmaxx8, _CMP_LE_OQ));
        ok = _mm256_and_ps(ok, _mm256_cmp_ps(qminy8, mx_y, _CMP_LE_OQ));
        ok = _mm256_and_ps(ok, _mm256_cmp_ps(mn_y, qmaxy8, _CMP_LE_OQ));
        ok = _mm256_and_ps(ok, _mm256_cmp_ps(qminz8, mx_z, _CMP_LE_OQ));
        ok = _mm256_and_ps(ok, _mm256_cmp_ps(mn_z, qmaxz8, _CMP_LE_OQ));

        int mask = _mm256_movemask_ps(ok);
        while (mask && w < out_cap) {
            unsigned long bit;
#if defined(_MSC_VER)
            _BitScanForward(&bit, (unsigned long)mask);
#else
            bit = (unsigned long)__builtin_ctz(mask);
#endif
            out_idx[w++] = i + (int)bit;
            mask &= mask - 1;
        }
    }
    for (; i < n && w < out_cap; ++i) {
        if (qminx <= maxx[i] && minx[i] <= qmaxx &&
            qminy <= maxy[i] && miny[i] <= qmaxy &&
            qminz <= maxz[i] && minz[i] <= qmaxz)
            out_idx[w++] = i;
    }
#else
    for (int i = 0; i < n && w < out_cap; ++i) {
        if (qminx <= maxx[i] && minx[i] <= qmaxx &&
            qminy <= maxy[i] && miny[i] <= qmaxy &&
            qminz <= maxz[i] && minz[i] <= qmaxz)
            out_idx[w++] = i;
    }
#endif
    return w;
}

} // namespace cork_simd

// Override hasIsct for BBox3d with SIMD path when included after bbox.h
#undef CORK_HAS_ISCT_BBOX3D_SIMD
#define CORK_HAS_ISCT_BBOX3D_SIMD 1

inline bool hasIsct_simd(const BBox3d &a, const BBox3d &b)
{
    return cork_simd::aabb_overlap(a, b);
}

#endif
