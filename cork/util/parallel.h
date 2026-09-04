// Thin parallel helpers: TBB when available, serial fallback otherwise.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

#if defined(CORK_USE_TBB)
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_invoke.h>
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#endif

namespace cork_par {

// Serial if n is below `threshold` (Manifold autoPolicy).  Default 0 = always
// honor grain / TBB.  Set CORK_AUTOPOLICY=1 to use 10000.
inline size_t auto_threshold() {
    static int t = -1;
    if (t < 0) t = std::getenv("CORK_AUTOPOLICY") ? 10000 : 0;
    return (size_t)t;
}

// f(begin, end) over [0, n) split in chunks of about `grain`
template<class F>
inline void for_range(size_t n, size_t grain, F &&f)
{
    if (n == 0) return;
    const size_t thr = auto_threshold();
    if (thr && n <= thr) { f((size_t)0, n); return; }
#if defined(CORK_USE_TBB)
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n, grain),
        [&](const tbb::blocked_range<size_t> &r) { f(r.begin(), r.end()); });
#else
    f((size_t)0, n);
#endif
}

// f(i) for each i in [0, n)
template<class F>
inline void for_each_idx(size_t n, size_t grain, F &&f)
{
    for_range(n, grain, [&](size_t b, size_t e) { for (size_t i = b; i < e; ++i) f(i); });
}

template<class It>
inline void sort(It b, It e)
{
#if defined(CORK_USE_TBB)
    tbb::parallel_sort(b, e);
#else
    std::sort(b, e);
#endif
}

template<class It, class Cmp>
inline void sort(It b, It e, Cmp cmp)
{
#if defined(CORK_USE_TBB)
    tbb::parallel_sort(b, e, cmp);
#else
    std::sort(b, e, cmp);
#endif
}

template<class... F>
inline void invoke(F &&...f)
{
#if defined(CORK_USE_TBB)
    tbb::parallel_invoke(std::forward<F>(f)...);
#else
    int dummy[] = { (f(), 0)... };
    (void)dummy;
#endif
}

// Per-thread accumulator with a combine step.
template<class T>
struct Local {
#if defined(CORK_USE_TBB)
    tbb::combinable<T> c;
    Local() {}
    template<class Init> explicit Local(Init init) : c(init) {}
    T &local() { return c.local(); }
    template<class F> void combine_each(F &&f) { c.combine_each(std::forward<F>(f)); }
#else
    T v;
    Local() {}
    template<class Init> explicit Local(Init init) : v(init()) {}
    T &local() { return v; }
    template<class F> void combine_each(F &&f) { f(v); }
#endif
};

// Build CSR (offsets, order) grouping m items by key(i) in [0, n).
// Items in each bucket are sorted ascending by item index (deterministic).
template<class KeyF>
inline void build_csr(size_t n, size_t m, KeyF key,
                      std::vector<unsigned> &offsets, std::vector<unsigned> &order)
{
    offsets.assign(n + 1, 0u);
    // histogram (serial is fine: streaming, ~ns per item)
    for (size_t i = 0; i < m; ++i) offsets[key(i) + 1]++;
    for (size_t k = 0; k < n; ++k) offsets[k + 1] += offsets[k];
    order.resize(m);
    std::vector<unsigned> cur(offsets.begin(), offsets.end() - 1);
    for (size_t i = 0; i < m; ++i) order[cur[key(i)]++] = (unsigned)i;
}

// Inclusive prefix of p[0..n): p[i] += p[i-1].  Blocked parallel for large n.
inline void prefix_sum_inplace(uint32_t *p, size_t n)
{
    if (n <= 1) return;
    if (n < 1 << 16) {
        for (size_t i = 1; i < n; ++i) p[i] += p[i - 1];
        return;
    }
    const size_t grain = 4096;
    const size_t nblk = (n + grain - 1) / grain;
    std::vector<uint32_t> blk(nblk);
    for_range(nblk, 1, [&](size_t b0, size_t b1) {
        for (size_t b = b0; b < b1; ++b) {
            size_t s = b * grain, e = std::min(n, s + grain);
            uint32_t acc = 0;
            for (size_t i = s; i < e; ++i) { acc += p[i]; p[i] = acc; }
            blk[b] = acc;
        }
    });
    for (size_t b = 1; b < nblk; ++b) blk[b] += blk[b - 1];
    for_range(nblk, 1, [&](size_t b0, size_t b1) {
        for (size_t b = b0; b < b1; ++b) {
            if (b == 0) continue;
            uint32_t add = blk[b - 1];
            size_t s = b * grain, e = std::min(n, s + grain);
            for (size_t i = s; i < e; ++i) p[i] += add;
        }
    });
}

// Stable LSD radix sort of n uint64 keys.  8-bit digits, ping-pong buffers.
inline void radix_sort_u64(uint64_t *a, size_t n)
{
    if (n < 2) return;
    if (n < 4096) { std::sort(a, a + n); return; }
    std::vector<uint64_t> tmp(n);
    uint64_t *src = a, *dst = tmp.data();
    const int BUCKETS = 256;
#if defined(CORK_USE_TBB)
    const size_t grain = 1 << 16;
    const size_t nblk = (n + grain - 1) / grain;
    std::vector<uint32_t> hist((size_t)nblk * BUCKETS);
    for (int pass = 0; pass < 8; ++pass) {
        const int shift = pass * 8;
        std::fill(hist.begin(), hist.end(), 0u);
        for_range(nblk, 1, [&](size_t b0, size_t b1) {
            for (size_t b = b0; b < b1; ++b) {
                size_t s = b * grain, e = std::min(n, s + grain);
                uint32_t *h = &hist[b * BUCKETS];
                for (size_t i = s; i < e; ++i)
                    h[(src[i] >> shift) & 0xffu]++;
            }
        });
        // exclusive prefix across (block, bucket) in bucket-major order
        uint32_t acc = 0;
        uint32_t off[BUCKETS];
        for (int k = 0; k < BUCKETS; ++k) {
            off[k] = acc;
            for (size_t b = 0; b < nblk; ++b) {
                uint32_t c = hist[b * BUCKETS + k];
                hist[b * BUCKETS + k] = acc;
                acc += c;
            }
        }
        (void)off;
        for_range(nblk, 1, [&](size_t b0, size_t b1) {
            for (size_t b = b0; b < b1; ++b) {
                size_t s = b * grain, e = std::min(n, s + grain);
                uint32_t *h = &hist[b * BUCKETS];
                for (size_t i = s; i < e; ++i) {
                    uint32_t d = (uint32_t)((src[i] >> shift) & 0xffu);
                    dst[h[d]++] = src[i];
                }
            }
        });
        std::swap(src, dst);
    }
    if (src != a) std::memcpy(a, src, n * sizeof(uint64_t));
#else
    uint32_t hist[BUCKETS];
    for (int pass = 0; pass < 8; ++pass) {
        const int shift = pass * 8;
        std::memset(hist, 0, sizeof(hist));
        for (size_t i = 0; i < n; ++i) hist[(src[i] >> shift) & 0xffu]++;
        uint32_t sum = 0;
        for (int k = 0; k < BUCKETS; ++k) { uint32_t c = hist[k]; hist[k] = sum; sum += c; }
        for (size_t i = 0; i < n; ++i) {
            uint32_t d = (uint32_t)((src[i] >> shift) & 0xffu);
            dst[hist[d]++] = src[i];
        }
        std::swap(src, dst);
    }
    if (src != a) std::memcpy(a, src, n * sizeof(uint64_t));
#endif
}

// Stable LSD radix sort of records by a uint64 key extractor.
template<class T, class KeyFn>
inline void radix_sort_by_key(T *a, size_t n, KeyFn keyfn)
{
    if (n < 2) return;
    if (n < 4096) {
        std::sort(a, a + n, [&](const T &x, const T &y) { return keyfn(x) < keyfn(y); });
        return;
    }
    std::vector<T> tmp(n);
    T *src = a, *dst = tmp.data();
    const int BUCKETS = 256;
#if defined(CORK_USE_TBB)
    const size_t grain = 1 << 15;
    const size_t nblk = (n + grain - 1) / grain;
    std::vector<uint32_t> hist((size_t)nblk * BUCKETS);
    for (int pass = 0; pass < 8; ++pass) {
        const int shift = pass * 8;
        std::fill(hist.begin(), hist.end(), 0u);
        for_range(nblk, 1, [&](size_t b0, size_t b1) {
            for (size_t b = b0; b < b1; ++b) {
                size_t s = b * grain, e = std::min(n, s + grain);
                uint32_t *h = &hist[b * BUCKETS];
                for (size_t i = s; i < e; ++i)
                    h[(keyfn(src[i]) >> shift) & 0xffu]++;
            }
        });
        uint32_t acc = 0;
        for (int k = 0; k < BUCKETS; ++k) {
            for (size_t b = 0; b < nblk; ++b) {
                uint32_t c = hist[b * BUCKETS + k];
                hist[b * BUCKETS + k] = acc;
                acc += c;
            }
        }
        for_range(nblk, 1, [&](size_t b0, size_t b1) {
            for (size_t b = b0; b < b1; ++b) {
                size_t s = b * grain, e = std::min(n, s + grain);
                uint32_t *h = &hist[b * BUCKETS];
                for (size_t i = s; i < e; ++i) {
                    uint32_t d = (uint32_t)((keyfn(src[i]) >> shift) & 0xffu);
                    dst[h[d]++] = src[i];
                }
            }
        });
        std::swap(src, dst);
    }
    if (src != a) std::memcpy(a, src, n * sizeof(T));
#else
    uint32_t hist[BUCKETS];
    for (int pass = 0; pass < 8; ++pass) {
        const int shift = pass * 8;
        std::memset(hist, 0, sizeof(hist));
        for (size_t i = 0; i < n; ++i) hist[(keyfn(src[i]) >> shift) & 0xffu]++;
        uint32_t sum = 0;
        for (int k = 0; k < BUCKETS; ++k) { uint32_t c = hist[k]; hist[k] = sum; sum += c; }
        for (size_t i = 0; i < n; ++i) {
            uint32_t d = (uint32_t)((keyfn(src[i]) >> shift) & 0xffu);
            dst[hist[d]++] = src[i];
        }
        std::swap(src, dst);
    }
    if (src != a) std::memcpy(a, src, n * sizeof(T));
#endif
}

} // namespace cork_par
