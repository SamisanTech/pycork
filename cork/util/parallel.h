// Thin parallel helpers: TBB when available, serial fallback otherwise.
#pragma once

#include <algorithm>
#include <cstddef>
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

// f(begin, end) over [0, n) split in chunks of about `grain`
template<class F>
inline void for_range(size_t n, size_t grain, F &&f)
{
#if defined(CORK_USE_TBB)
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n, grain),
        [&](const tbb::blocked_range<size_t> &r) { f(r.begin(), r.end()); });
#else
    if (n) f((size_t)0, n);
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

} // namespace cork_par
