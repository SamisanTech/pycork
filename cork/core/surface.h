// Templated triangle soup. File / Python dtypes live here.
// CorkMesh kernels stay double (exact predicates). Cast once at the boundary.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace cork {

template <class S>
struct Surface {
    using scalar_type = S;
    using index_type = std::uint32_t;

    std::vector<S> xyz;          // 3 * nV, xyzxyz...
    std::vector<index_type> f;   // 3 * nF

    std::size_t nV() const { return xyz.size() / 3; }
    std::size_t nF() const { return f.size() / 3; }

    template <class D>
    Surface<D> cast() const {
        Surface<D> o;
        o.xyz.resize(xyz.size());
        o.f = f;
        for (std::size_t i = 0; i < xyz.size(); ++i)
            o.xyz[i] = static_cast<D>(xyz[i]);
        return o;
    }
};

}  // namespace cork
