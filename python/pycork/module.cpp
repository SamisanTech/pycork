#include <string>

#include <Eigen/Eigen>

#include <cork/rawmesh/rawMesh.h>
#include <cork/mesh/mesh.h>
#include <cork/cork.h>
#include <cork/util/profile.h>
#include <cork/util/parallel.h>

#include <pybind11/pybind11.h>

#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>

#include <tuple>
#include <vector>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <algorithm>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace py = pybind11;

namespace pycork {

typedef Eigen::Matrix<double,Eigen::Dynamic,3> EigenVecX3d;
typedef Eigen::Matrix<uint64_t,Eigen::Dynamic,3> EigenVecX3i;
typedef std::tuple<EigenVecX3d, EigenVecX3i> MeshTuple;

bool isSolid(const EigenVecX3d &verts,
             const EigenVecX3i &tris) {

    CorkMesh mesh;

    eigenToCorkMesh(verts, tris, &mesh);

    bool result = true;

    if(mesh.isSelfIntersecting())
        result = false;

    if(!mesh.isClosed())
        result = false;

    return result;

}


MeshTuple booleanUnion(const EigenVecX3d &vertsA,
                       const EigenVecX3i &trisA,
                       const EigenVecX3d &vertsB,
                       const EigenVecX3i &trisB) {

    CorkMesh meshA, meshB;

    eigenToCorkMesh(vertsA, trisA, &meshA);
    eigenToCorkMesh(vertsB, trisB, &meshB);

    MeshTuple meshOut;

    meshA.boolUnion(meshB);

    corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

    return meshOut;
}

MeshTuple booleanDifference(const EigenVecX3d &vertsA,
                            const EigenVecX3i &trisA,
                            const EigenVecX3d &vertsB,
                            const EigenVecX3i &trisB) {

    CorkMesh meshA, meshB;

    eigenToCorkMesh(vertsA, trisA, &meshA);
    eigenToCorkMesh(vertsB, trisB, &meshB);

    MeshTuple meshOut;

    meshA.boolDiff(meshB);

    corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

    return meshOut;
}

MeshTuple booleanIntersection(const EigenVecX3d &vertsA,
                              const EigenVecX3i &trisA,
                              const EigenVecX3d &vertsB,
                              const EigenVecX3i &trisB) {

    CorkMesh meshA, meshB;

    eigenToCorkMesh(vertsA, trisA, &meshA);
    eigenToCorkMesh(vertsB, trisB, &meshB);

    MeshTuple meshOut;

    meshA.boolIsct(meshB);

    corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

    return meshOut;
}


MeshTuple booleanXor(const EigenVecX3d &vertsA,
                     const EigenVecX3i &trisA,
                     const EigenVecX3d &vertsB,
                     const EigenVecX3i &trisB) {

    CorkMesh meshA, meshB;

    eigenToCorkMesh(vertsA, trisA, &meshA);
    eigenToCorkMesh(vertsB, trisB, &meshB);

    MeshTuple meshOut;

    meshA.boolXor(meshB);

    corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

    return meshOut;
}


MeshTuple resolveIntersection(const EigenVecX3d &vertsA,
                              const EigenVecX3i &trisA) {
    CORK_PROF("py.resolveIntersection total");
    CorkMesh meshA;

    {
        CORK_PROF("py.eigenToCorkMesh");
        eigenToCorkMesh(vertsA, trisA, &meshA);
    }

    MeshTuple meshOut;

    {
        CORK_PROF("py.Mesh::resolveIntersections");
        meshA.resolveIntersections();
    }

    {
        CORK_PROF("py.corkMesh2Eigen");
        corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));
    }

    return meshOut;
}

// resolve self-intersections, then keep only the outer hull
// (faces with generalized winding number 0 on their outside).
// Returns (verts, tris, stats) with stats = {patches, kept, flipped, deleted, unresolved, rays}
std::tuple<EigenVecX3d, EigenVecX3i, py::dict> outerHull(const EigenVecX3d &vertsA,
                                                          const EigenVecX3i &trisA,
                                                          int raysPerPatch,
                                                          bool resolve) {
    CORK_PROF("py.outerHull total");
    CorkMesh meshA;
    {
        CORK_PROF("py.eigenToCorkMesh");
        eigenToCorkMesh(vertsA, trisA, &meshA);
    }
    if (resolve) {
        CORK_PROF("py.Mesh::resolveIntersections");
        meshA.resolveIntersections();
    }
    cork_hull::HullStats st;
    {
        CORK_PROF("py.Mesh::outerHull");
        meshA.outerHull(raysPerPatch, &st);
    }
    std::tuple<EigenVecX3d, EigenVecX3i, py::dict> out;
    {
        CORK_PROF("py.corkMesh2Eigen");
        corkMesh2Eigen(meshA, std::get<0>(out), std::get<1>(out));
    }
    py::dict d;
    d["patches"] = st.patches;
    d["kept"] = st.kept;
    d["flipped"] = st.flipped;
    d["deleted"] = st.deleted;
    d["unresolved_patches"] = st.unresolved;
    d["split_patches"] = st.split_patches;
    d["seeds"] = st.seeds;
    d["sliver_faces_pruned"] = st.pruned;
    d["rays"] = st.rays;
    std::get<2>(out) = d;
    return out;
}

#ifdef _WIN32
static std::wstring utf8_to_wide(const std::string &s) {
    if (s.empty()) return std::wstring();
    int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), nullptr, 0);
    std::wstring w((size_t)n, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), &w[0], n);
    return w;
}
#endif

#pragma pack(push, 1)
struct StlTriRec {
    float n[3];
    float v[3][3];
    uint16_t attr;
};
#pragma pack(pop)
static_assert(sizeof(StlTriRec) == 50, "binary STL triangle is 50 bytes");

static uint32_t hash3(uint32_t a, uint32_t b, uint32_t c) {
    uint32_t h = a * 0x9e3779b1u;
    h ^= b + 0x9e3779b1u + (h << 6) + (h >> 2);
    h ^= c + 0x9e3779b1u + (h << 6) + (h >> 2);
    return h;
}

static const char *map_read(const std::string &path, size_t &n, void **keep) {
#ifdef _WIN32
    std::wstring w = utf8_to_wide(path);
    HANDLE f = CreateFileW(w.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                           OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, nullptr);
    if (f == INVALID_HANDLE_VALUE)
        throw std::runtime_error("readSTL: cannot open " + path);
    LARGE_INTEGER sz;
    GetFileSizeEx(f, &sz);
    n = (size_t)sz.QuadPart;
    HANDLE m = CreateFileMappingW(f, nullptr, PAGE_READONLY, 0, 0, nullptr);
    const char *p = m ? (const char *)MapViewOfFile(m, FILE_MAP_READ, 0, 0, 0) : nullptr;
    CloseHandle(f);
    if (!p) {
        if (m) CloseHandle(m);
        throw std::runtime_error("readSTL: map failed " + path);
    }
    keep[0] = m;
    keep[1] = (void *)p;
    return p;
#else
    (void)keep;
    throw std::runtime_error("readSTL: mmap not implemented");
#endif
}

static char *map_write(const std::string &path, size_t n, void **keep) {
#ifdef _WIN32
    std::wstring w = utf8_to_wide(path);
    HANDLE f = CreateFileW(w.c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr,
                           CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (f == INVALID_HANDLE_VALUE)
        throw std::runtime_error("writeSTL: cannot open " + path);
    LARGE_INTEGER sz;
    sz.QuadPart = (LONGLONG)n;
    SetFilePointerEx(f, sz, nullptr, FILE_BEGIN);
    SetEndOfFile(f);
    HANDLE m = CreateFileMappingW(f, nullptr, PAGE_READWRITE, sz.HighPart, sz.LowPart, nullptr);
    char *p = m ? (char *)MapViewOfFile(m, FILE_MAP_WRITE, 0, 0, n) : nullptr;
    CloseHandle(f);
    if (!p) {
        if (m) CloseHandle(m);
        throw std::runtime_error("writeSTL: map failed " + path);
    }
    keep[0] = m;
    keep[1] = p;
    return p;
#else
    (void)path; (void)n; (void)keep;
    throw std::runtime_error("writeSTL: mmap not implemented");
#endif
}

static void unmap_keep(void **keep) {
#ifdef _WIN32
    if (keep[1]) UnmapViewOfFile(keep[1]);
    if (keep[0]) CloseHandle((HANDLE)keep[0]);
#endif
    keep[0] = keep[1] = nullptr;
}

std::tuple<py::array_t<double>, py::array_t<uint64_t>> readSTL(const std::string &path) {
    CORK_PROF("py.readSTL");
    size_t nbytes = 0;
    void *keep[2] = {nullptr, nullptr};
    const char *buf = map_read(path, nbytes, keep);
    if (nbytes < 84) {
        unmap_keep(keep);
        throw std::runtime_error("readSTL: file too small " + path);
    }
    uint32_t ntri = 0;
    std::memcpy(&ntri, buf + 80, 4);
    const size_t need = 84ull + (size_t)ntri * 50ull;
    if (nbytes < need) {
        unmap_keep(keep);
        throw std::runtime_error("readSTL: truncated binary STL " + path);
    }

    const StlTriRec *rec = reinterpret_cast<const StlTriRec *>(buf + 84);
    const size_t ncorner = (size_t)ntri * 3;

    // 8-byte slots (hash + id). Compare coordinates in the packed xyz
    // array so the table stays in L3 (21.stl: ~2M slots, 16 MB).
    size_t cap = 1;
    while (cap < ((size_t)ntri < 8 ? 8 : (size_t)ntri)) cap <<= 1;
    if (cap < 1024) cap = 1024;
    struct Slot { uint32_t h, id; };
    std::vector<Slot> tab(cap);
    std::memset(tab.data(), 0xff, cap * sizeof(Slot));
    uint32_t mask = (uint32_t)cap - 1;

    std::vector<uint32_t> xyz;
    xyz.reserve(((ncorner >> 1) + 16) * 3);
    auto F = py::array_t<uint64_t>({ (py::ssize_t)ntri, (py::ssize_t)3 });
    uint64_t *fd = F.mutable_data();

    uint32_t next = 0;
    auto grow = [&]() {
        size_t ncap = cap << 1;
        std::vector<Slot> ntab(ncap);
        std::memset(ntab.data(), 0xff, ncap * sizeof(Slot));
        const uint32_t nmask = (uint32_t)ncap - 1;
        for (size_t i = 0; i < cap; ++i) {
            if (tab[i].id == 0xffffffffu) continue;
            uint32_t j = tab[i].h & nmask;
            while (ntab[j].id != 0xffffffffu) j = (j + 1) & nmask;
            ntab[j] = tab[i];
        }
        tab.swap(ntab);
        cap = ncap;
        mask = nmask;
    };

    for (uint32_t t = 0; t < ntri; ++t) {
        for (int k = 0; k < 3; ++k) {
            uint32_t b0, b1, b2;
            std::memcpy(&b0, &rec[t].v[k][0], 4);
            std::memcpy(&b1, &rec[t].v[k][1], 4);
            std::memcpy(&b2, &rec[t].v[k][2], 4);
            const uint32_t h = hash3(b0, b1, b2);
            uint32_t i = h & mask;
            for (;;) {
                Slot &s = tab[i];
                if (s.id == 0xffffffffu) {
                    if (next * 4 > (uint32_t)cap * 3) { grow(); i = h & mask; continue; }
                    s.h = h; s.id = next;
                    xyz.push_back(b0); xyz.push_back(b1); xyz.push_back(b2);
                    fd[t * 3 + k] = next++;
                    break;
                }
                if (s.h == h) {
                    const uint32_t *p = &xyz[s.id * 3];
                    if (p[0] == b0 && p[1] == b1 && p[2] == b2) {
                        fd[t * 3 + k] = s.id;
                        break;
                    }
                }
                i = (i + 1) & mask;
            }
        }
    }
    unmap_keep(keep);

    auto V = py::array_t<double>({ (py::ssize_t)next, (py::ssize_t)3 });
    double *vd = V.mutable_data();
    const float *fp = reinterpret_cast<const float *>(xyz.data());
    for (uint32_t i = 0, n = next * 3; i < n; ++i)
        vd[i] = fp[i];
    return {V, F};
}

void writeSTL(const std::string &path,
              py::array_t<double, py::array::c_style | py::array::forcecast> V,
              py::array_t<uint64_t, py::array::c_style | py::array::forcecast> F) {
    CORK_PROF("py.writeSTL");
    if (V.ndim() != 2 || V.shape(1) != 3 || F.ndim() != 2 || F.shape(1) != 3)
        throw std::runtime_error("writeSTL: verts and tris must be (N,3)");
    const uint32_t ntri = (uint32_t)F.shape(0);
    const size_t nbytes = 84ull + (size_t)ntri * 50ull;
    void *keep[2] = {nullptr, nullptr};
    char *buf = map_write(path, nbytes, keep);
    std::memset(buf, 0, 84);
    std::memcpy(buf + 80, &ntri, 4);
    StlTriRec *rec = reinterpret_cast<StlTriRec *>(buf + 84);
    const double *vp = V.data();
    const uint64_t *fp = F.data();
    cork_par::for_range((size_t)ntri, 4096, [&](size_t b, size_t e) {
        for (size_t t = b; t < e; ++t) {
            const uint64_t a = fp[t * 3 + 0];
            const uint64_t bb = fp[t * 3 + 1];
            const uint64_t c = fp[t * 3 + 2];
            rec[t].n[0] = rec[t].n[1] = rec[t].n[2] = 0.f;
            rec[t].attr = 0;
            rec[t].v[0][0] = (float)vp[a * 3 + 0];
            rec[t].v[0][1] = (float)vp[a * 3 + 1];
            rec[t].v[0][2] = (float)vp[a * 3 + 2];
            rec[t].v[1][0] = (float)vp[bb * 3 + 0];
            rec[t].v[1][1] = (float)vp[bb * 3 + 1];
            rec[t].v[1][2] = (float)vp[bb * 3 + 2];
            rec[t].v[2][0] = (float)vp[c * 3 + 0];
            rec[t].v[2][1] = (float)vp[c * 3 + 1];
            rec[t].v[2][2] = (float)vp[c * 3 + 2];
        }
    });
    unmap_keep(keep);
}

} // end of namespace



PYBIND11_MODULE(pycork, m) {

    m.doc() = R"pbdoc(
        Pycork Module
        -----------------------
        .. currentmodule:: pycork
        .. autosummary::
           :toctree: _generate

    )pbdoc";

    m.def("isSolid", &pycork::isSolid, "Determines if the mesh is manifold",
                     py::arg("vertices"), py::arg("tris"))
     .def("union", &pycork::booleanUnion, "Computes boolean Union between two meshes",
                    py::arg("vertsA"), py::arg("trisA"),
                    py::arg("vertsB"), py::arg("trisB"))
     .def("difference", &pycork::booleanDifference, "Computes boolean differnece between two meshes",
                   py::arg("vertsA"), py::arg("trisA"),
                   py::arg("vertsB"), py::arg("trisB"))
     .def("intersection", &pycork::booleanIntersection, "Computes boolean intersection between two meshes",
                   py::arg("vertsA"), py::arg("trisA"),
                   py::arg("vertsB"), py::arg("trisB"))
     .def("intersection", &pycork::booleanIntersection, "Computes boolean xor between two meshes",
                           py::arg("vertsA"), py::arg("trisA"),
                           py::arg("vertsB"), py::arg("trisB"))
     .def("resolveIntersection", &pycork::resolveIntersection, "Computes the intersection between two meshes",
                                  py::arg("vertsA"), py::arg("trisA"))
     .def("outerHull", &pycork::outerHull,
          "Resolves self-intersections and keeps only the outer hull (winding-number 0 side). "
          "Returns (verts, tris, stats)",
          py::arg("verts"), py::arg("tris"), py::arg("raysPerPatch") = 5, py::arg("resolve") = true)
     .def("readSTL", &pycork::readSTL,
          "Fast binary STL load + exact-float vertex weld. Returns (verts, tris).",
          py::arg("path"))
     .def("writeSTL", &pycork::writeSTL,
          "Fast binary STL write.",
          py::arg("path"), py::arg("verts"), py::arg("tris"));


#ifdef PROJECT_VERSION
    m.attr("__version__") = "PROJECT_VERSION";
#else
    m.attr("__version__") = "dev";
#endif

}

