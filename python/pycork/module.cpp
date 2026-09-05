#include <string>

#include <Eigen/Eigen>

#include <cork/rawmesh/rawMesh.h>
#include <cork/mesh/mesh.h>
#include <cork/cork.h>
#include <cork/core/surface.h>
#include <cork/io/stl.h>
#include <cork/repair/options.h>
#include <cork/repair/pipeline.h>
#include <cork/util/profile.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>

namespace py = pybind11;

namespace pycork {

typedef Eigen::Matrix<double,Eigen::Dynamic,3> EigenVecX3d;
typedef Eigen::Matrix<uint64_t,Eigen::Dynamic,3> EigenVecX3i;
typedef std::tuple<EigenVecX3d, EigenVecX3i> MeshTuple;

static py::object g_manifold;

void set_manifold(py::object mod) { g_manifold = std::move(mod); }

static bool manifold_hull_requested(const std::string &backend) {
    return backend == "manifold" || backend == "exact" || backend == "union_all";
}

static bool has_manifold() {
    return bool(g_manifold) && !g_manifold.is_none() && py::hasattr(g_manifold, "union_all");
}

// Exact Manifold outer solid: Decompose + BatchBoolean Add (not convex Hull).
static bool apply_manifold_union_all(CorkMesh &mesh) {
    if (!has_manifold())
        return false;
    EigenVecX3d v;
    EigenVecX3i f;
    corkMesh2Eigen(mesh, v, f);
    py::object out = g_manifold.attr("union_all")(v, f);
    py::sequence t = out;
    EigenVecX3d v2 = t[0].cast<EigenVecX3d>();
    typedef Eigen::Matrix<int64_t, Eigen::Dynamic, 3> EigenVecX3i64;
    EigenVecX3i64 fi = t[1].cast<EigenVecX3i64>();
    EigenVecX3i f2 = fi.cast<uint64_t>();
    if (v2.rows() == 0 || f2.rows() == 0)
        return false;
    eigenToCorkMesh(v2, f2, &mesh);
    return true;
}

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
                                                          bool resolve,
                                                          const std::string &hull_backend,
                                                          bool exact) {
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
    std::string used = "ours";
    const bool wantManifold = !exact && manifold_hull_requested(hull_backend);
    if (wantManifold) {
        CORK_PROF("py.outerHull.manifold");
        bool ok = false;
        try {
            ok = apply_manifold_union_all(meshA);
        } catch (const py::error_already_set &) {
            ok = false;
        } catch (const std::exception &) {
            ok = false;
        }
        if (ok) {
            used = "manifold";
        } else {
            meshA.outerHull(raysPerPatch, &st);
            used = "ours";
        }
    } else if (exact) {
        CORK_PROF("py.Mesh::outerHullExact");
        meshA.outerHullExact(&st);
        used = "exact";
    } else {
        CORK_PROF("py.Mesh::outerHull");
        meshA.outerHull(raysPerPatch, &st);
    }
    std::tuple<EigenVecX3d, EigenVecX3i, py::dict> out;
    {
        CORK_PROF("py.corkMesh2Eigen");
        corkMesh2Eigen(meshA, std::get<0>(out), std::get<1>(out));
    }
    py::dict d;
    d["hull_backend"] = used;
    d["exactHull"] = exact;
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

template <class S>
std::tuple<py::array_t<double>, py::array_t<uint64_t>> surfaceToNumpy(const cork::Surface<S> &s) {
    auto V = py::array_t<double>({ (py::ssize_t)s.nV(), (py::ssize_t)3 });
    auto F = py::array_t<uint64_t>({ (py::ssize_t)s.nF(), (py::ssize_t)3 });
    double *vd = V.mutable_data();
    uint64_t *fd = F.mutable_data();
    for (size_t i = 0, n = s.xyz.size(); i < n; ++i)
        vd[i] = static_cast<double>(s.xyz[i]);
    for (size_t i = 0, n = s.f.size(); i < n; ++i)
        fd[i] = s.f[i];
    return {V, F};
}

std::tuple<py::array_t<double>, py::array_t<uint64_t>> readSTL(const std::string &path) {
    CORK_PROF("py.readSTL");
    return surfaceToNumpy(cork::io::read_stl(path));
}

void writeSTL(const std::string &path,
              py::array_t<double, py::array::c_style | py::array::forcecast> V,
              py::array_t<uint64_t, py::array::c_style | py::array::forcecast> F) {
    CORK_PROF("py.writeSTL");
    if (V.ndim() != 2 || V.shape(1) != 3 || F.ndim() != 2 || F.shape(1) != 3)
        throw std::runtime_error("writeSTL: verts and tris must be (N,3)");
    cork::io::write_stl(path, V.data(), F.data(), (std::uint32_t)F.shape(0));
}

// One CorkMesh build, then the optional repair stages. Existing
// resolveIntersection / outerHull stay as the thin one-kernel bindings.
std::tuple<EigenVecX3d, EigenVecX3i, py::dict> repair(const EigenVecX3d &verts,
                                                     const EigenVecX3i &tris,
                                                     bool resolve,
                                                     bool hull,
                                                     bool unify,
                                                     bool cluster,
                                                     bool puzzle,
                                                     bool clean,
                                                     bool perturb,
                                                     bool noise,
                                                     bool collapse,
                                                     bool si_subset,
                                                     int raysPerPatch,
                                                     int minFaces,
                                                     int hardDegree,
                                                     double perturbIntensity,
                                                     double collapseRel,
                                                     const std::string &hull_backend,
                                                     bool exactHull) {
    CORK_PROF("py.repair total");
    CorkMesh mesh;
    {
        CORK_PROF("py.eigenToCorkMesh");
        eigenToCorkMesh(verts, tris, &mesh);
    }
    const bool wantManifold = hull && !exactHull && manifold_hull_requested(hull_backend);
    cork::repair::Options opt;
    opt.resolve = resolve;
    opt.hull = hull;
    opt.exactHull = exactHull;
    opt.deferHull = wantManifold;
    opt.unify = unify;
    opt.cluster = cluster;
    opt.puzzle = puzzle;
    opt.clean = clean;
    opt.perturb = perturb;
    opt.noise = noise;
    opt.collapse = collapse;
    opt.si_subset = si_subset;
    opt.fill = !wantManifold;
    opt.raysPerPatch = raysPerPatch;
    opt.minFaces = minFaces;
    opt.hardDegree = hardDegree;
    opt.perturbIntensity = perturbIntensity;
    opt.collapseRel = collapseRel;
    cork::repair::Stats st;
    cork::repair::Pipeline::run(mesh, opt, &st);

    std::string used = !hull ? "off" : (exactHull ? "exact" : "ours");
    bool fallback = false;
    if (wantManifold) {
        CORK_PROF("repair.hull.manifold");
        bool ok = false;
        try {
            ok = apply_manifold_union_all(mesh);
        } catch (const py::error_already_set &) {
            ok = false;
        } catch (const std::exception &) {
            ok = false;
        }
        if (ok) {
            used = "manifold";
        } else {
            cork_hull::HullStats hs;
            cork::repair::extract_hull(mesh, opt, &hs);
            used = opt.exactHull ? "exact" : "ours";
            fallback = true;
        }
    }

    std::tuple<EigenVecX3d, EigenVecX3i, py::dict> out;
    {
        CORK_PROF("py.corkMesh2Eigen");
        corkMesh2Eigen(mesh, std::get<0>(out), std::get<1>(out));
    }
    py::dict d;
    d["resolve"] = resolve;
    d["hull"] = hull;
    d["hull_backend"] = used;
    d["exactHull"] = exactHull;
    d["hull_fallback"] = fallback;
    d["unify"] = unify;
    d["cluster"] = cluster;
    d["puzzle"] = puzzle;
    d["clean"] = clean;
    d["perturb"] = perturb;
    d["noise"] = noise;
    d["collapse"] = collapse;
    d["si_subset"] = si_subset;
    d["raysPerPatch"] = raysPerPatch;
    d["minFaces"] = minFaces;
    d["shells"] = st.shells;
    d["clusters"] = st.clusters;
    d["closed"] = st.closed;
    d["open"] = st.open;
    d["hard"] = st.hard;
    d["unified"] = st.unified;
    d["puzzled"] = st.puzzled;
    std::get<2>(out) = d;
    return out;
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
     .def("set_manifold", &pycork::set_manifold,
          "Inject Exact Manifold module (must expose union_all). Used when hull_backend='manifold'.",
          py::arg("mod"))
     .def("has_manifold", &pycork::has_manifold,
          "True if set_manifold() was given a module with union_all.")
     .def("outerHull", &pycork::outerHull,
          "Resolves self-intersections and keeps only the outer hull. "
          "hull_backend: 'ours' (winding) or 'manifold' (Exact Decompose+BatchBoolean Add). "
          "Returns (verts, tris, stats)",
          py::arg("verts"), py::arg("tris"), py::arg("raysPerPatch") = 5, py::arg("resolve") = true,
          py::arg("hull_backend") = "ours", py::arg("exact") = false)
     .def("readSTL", &pycork::readSTL,
          "Fast binary STL load + exact-float vertex weld. Returns (verts, tris).",
          py::arg("path"))
     .def("writeSTL", &pycork::writeSTL,
          "Fast binary STL write. verts any float dtype, tris any integer dtype.",
          py::arg("path"), py::arg("verts"), py::arg("tris"))
     .def("repair", &pycork::repair,
          "In-place repair pipeline (one CorkMesh). Extra Mira stages default off.",
          py::arg("verts"), py::arg("tris"),
          py::arg("resolve") = true,
          py::arg("hull") = true,
          py::arg("unify") = false,
          py::arg("cluster") = false,
          py::arg("puzzle") = false,
          py::arg("clean") = false,
          py::arg("perturb") = false,
          py::arg("noise") = true,
          py::arg("collapse") = false,
          py::arg("si_subset") = false,
          py::arg("raysPerPatch") = 5,
          py::arg("minFaces") = 5,
          py::arg("hardDegree") = 30,
          py::arg("perturbIntensity") = 1e-3,
          py::arg("collapseRel") = 0.02,
          py::arg("hull_backend") = "ours",
          py::arg("exactHull") = false);


#ifdef PROJECT_VERSION
    m.attr("__version__") = "PROJECT_VERSION";
#else
    m.attr("__version__") = "dev";
#endif

}

