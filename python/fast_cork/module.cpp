#include <cstdio>
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
#include <cstring>
#include <stdexcept>
#include <tuple>

#include "seh_guard.h"

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
// The body runs under an SEH guard: Triangle's exit() is now an exception,
// and an access violation returns si_aborted instead of killing Mira.
struct RepairJob {
    const EigenVecX3d *verts;
    const EigenVecX3i *tris;
    bool resolve, hull, unify, cluster, puzzle, clean, perturb, noise, collapse, si_subset, exactHull;
    int raysPerPatch, minFaces, hardDegree;
    double perturbIntensity, collapseRel;
    const std::string *hull_backend;
    EigenVecX3d *outV;
    EigenVecX3i *outF;
    bool produced;
    bool aborted;
    long long pairs;
    long long hits;
    std::string used;
    bool fallback;
    cork::repair::Stats st;
    char error[160];
};

static void repair_invoke(void *raw) {
    auto *job = static_cast<RepairJob *>(raw);
    job->produced = false;
    job->aborted = false;
    job->pairs = 0;
    job->hits = 0;
    job->fallback = false;
    job->outV = nullptr;
    job->outF = nullptr;
    job->error[0] = 0;
    try {
        job->used = !job->hull ? "off" : (job->exactHull ? "exact" : "ours");
        CorkMesh mesh;
        {
            CORK_PROF("py.eigenToCorkMesh");
            eigenToCorkMesh(*job->verts, *job->tris, &mesh);
        }
        const bool wantManifold = job->hull && !job->exactHull && manifold_hull_requested(*job->hull_backend);
        cork::repair::Options opt;
        opt.resolve = job->resolve;
        opt.hull = job->hull;
        opt.exactHull = job->exactHull;
        opt.deferHull = wantManifold;
        opt.unify = job->unify;
        opt.cluster = job->cluster;
        opt.puzzle = job->puzzle;
        opt.clean = job->clean;
        opt.perturb = job->perturb;
        opt.noise = job->noise;
        opt.collapse = job->collapse;
        opt.si_subset = job->si_subset;
        opt.fill = !wantManifold;
        opt.raysPerPatch = job->raysPerPatch;
        opt.minFaces = job->minFaces;
        opt.hardDegree = job->hardDegree;
        opt.perturbIntensity = job->perturbIntensity;
        opt.collapseRel = job->collapseRel;
        cork::repair::Pipeline::run(mesh, opt, &job->st);
        if (mesh.siUnbounded) {
            job->aborted = true;
            job->pairs = mesh.siAbortPairs;
            job->hits = mesh.siAbortHits;
            job->produced = true;
            return;
        }
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
                job->used = "manifold";
            } else {
                cork_hull::HullStats hs;
                cork::repair::extract_hull(mesh, opt, &hs);
                job->used = opt.exactHull ? "exact" : "ours";
                job->fallback = true;
            }
        }
        auto *ov = new EigenVecX3d();
        auto *of = new EigenVecX3i();
        {
            CORK_PROF("py.corkMesh2Eigen");
            corkMesh2Eigen(mesh, *ov, *of);
        }
        job->outV = ov;
        job->outF = of;
        job->produced = true;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "fast_cork: repair exception: %s\n", e.what());
        std::fflush(stderr);
        std::snprintf(job->error, sizeof(job->error), "%s", e.what());
        job->aborted = true;
        job->produced = true;
    } catch (...) {
        std::fprintf(stderr, "fast_cork: repair unknown exception\n");
        std::fflush(stderr);
        std::snprintf(job->error, sizeof(job->error), "unknown");
        job->aborted = true;
        job->produced = true;
    }
}

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
    std::fprintf(stderr, "fast_cork: repair enter\n");
    std::fflush(stderr);
    RepairJob job;
    job.verts = &verts;
    job.tris = &tris;
    job.resolve = resolve;
    job.hull = hull;
    job.unify = unify;
    job.cluster = cluster;
    job.puzzle = puzzle;
    job.clean = clean;
    job.perturb = perturb;
    job.noise = noise;
    job.collapse = collapse;
    job.si_subset = si_subset;
    job.exactHull = exactHull;
    job.raysPerPatch = raysPerPatch;
    job.minFaces = minFaces;
    job.hardDegree = hardDegree;
    job.perturbIntensity = perturbIntensity;
    job.collapseRel = collapseRel;
    job.hull_backend = &hull_backend;
    job.outV = nullptr;
    job.outF = nullptr;
    job.produced = false;
    job.aborted = false;
    job.pairs = 0;
    job.hits = 0;
    job.fallback = false;

    const int native = cork_seh_invoke(repair_invoke, &job);
    const bool crashed = native != 0 || !job.produced;
    if (crashed) {
        std::fprintf(stderr, "fast_cork: repair aborted (native 0x%08X); returning unchanged\n",
                     (unsigned)native);
        std::fflush(stderr);
        delete job.outV;
        delete job.outF;
        job.outV = nullptr;
        job.outF = nullptr;
        job.aborted = true;
        job.pairs = 0;
        job.hits = 0;
    }

    std::tuple<EigenVecX3d, EigenVecX3i, py::dict> out;
    if (job.outV && job.outF) {
        std::get<0>(out) = std::move(*job.outV);
        std::get<1>(out) = std::move(*job.outF);
        delete job.outV;
        delete job.outF;
    }
    py::dict d;
    d["resolve"] = resolve;
    d["hull"] = hull;
    d["hull_backend"] = job.used;
    d["exactHull"] = exactHull;
    d["hull_fallback"] = job.fallback;
    d["unify"] = unify;
    d["cluster"] = cluster;
    d["puzzle"] = puzzle;
    d["clean"] = clean;
    d["perturb"] = perturb;
    d["noise"] = noise;
    d["collapse"] = collapse;
    d["si_subset"] = si_subset;
    d["si_aborted"] = job.aborted || crashed;
    d["si_cross_pairs"] = job.pairs;
    d["si_hits"] = job.hits;
    d["si_native"] = (unsigned)(crashed ? native : 0);
    d["raysPerPatch"] = raysPerPatch;
    d["minFaces"] = minFaces;
    d["shells"] = job.st.shells;
    d["clusters"] = job.st.clusters;
    d["closed"] = job.st.closed;
    d["open"] = job.st.open;
    d["hard"] = job.st.hard;
    d["unified"] = job.st.unified;
    d["puzzled"] = job.st.puzzled;
    if (job.error[0]) d["si_error"] = job.error;
    std::get<2>(out) = d;
    return out;
}

} // end of namespace



PYBIND11_MODULE(fast_cork, m) {

    m.doc() = R"pbdoc(
        Pycork Module
        -----------------------
        .. currentmodule:: fast_cork
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

