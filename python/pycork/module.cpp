#include <string>

#include <Eigen/Eigen>

#include <cork/rawmesh/rawMesh.h>
#include <cork/mesh/mesh.h>
#include <cork/cork.h>
#include <cork/util/profile.h>

#include <pybind11/pybind11.h>

#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>

#include <tuple>

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
          py::arg("verts"), py::arg("tris"), py::arg("raysPerPatch") = 5, py::arg("resolve") = true);


#ifdef PROJECT_VERSION
    m.attr("__version__") = "PROJECT_VERSION";
#else
    m.attr("__version__") = "dev";
#endif

}

