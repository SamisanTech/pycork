#include <string>

#include <Eigen/Eigen>

#include <cork/rawmesh/rawMesh.h>
#include <cork/mesh/mesh.h>
#include <cork/cork.h>
#include <cork/cancel.h>

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

[[noreturn]] void throw_cancelled()
{
    throw std::runtime_error("cancelled");
}

bool isSolid(const EigenVecX3d &verts,
             const EigenVecX3i &tris) {

    cork_reset_cancel();
    try {
        py::gil_scoped_release release;
        CorkMesh mesh;

        eigenToCorkMesh(verts, tris, &mesh);

        bool result = true;

        if(mesh.isSelfIntersecting())
            result = false;

        if(!mesh.isClosed())
            result = false;

        return result;
    } catch (const CorkCancelled &) {
        throw_cancelled();
    }
}


MeshTuple booleanUnion(const EigenVecX3d &vertsA,
                       const EigenVecX3i &trisA,
                       const EigenVecX3d &vertsB,
                       const EigenVecX3i &trisB) {

    cork_reset_cancel();
    try {
        py::gil_scoped_release release;
        CorkMesh meshA, meshB;

        eigenToCorkMesh(vertsA, trisA, &meshA);
        eigenToCorkMesh(vertsB, trisB, &meshB);

        MeshTuple meshOut;

        meshA.boolUnion(meshB);

        corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

        return meshOut;
    } catch (const CorkCancelled &) {
        throw_cancelled();
    }
}

MeshTuple booleanDifference(const EigenVecX3d &vertsA,
                            const EigenVecX3i &trisA,
                            const EigenVecX3d &vertsB,
                            const EigenVecX3i &trisB) {

    cork_reset_cancel();
    try {
        py::gil_scoped_release release;
        CorkMesh meshA, meshB;

        eigenToCorkMesh(vertsA, trisA, &meshA);
        eigenToCorkMesh(vertsB, trisB, &meshB);

        MeshTuple meshOut;

        meshA.boolDiff(meshB);

        corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

        return meshOut;
    } catch (const CorkCancelled &) {
        throw_cancelled();
    }
}

MeshTuple booleanIntersection(const EigenVecX3d &vertsA,
                              const EigenVecX3i &trisA,
                              const EigenVecX3d &vertsB,
                              const EigenVecX3i &trisB) {

    cork_reset_cancel();
    try {
        py::gil_scoped_release release;
        CorkMesh meshA, meshB;

        eigenToCorkMesh(vertsA, trisA, &meshA);
        eigenToCorkMesh(vertsB, trisB, &meshB);

        MeshTuple meshOut;

        meshA.boolIsct(meshB);

        corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

        return meshOut;
    } catch (const CorkCancelled &) {
        throw_cancelled();
    }
}


MeshTuple booleanXor(const EigenVecX3d &vertsA,
                     const EigenVecX3i &trisA,
                     const EigenVecX3d &vertsB,
                     const EigenVecX3i &trisB) {

    cork_reset_cancel();
    try {
        py::gil_scoped_release release;
        CorkMesh meshA, meshB;

        eigenToCorkMesh(vertsA, trisA, &meshA);
        eigenToCorkMesh(vertsB, trisB, &meshB);

        MeshTuple meshOut;

        meshA.boolXor(meshB);

        corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

        return meshOut;
    } catch (const CorkCancelled &) {
        throw_cancelled();
    }
}


MeshTuple resolveIntersection(const EigenVecX3d &vertsA,
                              const EigenVecX3i &trisA) {

    cork_reset_cancel();
    try {
        py::gil_scoped_release release;
        CorkMesh meshA;

        eigenToCorkMesh(vertsA, trisA, &meshA);

        MeshTuple meshOut;

        meshA.resolveIntersections();

        corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

        return meshOut;
    } catch (const CorkCancelled &) {
        throw_cancelled();
    }
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

    m.def("request_cancel", &cork_request_cancel,
          "Request cooperative cancel of an in-flight Cork operation")
     .def("reset_cancel", &cork_reset_cancel,
          "Clear the cooperative cancel flag before starting a new Cork op")
     .def("is_cancelled", &cork_is_cancelled,
          "True if cancel was requested for the current Cork op")
     .def("isSolid", &pycork::isSolid, "Determines if the mesh is manifold",
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
     .def("xor", &pycork::booleanXor, "Computes boolean xor between two meshes",
                   py::arg("vertsA"), py::arg("trisA"),
                   py::arg("vertsB"), py::arg("trisB"))
     .def("resolveIntersection", &pycork::resolveIntersection, "Computes the intersection between two meshes",
                                  py::arg("vertsA"), py::arg("trisA"));


#ifdef PROJECT_VERSION
    m.attr("__version__") = "PROJECT_VERSION";
#else
    m.attr("__version__") = "dev";
#endif

}
