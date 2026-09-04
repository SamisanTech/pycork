cmake_minimum_required(VERSION 3.12)

project(pycork VERSION 0.1.3)

# C++14 for OpenMP-friendly lambdas / better codegen
set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cMake")

# Copy the include directory if specified
set(INCLUDE_INSTALL_DIR "${INCLUDE_INSTALL_DIR}/pycork")

# Generate DLL Pragma for creating shared libraries in Windows
include (GenerateExportHeader)

add_definitions("-DPROJECT_VERSION=\"${PROJECT_VERSION}\"" )

# NOTE: CORK_DEBUG removed — it was a major performance footgun.
# CORK_FAST disables ENSURE checks in release hot paths.
add_definitions("-DNOMINMAX -DFMT_HEADER_ONLY -DCORK_FAST -DNDEBUG -DANSI_DECLARATORS -D_USE_MATH_DEFINES -DTRILIBRARY
                 -DNO_TIMER -DTRIANGLE_INTEL64 -DREDUCED -D_CRT_SECURE_NO_WARNINGS -DCDT_ONLY")

option(BUILD_PYTHON "Builds a python extension" ON)
option(CORK_OPENMP "Enable OpenMP for SI detection" ON)
option(CORK_TBB "Fetch and enable oneTBB for SI detection" ON)
option(CORK_ISPC "Compile ISPC AVX2 bbox batch kernel" ON)

include(FetchContent)

if(CORK_TBB)
    # Use a locally vendored oneTBB tree (no network / remote fetch).
    # Expected at external/oneTBB (copy from a prior FetchContent tbb-src if needed).
    set(_CORK_TBB_SRC "${CMAKE_SOURCE_DIR}/external/oneTBB")
    if(NOT EXISTS "${_CORK_TBB_SRC}/CMakeLists.txt")
        message(FATAL_ERROR "CORK_TBB enabled but external/oneTBB is missing. Place a local oneTBB source tree there.")
    endif()
    FetchContent_Declare(
        tbb
        SOURCE_DIR ${_CORK_TBB_SRC}
    )
    set(TBB_TEST OFF CACHE BOOL "" FORCE)
    set(TBB_STRICT OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(tbb)
    message(STATUS "oneTBB enabled for cork")
    add_compile_definitions(CORK_USE_TBB=1)
endif()

# ISPC toolchain
set(ISPC_EXECUTABLE "${CMAKE_SOURCE_DIR}/external/ispc/ispc-v1.25.3-windows/bin/ispc.exe")
if(CORK_ISPC AND EXISTS "${ISPC_EXECUTABLE}")
    message(STATUS "ISPC found: ${ISPC_EXECUTABLE}")
    set(ISPC_SRC ${CMAKE_SOURCE_DIR}/cork/accel/bbox_batch.ispc)
    set(ISPC_OBJ ${CMAKE_BINARY_DIR}/bbox_batch.ispc.obj)
    set(ISPC_HDR ${CMAKE_BINARY_DIR}/bbox_batch_ispc.h)
    add_custom_command(
        OUTPUT ${ISPC_OBJ} ${ISPC_HDR}
        COMMAND "${ISPC_EXECUTABLE}" -O3 --target=avx2-i32x8 --arch=x86-64
                -h ${ISPC_HDR} -o ${ISPC_OBJ} ${ISPC_SRC}
        DEPENDS ${ISPC_SRC}
        COMMENT "ISPC compiling bbox_batch.ispc (avx2-i32x8)"
        VERBATIM
    )
    add_custom_target(cork_ispc_gen DEPENDS ${ISPC_OBJ} ${ISPC_HDR})
    add_compile_definitions(CORK_USE_ISPC=1)
    include_directories(${CMAKE_BINARY_DIR})
    set(CORK_ISPC_OBJ ${ISPC_OBJ})
else()
    message(STATUS "ISPC skipped")
    set(CORK_ISPC OFF)
    set(CORK_ISPC_OBJ "")
endif()

if(WIN32)
    # Remove Security definitions for the library
    # Remove run time checks for windows

    add_subdirectory(mpir)

    IF(MSVC)

        set(COMMON_LANGUAGE_RUNTIME "")
        set(EIGEN3_INCLUDE_DIR external/eigen)

        # Aggressive Release codegen for cork
        add_compile_options(
            "/GS-"
            "/INCREMENTAL:NO"
            "$<$<CONFIG:Release>:/O2>"
            "$<$<CONFIG:Release>:/Ob2>"
            "$<$<CONFIG:Release>:/Oi>"
            "$<$<CONFIG:Release>:/Ot>"
            "$<$<CONFIG:Release>:/GL>"
            "$<$<CONFIG:Release>:/arch:AVX2>"
            "$<$<CONFIG:Release>:/fp:precise>"
            "$<$<CONFIG:Release>:/Gy>"
        )
        add_link_options("$<$<CONFIG:Release>:/LTCG>")

        if(CORK_OPENMP AND NOT CORK_TBB)
            find_package(OpenMP)
            if(OpenMP_CXX_FOUND)
                message(STATUS "OpenMP enabled for cork")
                add_compile_definitions(CORK_USE_OPENMP=1)
            else()
                message(WARNING "OpenMP requested but not found")
            endif()
        endif()

        set(GMP_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/mpir)

        if(BUILD_PYTHON)
            set(GMP_LIBRARIES mpir_static)
        else(BUILD_PYTHON)
            set(GMP_LIBRARIES mpir)
        endif(BUILD_PYTHON)

    endif(MSVC)

    # Copy the shared libraries and libs to the main directory
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

else(WIN32)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

    set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

    # Find Required packages
    find_package(Eigen3 3.3 REQUIRED)
    find_package(GMP REQUIRED)

endif(WIN32)


# Use the replacement of Boost::filesystem from a git submodule provided by WJakob
# in order to reduce compile time dependencies
# filesystem (vendored under external/filesystem)

set(CPP_FILESYSTEM
    external/filesystem
)

## A replacement library for std::Format
set(CPP_FMT
    external/fmt/include
)

#add_definitions("-DFMT_HEADER_ONLY=1 -DNO_TIMER -DANSI_DECLARATORS -DREDUCED -DTRIANGLE_INTEL64 -CDT_ONLY -DCPU86 -DTRILIBRARY -DANSI_DECLARATORS -D_USE_MATH_DEFINES - -")

#add_compile_definitions(

link_directories(
    ${GMP_LIBRARY_DIR}
)

if(BUILD_PYTHON)
    message(STATUS "Building PyCork Python Module")

    add_subdirectory(external/pybind11)

    set(PYBIND11_CPP_STANDARD /std:c++14)


endif(BUILD_PYTHON)

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_BINARY_DIR}
    ${GMP_INCLUDE_DIRS}

    ${EIGEN3_INCLUDE_DIR}
    ${CPP_FILESYSTEM}
    ${CPP_FMT}
)

link_directories(
    ${GMP_LIBRARIES}
)

set(APP_H_SRCS
    cork/cork.h
    cork/accel/aabvh.h
    cork/file_formats/files.h
    cork/isct/absext4.h
    cork/isct/empty3d.h
    cork/isct/ext4.h
    cork/isct/fixext4.h
    cork/isct/fixint.h
    cork/isct/gmpext4.h
    cork/isct/quantization.h
    cork/isct/triangle.h
    cork/isct/unsafeRayTriIsct.h
    cork/math/bbox.h
    cork/math/ray.h
    cork/math/vec.h
    cork/mesh/mesh.h
    cork/rawmesh/rawMesh.h
    cork/util/iterPool.h
    cork/util/memPool.h
    cork/util/prelude.h
    cork/util/shortVec.h
    cork/util/unionFind.h
    cork/mesh/mesh.remesh.tpp
    cork/mesh/mesh.isct.tpp
    cork/mesh/mesh.bool.tpp
    cork/mesh/mesh.cpp
)

set(APP_CPP_SRCS
    cork/cork.cpp
    cork/file_formats/files.cpp
    cork/file_formats/ifs.cpp
    cork/file_formats/off.cpp
    cork/isct/empty3d.cpp
    cork/isct/quantization.cpp
    cork/isct/triangle.cpp
    cork/util/log.cpp
    cork/util/timer.cpp

)

set(APP_SRCS
    ${APP_H_SRCS}
    ${APP_CPP_SRCS}
)

SOURCE_GROUP("App" FILES
    ${APP_SRCS}
)

set(LIBCORK_SRCS
    ${APP_SRCS}
)

set(PYCORK_LIBS
    ${GMP_LIBRARIES}
    )

if(BUILD_PYTHON)
    # Add the library

    message(STATUS "Building PyCork Python-  Static Library")

    add_definitions("-DCORK_BUILT_AS_STATIC")
    add_library(cork_static STATIC ${LIBCORK_SRCS})

    target_link_libraries(cork_static PRIVATE ${PYCORK_LIBS})
    if(TARGET tbb)
        target_link_libraries(cork_static PRIVATE TBB::tbb)
    endif()
    if(OpenMP_CXX_FOUND)
        target_link_libraries(cork_static PRIVATE OpenMP::OpenMP_CXX)
    endif()
    if(CORK_ISPC AND CORK_ISPC_OBJ)
        set_source_files_properties(${CORK_ISPC_OBJ} PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)
        target_sources(cork_static PRIVATE ${CORK_ISPC_OBJ})
        add_dependencies(cork_static cork_ispc_gen)
    endif()
    GENERATE_EXPORT_HEADER(cork_static
                 BASE_NAME CORK
                 EXPORT_MACRO_NAME CORK_EXPORT
                 EXPORT_FILE_NAME CORK_Export.h
                 STATIC_DEFINE CORK_BUILT_AS_STATIC)


else(BUILD_PYTHON)
    message(STATUS "Building Cork Module Standalone - Dynamic Library")

    add_library(cork SHARED ${LIBCORK_SRCS})
    target_link_libraries(cork ${PYCORK_LIBS})

    GENERATE_EXPORT_HEADER(cork
                 BASE_NAME CORK
                 EXPORT_MACRO_NAME CORK_EXPORT
                 EXPORT_FILE_NAME CORK_Export.h
                 STATIC_DEFINE CORK_BUILT_AS_STATIC)

endif(BUILD_PYTHON)

set(App_SRCS
    cork/main.cpp
)

SOURCE_GROUP("App" FILES ${App_SRCS})

if(BUILD_PYTHON)

    set(PYCORK_SRCS

        python/pycork/module.cpp
    )

    SOURCE_GROUP("Python" FILES ${PYCORK_SRCS})

    pybind11_add_module(pycork ${PYCORK_SRCS})

    #add_library(example MODULE main.cpp)

    target_link_libraries(pycork PRIVATE pybind11::module cork_static ${PYCORK_LIBS})
    if(TARGET tbb)
        target_link_libraries(pycork PRIVATE TBB::tbb)
    endif()
    set_target_properties(pycork PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                          SUFFIX "${PYTHON_MODULE_EXTENSION}")


     add_executable(main ${App_SRCS})
     target_link_libraries(main cork_static ${PYCORK_LIBS})

    install(TARGETS pycork DESTINATION lib/pycork)

    #install(TARGETS mpir_static DESTINATION )
else(BUILD_PYTHON)


    add_executable(main ${App_SRCS})
    target_link_libraries(main cork ${PYCORK_LIBS})

    #install(TARGETS mpir DESTINATION )

endif(BUILD_PYTHON)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CORK_Export.h
    ${APP_H_SRCS}
    DESTINATION include/pycork
)

