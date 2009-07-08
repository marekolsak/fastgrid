CMAKE_MINIMUM_REQUIRED(VERSION 2.6.3)

PROJECT(autogrid)
SET(CMAKE_BUILD_TYPE Release)
SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};--ptxas-options=-v)

SET(AUTOGRID_SOURCES
    ./ann/src/perf.cpp
    ./ann/src/kd_dump.cpp
    ./ann/src/bd_tree.cpp
    ./ann/src/ANN.cpp
    ./ann/src/kd_search.cpp
    ./ann/src/brute.cpp
    ./ann/src/kd_tree.cpp
    ./ann/src/bd_search.cpp
    ./ann/src/bd_fix_rad_search.cpp
    ./ann/src/kd_util.cpp
    ./ann/src/kd_pr_search.cpp
    ./ann/src/kd_fix_rad_search.cpp
    ./ann/src/bd_pr_search.cpp
    ./ann/src/kd_split.cpp

    ./math/Geometric/Sphere3.cpp
    ./math/Geometric/Matrix3.cpp
    ./math/Geometric/Ray2.cpp
    ./math/Geometric/AxisAlignedBox3.cpp
    ./math/Geometric/Vertex3.cpp
    ./math/Geometric/Matrix4.cpp
    ./math/Geometric/Vec4.cpp
    ./math/Geometric/Intersection2.cpp
    ./math/Geometric/Frustum3.cpp
    ./math/Geometric/Plane3.cpp
    ./math/Geometric/Line2.cpp
    ./math/Geometric/Ray3.cpp
    ./math/Geometric/Quaternion.cpp
    ./math/Geometric/Intersection3.cpp

    ./times.cpp
    ./desolvexpfunc.cpp
    ./inputdataloader.cpp
    ./parameterlibrary.cpp
    ./pairwiseinteractionenergies.cpp
    ./main.cpp
    ./bondvectors.cpp
    ./logfile.cpp
    ./calculategridmaps.cpp
    ./programparameters.cpp
    ./utils.cpp
    ./gridmap.cpp

    ./electrostatics/cpu.cpp
)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/cuda" ${CMAKE_MODULE_PATH})

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    ADD_DEFINITIONS(-DAG_OPENMP)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF(OPENMP_FOUND)

FIND_PACKAGE(CUDA)
#IF(CUDA_FOUND)
    ADD_DEFINITIONS(-DAG_CUDA)
    INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})

    SET(AUTOGRID_SOURCES
        ${AUTOGRID_SOURCES}

        ./electrostatics/cuda.cpp
        ./electrostatics/cuda_internal.cu
    )

    CUDA_ADD_EXECUTABLE(autogrid4 ${AUTOGRID_SOURCES})

    CUDA_BUILD_CLEAN_TARGET()
#ELSE(CUDA_FOUND)
#    ADD_EXECUTABLE(autogrid4 ${AUTOGRID_SOURCES})
#ENDIF(CUDA_FOUND)

