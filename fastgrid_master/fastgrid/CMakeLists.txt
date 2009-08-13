CMAKE_MINIMUM_REQUIRED(VERSION 2.6.3)

PROJECT(fastgrid)
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
    ./math/Geometric/AxisAlignedBox3.cpp
    ./math/Geometric/Vec4.cpp
    ./math/Geometric/Intersection3.cpp

    ./openthreads/pthread_src/PThread.c++
    ./openthreads/pthread_src/PThreadBarrier.c++
    ./openthreads/pthread_src/PThreadCondition.c++
    ./openthreads/pthread_src/PThreadMutex.c++

    ./times.cpp
    ./DesolvExpFunc.cpp
    ./InputDataLoader.cpp
    ./ParameterLibrary.cpp
    ./PairwiseInteractionEnergies.cpp
    ./main.cpp
    ./BondVectors.cpp
    ./LogFile.cpp
    ./CalculateGridmaps.cpp
    ./ProgramParameters.cpp
    ./Utils.cpp
    ./GridMap.cpp

    ./electrostatics/Cpu.cpp
)

ADD_DEFINITIONS(-DLinux -DFASTGRID)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread")

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

        ./electrostatics/Cuda.cpp
        ./electrostatics/CudaConstantMemory.cpp
        ./electrostatics/CudaEvents.cpp
        ./electrostatics/CudaFloatTexture1D.cpp
        ./electrostatics/CudaGridMap.cpp
        ./electrostatics/cuda_internal/Interface.cpp
        ./electrostatics/cuda_internal/StandardKernels.cu
        ./electrostatics/cuda_internal/DDDConstMemKernels.cu
    )

    CUDA_ADD_EXECUTABLE(fastgrid4 ${AUTOGRID_SOURCES})

    CUDA_BUILD_CLEAN_TARGET()
#ELSE(CUDA_FOUND)
#    ADD_EXECUTABLE(fastgrid4 ${AUTOGRID_SOURCES})
#ENDIF(CUDA_FOUND)
