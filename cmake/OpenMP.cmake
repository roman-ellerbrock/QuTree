
if (NOT openmp)
    set(openmp OFF)
    message("Setting openMP off (${openmp})")
else ()
    message("openmp=${openmp}")
endif()

if (openmp)
    find_package(OpenMP REQUIRED)
    target_link_libraries(QuTree PRIVATE OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran)
    set(QuTree_RELEASE_FLAGS "-openmp ${QuTree_RELEASE_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp")
endif()

