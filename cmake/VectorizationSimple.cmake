
set(FLAGS_VECTORIZE "-fopenmp-simd -march=native")
set(QuTree_RELEASE_FLAGS "${FLAGS_VECTORIZE} -O3 -ffast-math")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp -Ofast -ffast-math")
