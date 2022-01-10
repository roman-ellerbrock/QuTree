
set(QuTree_benchmarks
        benchmark_tensor.cpp
        benchmark_tensor.h
        Benchmark.cpp
        Benchmark.h
        )

add_executable(benchmarks benchmarks.cpp ${QuTree_benchmarks})
target_link_libraries(benchmarks QuTree /usr/local/opt/openblas/lib/libblas.dylib)
target_compile_features(benchmarks PRIVATE cxx_std_14)
target_compile_options(benchmarks PRIVATE "-O3")

