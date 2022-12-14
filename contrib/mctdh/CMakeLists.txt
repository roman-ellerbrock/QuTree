
find_package(QuTree REQUIRED)
find_package(yaml-cpp REQUIRED)

include_directories(${YAML_CPP_INCLUDE_DIR})
add_library(Yaml INTERFACE IMPORTED)

add_executable(mctdh
        main.cpp
        mctdh.cpp
        mctdh.h
        mctdh_state.h
        )

target_compile_features(mctdh PRIVATE cxx_std_14)
target_compile_options(mctdh PUBLIC "$<$<CONFIG:DEBUG>:${mctdh_DEBUG_FLAGS}>")
target_compile_options(mctdh PRIVATE "$<$<CONFIG:RELEASE>:${mctdh_RELEASE_FLAGS}>")
target_include_directories(mctdh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(mctdh PUBLIC "-ffast-math" "-O3")
target_link_libraries(mctdh QuTree::QuTree yaml-cpp Hamiltonians)

