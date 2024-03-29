
# Compile the Operators
add_subdirectory("Operators")

# Compile the Potential Energy Surfaces
add_subdirectory("PESs")

find_package(yaml-cpp REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})

# Link the libraries "Operators" and "PESs" into a single library "Hamiltonians"
add_library(Hamiltonians SHARED
        Hamiltonians.cpp
        Hamiltonians.h
        )

target_compile_features(Hamiltonians PRIVATE cxx_std_14)
target_link_libraries(Hamiltonians LINK_PUBLIC Operators PESs)
target_link_libraries(Hamiltonians LINK_PUBLIC QuTree::QuTree)
target_link_libraries(Hamiltonians LINK_PUBLIC yaml-cpp)
target_compile_options(Hamiltonians PUBLIC "-ffast-math" "-O3")
target_include_directories(Hamiltonians PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS Hamiltonians
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
)
