# Project
cmake_minimum_required(VERSION 3.12)
project(
  QuTree
  VERSION 2.0
  LANGUAGES CXX)

option(BUILD_TESTS "Build testing framework" ON)

# If you set any CMAKE_ variables, that can go here.
# (But usually don't do this, except maybe for C++ standard)

# Find packages go here.
find_package(Torch REQUIRED)
find_package(yaml-cpp REQUIRED)
# Add this, because yaml-cpp.cmake sets empty include & library dirs (yeah, great..)
get_property(YAML_CPP_INCLUDE_DIR  TARGET yaml-cpp  PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(YAML_CPP_LIBRARIES yaml-cpp IMPORTED_LOCATION_RELEASE)
find_package(nlohmann_json REQUIRED)

# General flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# This is a "default" library, and will match the *** variable setting.
# Other common choices are STATIC, SHARED, and MODULE
# Including header files here helps IDEs but is not required.
# Output libname matches target name, with the usual extensions on your system
add_subdirectory(src)

# Link each target with other targets or add options, etc.

# The project's main target
add_executable(qutree ${PROJECT_SOURCE_DIR}/src/main.cpp)

set_target_properties(qutree PROPERTIES CXX_STANDARD 17)

# Make sure you link your targets with this command. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.
target_link_libraries(qutree PRIVATE 
    qutreelib
    ${YAML_CPP_LIBRARIES}
    ${TORCH_LIBRARIES}
    nlohmann_json::nlohmann_json
)

if (BUILD_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
endif()


