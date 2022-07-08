
# Option to use internal eigen

option(USE_INTERNAL_EIGEN "Use internal Eigen module" OFF)

if (NOT USE_INTERNAL_EIGEN)
    find_package(Eigen3 NO_MODULE)
endif()

if (USE_INTERNAL_EIGEN OR NOT TARGET Eigen3::Eigen)
    message ( STATUS "Eigen not found, switching to internal submodule" )
    set(USE_INTERNAL_EIGEN ON)

    execute_process(COMMAND git submodule update --init -- external/eigen
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    # Set env and install headers as subdir of QuTree to avoid clashes
    set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen
        CACHE PATH "Eigen include directory")
    install(DIRECTORY ${EIGEN3_INCLUDE_DIR}/Eigen
        DESTINATION ${QuTree_INSTALL_INC_DIR})

    # Convenience target for exports
    add_library(Eigen INTERFACE)
    add_library(Eigen3::Eigen ALIAS Eigen)
    target_include_directories(Eigen INTERFACE
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:${QuTree_INSTALL_INC_DIR}>)
    install(TARGETS Eigen
        EXPORT QuTree-export
        DESTINATION ${QuTree_INSTALL_INC_DIR})
endif()
