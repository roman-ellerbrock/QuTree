
#####################################################################
# Target Installation setup
#####################################################################

set(QuTree_INSTALL_INC_DIR ${CMAKE_INSTALL_PREFIX}/include/QuTree)
set(QuTree_INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
set(QuTree_INSTALL_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/lib/cmake/QuTree)

# Set up headers in include/ for install & build and src/ for just building
target_include_directories(QuTree
    PUBLIC
        $<INSTALL_INTERFACE:${QuTree_INSTALL_INC_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        ${BLAS_INCLUDE}
	${CUDA_INCLUDE_DIRS}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/deprecated
    )

# Install rules
install(TARGETS QuTree
    EXPORT QuTree-export
    LIBRARY DESTINATION ${QuTree_INSTALL_LIB_DIR}
    ARCHIVE DESTINATION ${QuTree_INSTALL_LIB_DIR}
    )

install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/
    DESTINATION ${QuTree_INSTALL_INC_DIR}
    FILES_MATCHING PATTERN "*.h*")

install(EXPORT QuTree-export
    FILE QuTreeTargets.cmake
    NAMESPACE QuTree::
    DESTINATION ${QuTree_INSTALL_CMAKE_DIR})

# Exports for external find_package
include(CMakePackageConfigHelpers)
configure_file(QuTreeConfig.cmake.in QuTreeConfig.cmake @ONLY)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/QuTreeConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    )

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/QuTreeConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/QuTreeConfigVersion.cmake
    DESTINATION ${QuTree_INSTALL_CMAKE_DIR})

