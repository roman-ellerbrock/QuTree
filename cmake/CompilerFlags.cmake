
target_compile_features(QuTree PUBLIC cxx_std_17)

if (APPLE)
    message("On OSX")
    set(QuTree_DEBUG_FLAGS "-mmacosx-version-min=10.14 ${QuTree_DEBUG_FLAGS}")
    set(QuTree_RELEASE_FLAGS "-mmacosx-version-min=10.14 ${QuTree_RELEASE_FLAGS}")
    set(CMAKE_MACOSX_RPATH 1)
elseif (UNIX)
    message("On UNIX")
else ()
    message("Other OS")
endif ()

# Get command-line ready options for target_compile: https://stackoverflow.com/a/27651464/3052876
separate_arguments(QuTree_DEBUG_FLAGS UNIX_COMMAND "${QuTree_DEBUG_FLAGS}")
separate_arguments(QuTree_RELEASE_FLAGS UNIX_COMMAND "${QuTree_RELEASE_FLAGS}")

# Set debug/release flags: https://stackoverflow.com/a/23995391/3052876
target_compile_options(QuTree PUBLIC "$<$<CONFIG:DEBUG>:${QuTree_DEBUG_FLAGS}>")
target_compile_options(QuTree PRIVATE "$<$<CONFIG:RELEASE>:${QuTree_RELEASE_FLAGS}>")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

