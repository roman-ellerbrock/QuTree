

execute_process(COMMAND git submodule update --init -- external/PolymorphicMemory)
add_subdirectory("external/PolymorphicMemory")
set(PolymorphicMemory_DIR "${PROJECT_SOURCE_DIR}/external/PolymorphicMemory/src")

message("Added PolymorphicMemory package: ${PolymorphicMemory_DIR}")

