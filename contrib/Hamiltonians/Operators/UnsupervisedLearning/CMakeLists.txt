
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(UnsupervisedLearning MNISTdata.cpp MNISTdata.h)

target_include_directories (UnsupervisedLearning PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (UnsupervisedLearning MCTDHlib)

