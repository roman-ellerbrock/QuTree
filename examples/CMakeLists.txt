#add_executable(purification purification.cpp)
#target_include_directories(purification 
#	PUBLIC
#	$<BUILD_INTERFACE:${CMAKE_PROJECT_DIR}/src>
#	${BLAS_DIRECTORIES}
#	${PolymorphicMemory_DIR}
#        ${BLAS_INCLUDE}
#	${CUDA_INCLUDE_DIRS}
#	)
#
#    message("PROJECT NAME:         ${CMAKE_PROJECT_DIR}")
#
#target_link_libraries(purification QuTree)
#target_compile_features(purification INTERFACE cxx_std_17)


#add_executable(example_1 example_1.cpp)
#target_link_libraries(example_1 QuTree)
#
#add_executable(example_2 example_2.cpp)
#target_link_libraries(example_2 QuTree::QuTree)
#
#add_executable(example_3 example_3.cpp)
#target_link_libraries(example_3 QuTree)
#
#add_executable(example_4 multilayer_example.cpp)
#target_link_libraries(example_4 QuTree)
#
#add_executable(demo demo.cpp)
#target_link_libraries(demo QuTree)

# run unit tests as post build step
#add_custom_target(TestExamples
#    COMMAND "./example_1"
#    COMMAND "./example_2"
#    COMMAND "./example_3"
#    DEPENDS example_1 example_2 example_3
#    COMMENT "Run all examples")
#target_link_libraries(example_3 QuTree::QuTree)
