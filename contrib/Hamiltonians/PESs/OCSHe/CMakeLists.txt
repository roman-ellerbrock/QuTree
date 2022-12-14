
add_library(OCSHe_PES ocsheextpot-bdheocs.f OCSHePotential.cpp ocsheplm.f ocshesumleg.f ocshepot.f)
target_include_directories(OCSHe_PES PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(OCSHe_PES MCTDHlib)

