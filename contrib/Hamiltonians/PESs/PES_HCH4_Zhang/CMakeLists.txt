
add_library(PES_HCH4_Zhang PES_HCH4_Zhang.f PES_HCH4_Zhang.cpp)
target_include_directories(PES_HCH4_Zhang PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(PES_HCH4_Zhang MCTDHlib)

