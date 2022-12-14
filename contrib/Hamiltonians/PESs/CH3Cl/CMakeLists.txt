
set (FFLAGS "-mcmodel=large")
add_library(PES_CH3Cl PES_CH3Cl.cpp coordToEn.f ctrans_MeX.f
	  rs.f model_MeX.f vwz_MeX.f)
target_include_directories(PES_CH3Cl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(PES_CH3Cl MCTDHlib)

