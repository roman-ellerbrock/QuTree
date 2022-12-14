
# Link
add_library(PESs PESs.h
        CDVRModelV/CDVRModelV.h
        CDVRModelV/CDVRModelV.cpp
        CH3_1/CH3Potential.cpp
        CH3_1/CH3Potential.h
        CH3_1/CH3Potential.f
        NOCl/NOClPotential.h
        NOCl/NOClPotential.cpp
        liu_ch4cl/pipNN/potential.f
        liu_ch4cl/pipNN/ch4cl-pipnn.f
        liu_ch4cl/liuch4cl.cpp
        liu_ch4cl/liuch4cl.h)

#add_subdirectory("CDVRModelV")
#target_link_libraries(PESs LINK_PUBLIC
#        CDVRModelV
#        )

target_include_directories(PESs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src)
target_compile_options(PESs PUBLIC "-Ofast")
target_link_libraries(PESs LINK_PUBLIC QuTree::QuTree)

######################################################
# List the active PES's subdirectories here
# Note: Remember to link them to the libraries
######################################################
# Methyl PES
#add_subdirectory("NOCl")
#target_link_libraries(PESs LINK_PUBLIC NOCl)

#add_subdirectory("CH3_1")
#add_subdirectory("CH4_schwenke")
#add_subdirectory("CH4_SM")
#add_subdirectory("CH5+_Bowman")
#add_subdirectory("PES_HCH4_Zhang")
#add_subdirectory("hehe")
#add_subdirectory("OCSHe")
#add_subdirectory("CDVRModelV")
#add_subdirectory("CDVRModelV2")
#add_subdirectory("CDVRModelV3")
#add_subdirectory("CH3Cl")
#add_subdirectory("LinearVibronicCoupling")
#add_subdirectory("Model1")
#
#target_link_libraries(PESs LINK_PUBLIC PESCH3_1
#        PES_CH4
#        PES_CH3Cl
#        PES_CH5P
#        PES_CH4SM
#        PES_HCH4_Zhang
#        NOCl_PES
#        hehe_PES
#        OCSHe_PES
#        CDVRModelV
#        CDVRModelV2
#        CDVRModelV3
#        LinearVibronicCoupling
#        Model1
#        )

