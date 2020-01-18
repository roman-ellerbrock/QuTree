
set(QDlib_SOURCE_FILES
        TensorDim.cpp
        TensorDim_Extension.cpp
        TemplateInstantiations_1.cpp
        TemplateInstantiations_2.cpp
        FFT.cpp
        JacobiRotationFramework.cpp
        SimultaneousDiagonalization.cpp
        stdafx.cpp
        CWrapper.f
        long_integer.cpp
        )


add_library(QDlib ${QDlib_SOURCE_FILES})
install(TARGETS QDlib DESTINATION lib)

