
set(QDlib_SOURCE_FILES
        TensorDim.cpp
        TensorDim_Extension.cpp
        FFT.cpp
        JacobiRotationFramework.cpp
        SimultaneousDiagonalization.cpp
        CWrapper.f
        long_integer.cpp
        Vector_Instantiations.cpp Matrix_Instantiations.cpp Tensor_Instantiations.cpp Tensor_Extension_Instantiations.cpp FactorMatrix_Instantiations.cpp)


add_library(QDlib ${QDlib_SOURCE_FILES})
install(TARGETS QDlib DESTINATION lib)

