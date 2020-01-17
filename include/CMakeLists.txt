
set(QDlib_INCLUDE_FILES
        BS_integrator.h
        DoubleParticleOperator.h
        FFT.h
        FFTCooleyTukey.h
        JacobiRotationFramework.h
        Lanzcos.h
        Matrix.h
        Matrix_Implementation.h
        MultiIndex.h
        QDlib.h
        QMConstants.h
        SimultaneousDiagonalization.h
        SingleParticleOperator.h
        Tensor.h
        TensorDim.h
        TensorDim_Extension.h
        Tensor_Extension.h
        Tensor_Extension_Implementation.h
        Tensor_Implementation.h
        Tree.h
        TreeNode.h
        Vector.h
        Vector_Implementation.h
        long_integer.h
        stdafx.h
        )

install(FILES ${QDlib_INCLUDE_FILES} DESTINATION include/QDlib)
