# Easily regenerated with find include -name '*.h' | sort
set(QuTree_INCLUDE_FILES
        src/QuTree.h
        src/Check.h
        src/stdafx.h

        # Tensor
        src/Tensor/Tensor
        src/Tensor/Tensor.h
        src/Tensor/cuTensor.h
        src/Tensor/mxpTensor.h
        src/Tensor/TensorAllocator.h
        src/Tensor/TensorBLAS1.h
        src/Tensor/TensorBLAS2.h
        src/Tensor/TensorLapack.h
        src/Tensor/TensorShape.h
        src/Tensor/Tensor_Implementation.h

        # Tree
        src/Tree/PrimitiveBasis/BasisAPI.h
        src/Tree/PrimitiveBasis/BasisParameters.h
        src/Tree/PrimitiveBasis/DVR.h
        src/Tree/PrimitiveBasis/FFTGrid.h
        src/Tree/PrimitiveBasis/LegendrePolynomials.h
        src/Tree/PrimitiveBasis/PrimitiveBasis.h
        src/Tree/PrimitiveBasis/SpinGroup.h
        src/Tree/Leaf.h
        src/Tree/NodePosition.h
        src/Tree/Node.h
        src/Tree/NodeArray.h
        src/Tree/Edge.h
        src/Tree/EdgeArray.h
        src/Tree/Tree.h
        src/Tree/Factory.h

        # TreeNetwork
        src/TensorNetwork/TensorTree.h
        src/TensorNetwork/MatrixTree.h

        # Operators
        src/Operator/LeafMatrix.h
        src/Operator/LeafFunction.h
        src/Operator/LeafOperator.h
        src/Operator/ProductOperator.h

        # Utilities
        src/Util/blaspp_device_extension.h
        )

set(QuTree_SOURCE_FILES
        # Tensor/
        src/Tensor/SVD.cpp
        src/Tensor/TensorShape.cpp
        src/Tensor/TensorBLAS1.cpp
        src/Tensor/TensorBLAS2.cpp
        src/Tensor/TensorLapack.cpp
        src/Tensor/TensorSlow.cpp
        src/Tensor/Tensor_Instantiations.cpp
        src/Tensor/cuTensor_Instantiations.cpp
        src/Tensor/cuTensor.cpp
        src/Tensor/mxpTensor.cpp
        src/Tensor/stdafx.cpp

        # Primitve Basis
        src/Tree/PrimitiveBasis/DVR.cpp
        src/Tree/PrimitiveBasis/HarmonicOscillator.cpp
        src/Tree/PrimitiveBasis/LegendrePolynomials.cpp
        src/Tree/PrimitiveBasis/FFTGrid.cpp
        src/Tree/PrimitiveBasis/SpinGroup.cpp
        src/Tree/PrimitiveBasis/BasisAPI.cpp

        # Tree
        src/Tree/NodePosition.cpp
        src/Tree/Node.cpp
        src/Tree/LeafArray.cpp
        src/Tree/Edge.cpp
        src/Tree/SubTree.cpp
        src/Tree/Tree.cpp
        src/Tree/TreeFactory.cpp

        # TensorNetwork
        src/TensorNetwork/ApplySCF.cpp
        src/TensorNetwork/contractions.cpp
#        src/TensorNetwork/cuTensorTree.cpp
        src/TensorNetwork/TensorTree.cpp
        src/TensorNetwork/TensorTreeFactory.cpp

        # Operators
        src/Operator/LeafFunction.cpp
        src/Operator/LeafMatrix.cpp
        src/Operator/ProductOperator.cpp
        src/Operator/SumOfProductsOperator.cpp

        # Utilities
        src/Util/blaspp_device_extension.cpp
        src/Util/GateOperators.cpp
        src/Util/JacobiRotationFramework.cpp
        src/Util/SimultaneousDiagonalization.cpp
        )

if (CUDA_FEATS)
	set (QuTree_INCLUDE_FILES_CUDA 
	)

	set (QuTree_SOURCE_FILES_CUDA 
	)
endif()
