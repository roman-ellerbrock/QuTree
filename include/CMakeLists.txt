# Easily regenerated with find include -name '*.h' | sort
set(QuTree_INCLUDE_FILES
        QuTree.h

        include/Applications/CMFIntegrator.h
        include/Applications/Eigenstates.h
        include/Applications/SCF.h

        include/Core/Matrix.h
        include/Core/Matrix_Implementation.h
        include/Core/Tensor.h
        include/Core/TensorMath.h
        include/Core/TensorShape.h
        include/Core/Tensor_Extension.h
        include/Core/Tensor_Extension_Implementation.h
        include/Core/Tensor_Implementation.h
        include/Core/Vector.h
        include/Core/Vector_Implementation.h
        include/Core/stdafx.h
        include/QuTree.h

        include/CDVR/CDVR.h
        include/CDVR/cdvr_functions.h
        include/CDVR/DeltaVTree.h
        include/CDVR/MatrixTensorTree.h
        include/CDVR/MatrixTensorTreeFunctions.h
        include/CDVR/SymXMatrixTrees.h
        include/CDVR/TDDVR.h
        include/CDVR/TreeGrids.h
        include/CDVR/XMatrixTrees.h

        include/TreeClasses/EdgeAttribute.h
        include/TreeClasses/IntegratorVariables.h
        include/TreeClasses/LayerInterface.h
        include/TreeClasses/MatrixTensorTree.h
        include/TreeClasses/MatrixTensorTreeFunctions.h
        include/TreeClasses/MatrixTree.h
        include/TreeClasses/MatrixTreeFunctions.h
        include/TreeClasses/MatrixTreeFunctions_Implementation.h
        include/TreeClasses/MatrixTreeTransformations.h
        include/TreeClasses/MatrixTreeTransformations_Implementation.h
        include/TreeClasses/NodeAttribute.h
        include/TreeClasses/SparseMatrixTree.h
        include/TreeClasses/SparseMatrixTreeFunctions.h
        include/TreeClasses/SparseMatrixTreeFunctions_Implementation.h
        include/TreeClasses/SparseNodeAttribute.h
        include/TreeClasses/SparseTensorTree.h
        include/TreeClasses/SparseTensorTree_Implementation.h
        include/TreeClasses/SparseTree.h
        include/TreeClasses/SOPMatrixTrees.h
        include/TreeClasses/SpectralDecompositionTree.h
        include/TreeClasses/SymTensorTree.h
        include/TreeClasses/TensorTree.h
        include/TreeClasses/TensorTree_Implementation.h
        include/TreeClasses/TreeIO.h

        include/TreeClasses/Discrete/SymmetricSCF.h
        include/TreeClasses/Discrete/BlockTree.h

        # could move Fortran support to separate target
        include/TreeOperators/FortranOperator.h
        include/TreeOperators/FortranSOP.h

        include/TreeOperators/LeafFunction.h
        include/TreeOperators/LeafMatrix.h
        include/TreeOperators/LeafOperator.h
        include/TreeOperators/MultiLeafOperator.h
        include/TreeOperators/Hamiltonian.h
        include/TreeOperators/PotentialOperator.h
        include/TreeOperators/SOPVector.h
        include/TreeOperators/SumOfProductsOperator.h
        include/TreeOperators/SumOfProductsOperator_Implementation.h
        include/TreeOperators/TensorOperators/TensorTreeOperator.h
        include/TreeOperators/TensorOperators/TensorOperatorRepresentation.h
        include/TreeOperators/TensorOperators/TTOrepresentation.h
        include/TreeOperators/TensorOperators/TTOcontraction.h
        include/TreeOperators/TensorOperators/TTOHoleTree.h
        include/TreeOperators/TensorOperators/TTOMatrixTree.h

        include/TreeShape/AbstractNode.h
        include/TreeShape/Edge.h
        include/TreeShape/Leaf.h
        include/TreeShape/LeafTypes/DVRBasis.h
        include/TreeShape/LeafTypes/FFTGrid.h
        include/TreeShape/LeafTypes/HO_Basis.h
        include/TreeShape/LeafTypes/LeafInterface.h
        include/TreeShape/LeafTypes/LegendrePolynomials.h
        include/TreeShape/LeafTypes/SpinGroup.h
        include/TreeShape/LinearizedLeaves.h
        include/TreeShape/Node.h
        include/TreeShape/NodePosition.h
        include/TreeShape/Tree.h
        include/TreeShape/TreeFactory.h

        include/Util/BS_integrator.h
        include/Util/FFT.h
        include/Util/FFTCooleyTukey.h
        include/Util/GradientDescent.h
        include/Util/GradientDescent_Implementation.h
        include/Util/JacobiRotationFramework.h
        include/Util/Lanzcos.h
        include/Util/MultiIndex.h
        include/Util/normal_modes.h
        include/Util/Overlaps.h
        include/Util/QMConstants.h
        include/Util/RandomMatrices.h
        include/Util/RandomProjector.h
        include/Util/RandomProjector_Implementation.h
        include/Util/SimultaneousDiagonalization.h
        include/Util/string_ext.h
        include/Util/Tree.h
        include/Util/TreeNode.h
        include/Util/WeightedSimultaneousDiagonalization.h
        )

