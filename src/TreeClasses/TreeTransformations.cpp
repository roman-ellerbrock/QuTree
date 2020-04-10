//
// Created by Roman Ellerbrock on 4/9/20.
//

#include "TreeClasses/TreeTransformation_Implementation.h"

typedef complex<double> cd;
typedef double d;

namespace TreeFunctions {

template void Transform(TensorTree<cd>& Chi, const TensorTree<cd>& Psi, const MatrixTree<cd>& M, const MatrixTree<cd>& M_inv, const Tree& tree);
template void Transform(TensorTree<d>& Chi, const TensorTree<d>& Psi, const MatrixTree<d>& M, const MatrixTree<d>& M_inv, const Tree& tree);
template TensorTree<cd> ContractionNormalization(TensorTree<cd> Psi, const Tree& tree, bool orthogonal);
template TensorTree<d> ContractionNormalization(TensorTree<d> Psi, const Tree& tree, bool orthogonal);
template TensorTree<cd> DotProductNormalization(TensorTree<cd> Psi, const Tree& tree);
template TensorTree<d> DotProductNormalization(TensorTree<d> Psi, const Tree& tree);

}
