//
// Created by Roman Ellerbrock on 2/1/20.
//

#ifndef SPECTRALDECOMPOSITIONTREE_H
#define SPECTRALDECOMPOSITIONTREE_H
#include "NodeAttribute.h"
#include "Core/Matrix.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "Core/Vector.h"

template <typename T>
class SpectralDecompositionTree : public NodeAttribute<SpectralDecomposition<T>> {
public:
	using NodeAttribute<SpectralDecomposition<T>>::attributes_;

	SpectralDecompositionTree() = default;

	explicit SpectralDecompositionTree(const Tree& tree);

	SpectralDecompositionTree(const MatrixTree<T>& H, const Tree& tree);

	~SpectralDecompositionTree() = default;

	void Initialize(const Tree& tree);

	void Calculate(const MatrixTree<T>& H, const Tree& tree);

	MatrixTree<T> Invert(const Tree& tree, double eps = 1e-7);

	/// I/O
	void print(const Tree& tree) const;
};

template<typename T>
void CanonicalTransformation(TensorTree<T>& Psi, const Tree& tree, bool orthogonal = false);

template <typename T>
SpectralDecompositionTree<T> sqrt(SpectralDecompositionTree<T> X);

template <typename T>
MatrixTree<T> sqrt(MatrixTree<T> X, const Tree& tree);

template <typename T>
SpectralDecompositionTree<T> inverse(SpectralDecompositionTree<T> X, double eps = 1e-7);

template <typename T>
MatrixTree<T> inverse(MatrixTree<T> X, const Tree& tree, double eps = 1e-7);

template <typename T>
MatrixTree<T> to_matrixtree(const SpectralDecompositionTree<T>& X, const Tree& tree);

typedef SpectralDecompositionTree<complex<double>> SpectralDecompositionTreecd;
typedef SpectralDecompositionTree<double> SpectralDecompositionTreed;


#endif //SPECTRALDECOMPOSITIONTREE_H
