//
// Created by Roman Ellerbrock on 2/1/20.
//

#ifndef SPECTRALDECOMPOSITIONTREE_H
#define SPECTRALDECOMPOSITIONTREE_H
#include "NodeAttribute.h"
#include "Core/Matrix.h"
#include "Tree/MatrixTreeFunctions.h"
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

typedef SpectralDecompositionTree<complex<double>> SpectralDecompositionTreecd;
typedef SpectralDecompositionTree<double> SpectralDecompositionTreed;


#endif //SPECTRALDECOMPOSITIONTREE_H
