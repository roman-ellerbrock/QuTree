//
// Created by Roman Ellerbrock on 2/1/20.
//

#ifndef SPECTRALDECOMPOSITIONTREE_H
#define SPECTRALDECOMPOSITIONTREE_H
#include "TreeStructuredObject.h"
#include "Core/Matrix.h"
#include "Tree/MatrixTreeFunctions.h"
#include "Core/Vector.h"

template <typename T>
class SpectralDecompositionTree : public TreeStructuredObject<SpectralDecomposition<T>> {
public:
	using TreeStructuredObject<SpectralDecomposition<T>>::attributes_;

	SpectralDecompositionTree() = default;

	explicit SpectralDecompositionTree(const TTBasis& tree);

	SpectralDecompositionTree(const MatrixTree<T>& H, const TTBasis& tree);

	~SpectralDecompositionTree() = default;

	void Initialize(const TTBasis& tree);

	void Calculate(const MatrixTree<T>& H, const TTBasis& tree);

	MatrixTree<T> Invert(const TTBasis& tree, double eps = 1e-7);

	/// I/O
	void print(const TTBasis& tree) const;
};

typedef SpectralDecompositionTree<complex<double>> SpectralDecompositionTreecd;
typedef SpectralDecompositionTree<double> SpectralDecompositionTreed;


#endif //SPECTRALDECOMPOSITIONTREE_H
