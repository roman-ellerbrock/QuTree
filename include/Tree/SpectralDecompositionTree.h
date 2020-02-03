//
// Created by Roman Ellerbrock on 2/1/20.
//

#ifndef SPECTRALDECOMPOSITIONTREE_H
#define SPECTRALDECOMPOSITIONTREE_H
#include "TreeStructuredObject.h"
#include "Core/Matrix.h"
#include "HoleMatrixTree.h"
#include "Core/Vector.h"

template <typename T>
class SpectralDecompositionTree : public TreeStructuredObject<SpectralDecomposition<T>> {
public:
	using TreeStructuredObject<SpectralDecomposition<T>>::attributes_;

	SpectralDecompositionTree() = default;

	explicit SpectralDecompositionTree(const TTBasis& basis);

	SpectralDecompositionTree(const HoleMatrixTree<T>& H, const TTBasis& basis);

	~SpectralDecompositionTree() = default;

	void Initialize(const TTBasis& basis);

	void Calculate(const HoleMatrixTree<T>& H, const TTBasis& basis);

	HoleMatrixTree<T> Invert(const TTBasis& basis, double eps = 1e-7);

	/// I/O
	void print(const TTBasis& basis) const;
};

typedef SpectralDecompositionTree<complex<double>> SpectralDecompositionTreecd;
typedef SpectralDecompositionTree<double> SpectralDecompositionTreed;


#endif //SPECTRALDECOMPOSITIONTREE_H
