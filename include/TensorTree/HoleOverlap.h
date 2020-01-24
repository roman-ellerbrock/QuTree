//
// Created by Roman on 3/29/2019.
//

#ifndef MCTDH_HOLEOVERLAP_H
#define MCTDH_HOLEOVERLAP_H
#include "DenseOverlap.h"

template<typename T>
class HoleOverlap: public TreeStructuredObject<Matrix<T>> {
public:
	using TreeStructuredObject<Matrix<T>>::attributes;

	HoleOverlap() = default;
	~HoleOverlap() = default;

	HoleOverlap(const TTBasis& basis);

	HoleOverlap(
		const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const DenseOverlap<T>& S, const TTBasis& basis);

	void Initialize(const TTBasis& basis);

	void CalculateLayer(const TensorTree<T>& Psi,
		const TensorTree<T>& Chi, const DenseOverlap<T>& S,
		const Node& node);

	void Calculate(
		const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const DenseOverlap<T>& S, const TTBasis& basis);

	void print(const TTBasis& basis, ostream& os = cout) const;
};

typedef HoleOverlap<complex<double>> HoleOverlapcd;
typedef HoleOverlap<double> HoleOverlapd;

#endif //MCTDH_HOLEOVERLAP_H
