#pragma once
#include "TensorTree.h"
#include "FactorMatrix.h"
#include "TreeStructuredObject.h"

//! Calculate Overlap between two TensorTree objects with the same basis
/*!
This class can be used to calculate the overlap between two TensorTrees
<Psi_1|Psi_2>. The overlap matrix at each node are calculated, whether they
are perpendicular or not. Thus, the class is called "DenseOverlap".
If you look for Sparse-Overlaps, take a look at the H-Matrices class.
*/

template<typename T>
class DenseOverlap
	: public TreeStructuredObject<FactorMatrix<T>> {
public:
	using TreeStructuredObject<FactorMatrix<T>>::attributes;
	DenseOverlap() = default;

	DenseOverlap(istream& is);

	DenseOverlap(const string& filename);

	explicit DenseOverlap(const TTBasis& basis);

	DenseOverlap(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const TTBasis& basis);

	~DenseOverlap() = default;

	// Initialization routine
	void Initialize(const TTBasis& basis);

	// Calculate the overlap of two wavefunctions with the same basis
	FactorMatrix<T> Calculate(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const TTBasis& basis);

	// Calculate the overlap at a node
	void CalculateLayer(const Tensor<T>& Phi,
		Tensor<T> Chi, const Node& node);

	// Transform a Tensor using the sublying overlap matrices
	Tensor<T> TransformTensor(const Tensor<T>& Phi, const Node& node) const;

	// Get total overlap of the wavefunction (Calculate it first!)
	FactorMatrix<T>& Get() {
		assert(attributes.size() > 0);
		return attributes.back();
	}

	/// I/O
	void print(const TTBasis& basis, ostream& os = cout) const;
	void print(ostream& os = cout) const;

	void Write(ostream& os) const;
	void Read(istream& is);
};

template <typename T>
ostream& operator<<(ostream& os, const DenseOverlap<T>& S);

template <typename T>
istream& operator>>(istream& is, DenseOverlap<T>& S);

typedef DenseOverlap<complex<double>> DenseOverlapcd;
typedef DenseOverlap<complex<double>> DenseOverapd;

