#pragma once
#include "TensorTree.h"
#include "Core/FactorMatrix.h"
#include "TreeStructuredObject.h"

template<typename T>
class DenseOverlap
	: public TreeStructuredObject<FactorMatrix<T>>
	/**
	 * \class DenseOverlap
	 * \ingroup Tree
	 * \brief Calculate the Overlap between non-orthogonal TensorTrees
	 *
	 *
	 */
{
public:
	using TreeStructuredObject<FactorMatrix<T>>::attributes;
	/// Default constructor
	DenseOverlap() = default;

	/// Construct from stream
	explicit DenseOverlap(istream& is);

	/// Construct from file
	explicit DenseOverlap(const string& filename);

	/// Construct and allocate memory for every node
	explicit DenseOverlap(const TTBasis& basis);

	/// Construct, allocate and calculate
	DenseOverlap(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const TTBasis& basis);

	/// Default destructor
	~DenseOverlap() = default;

	/// Allocate a FactorMatrix at every node
	void Initialize(const TTBasis& basis);

	/// Calculate the tensor tree dot-product (Psi, Chi)_p
	FactorMatrix<T> Calculate(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const TTBasis& basis);

	/// Calculate the local tensor tree dot-product
	void CalculateLayer(const Tensor<T>& Phi,
		Tensor<T> Chi, const Node& node);

	/// Perform the (local) FactorMatrix tree - tensor tree product
	Tensor<T> TransformTensor(const Tensor<T>& Phi, const Node& node) const;

	/// Get FactorMatrix at Toplayer
	FactorMatrix<T>& Get() {
		assert(attributes.size() > 0);
		return attributes.back();
	}

	/// I/O
	/// Print human readable
	void print(const TTBasis& basis, ostream& os = cout) const;
	void print(ostream& os = cout) const;

	/// Write in binary format
	void Write(ostream& os) const;
	/// Read in binary format
	void Read(istream& is);
};

template<typename T>
ostream& operator<<(ostream& os, const DenseOverlap<T>& S);

template<typename T>
istream& operator>>(istream& is, DenseOverlap<T>& S);

typedef DenseOverlap<complex<double>> DenseOverlapcd;

typedef DenseOverlap<complex<double>> DenseOverapd;

