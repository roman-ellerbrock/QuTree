//
// Created by Roman Ellerbrock on 2020-01-23.
//

#ifndef TENSORTREE_H
#define TENSORTREE_H
#include "NodeAttribute.h"
#include "TreeShape/Tree.h"
#include "Core/Tensor.h"
#include <random>

/**
 * \defgroup Tree
 * \brief This group contains all functionality pertaining to Trees.
 */

template<typename T>
class TensorTree:
	public NodeAttribute<Tensor<T>>
	/**
	 * \class TensorTree
	 * \ingroup Tree
	 * \brief This class represents tensor trees.
	 *
	 * Usage:
	 * TRBasis tree(12, 2, 2)
	 * // Create Tensor with Zero-entry tensors at every node
	 * TensorTreecd Psi(tree);
	 * Psi.Write("filename.dat");
	 * TensorTreecd Chi("filename.dat");
	 */
{
public:
	using NodeAttribute<Tensor<T>>::attributes_;

	/// Default constructor without memory allocation
	TensorTree() = default;
	/// Constructor with allocation of memory
	explicit TensorTree(const Tree& tree);

	/// Construct TensorTree from stream
	explicit TensorTree(istream& is);

	/// Construct TensorTree from file
	explicit TensorTree(const string& filename);

	/// Create tensor tree and occupy the coefficients
	TensorTree(mt19937& gen, const Tree& tree, bool delta_lowest = true);

	/// Default destructor
	~TensorTree() = default;

	/// Create Tensors for all nodes
	virtual void Initialize(const Tree& tree);

	/// Generate TTs
	void FillRandom(mt19937& gen, const Tree& tree, bool delta_lowest = true);

	/// (File) I/O
	/// Read TensorTree from stream (binary format)
	void Read(istream& is);

	/// Write TensorTree to stream (binary format)
	void Write(ostream& os) const;

	/// Write TensorTree to file (binary format)
	void Write(const string& filename) const;

	/// Print info in human readable format
	void print(const Tree& tree, ostream& os = cout) const;

	////////////////////////////////////////////////////////////////////////
	/// Arithmetic operators
	////////////////////////////////////////////////////////////////////////
	TensorTree& operator+=(const TensorTree<T>& R) {
		for (size_t n = 0; n < attributes_.size(); ++n) {
			attributes_[n] += R.attributes_[n];
		}
		return *this;
	}

	TensorTree& operator-=(const TensorTree<T>& R) {
		for (size_t n = 0; n < attributes_.size(); ++n) {
			attributes_[n] -= R.attributes_[n];
		}
		return *this;
	}

	TensorTree operator+(const TensorTree<T>& R) {
		TensorTree<T> S(*this);
		S += R;
		return S;
	}

	TensorTree operator-(const TensorTree<T>& R) {
		TensorTree<T> S(*this);
		S -= R;
		return S;
	}

	void operator*=(T c) {
		for (auto& A : *this) {
			A *= c;
		}
	}

	void operator/=(T c) {
		for (auto& A : *this) {
			A /= c;
		}
	}

	TensorTree operator*(T c) {
		TensorTree<T> S(*this);
		S *= c;
		return S;
	}

	TensorTree operator/(T c) {
		TensorTree<T> S(*this);
		S /= c;
		return S;
	}

protected:
	void FillBottom(Tensor<T>& Phi, const Node& node);
	void FillUpper(Tensor<T>& Phi, std::mt19937& gen,
		const Node& node, bool delta_lowest = true);
};

typedef TensorTree<complex<double>> TensorTreecd;

typedef TensorTree<double> TensorTreed;

template <typename T>
void Orthogonal(TensorTree<T>& Psi, const Tree& tree);

template <typename T>
void QROrthogonal(TensorTree<T>& Psi, const Tree& tree);

template <typename T>
void Orthonormal(TensorTree<T>& Psi, const Tree& tree);

template<typename T>
ostream& operator<<(ostream& os, const TensorTree<T>& t);
template<typename T>
istream& operator>>(istream& is, TensorTree<T>& t);

template <typename T>
TensorTree<T> operator*(T c, TensorTree<T> R) {
	R *= c;
	return R;
}

template <typename T>
TensorTree<T> operator/(T c, TensorTree<T> R) {
	R /= c;
	return R;
}

#endif //TENSORTREE_H
