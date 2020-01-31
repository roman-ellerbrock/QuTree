//
// Created by Roman on 3/29/2019.
//

#ifndef MCTDH_HOLEOVERLAP_H
#define MCTDH_HOLEOVERLAP_H
#include "DenseOverlap.h"

template<typename T>
class HoleOverlap: public TreeStructuredObject<Matrix<T>>
/**
 * \class HoleOverlap
 * \ingroup Tree
 * \brief Calculate Hole-Overlaps (A, B)_(p)
 *
 * This class calculates hole-overlaps and holds the resulting
 * Matrices. HoleOverlaps can be calculated from TensorTree
 * objects and corresponding DenseOverlap objects. Similar to
 * the DenseOverlap class, it works also for
 * non-orthogonal TensorTrees.
 */
{
public:
	using TreeStructuredObject<Matrix<T>>::attributes;

	/// Default Constructor
	HoleOverlap() = default;

	/// Construct from stream
	explicit HoleOverlap(istream& is);

	/// Construct from file
	explicit HoleOverlap(const string& filename);

	/// Construct and allocate memory for every node
	explicit HoleOverlap(const TTBasis& basis);

	/// Construct, allocate and calculate
	HoleOverlap(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const DenseOverlap<T>& S, const TTBasis& basis);

	/// Default destructor
	~HoleOverlap() = default;

	/// (Re-)allocate memory for every node
	void Initialize(const TTBasis& basis);

	/// Calculate hole-matrices locally
	void CalculateLayer(const TensorTree<T>& Psi,
		const TensorTree<T>& Chi, const DenseOverlap<T>& S,
		const Node& node);

	/// Calculate hole-matrices for the whole tree
	void Calculate(
		const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const DenseOverlap<T>& S, const TTBasis& basis);

	/// I/O
	/// Print in human readable format
	void print(const TTBasis& basis, ostream& os = cout) const;
	void print(ostream& os = cout) const;

	/// Write in binary format
	void Write(ostream& os) const;
	void Write(const string& filename) const;

	/// Read in binary format
	void Read(istream& is);
	void Read(const string& filename);
};

template<typename T>
ostream& operator<<(ostream& os, const HoleOverlap<T>& S);

template<typename T>
istream& operator>>(istream& is, HoleOverlap<T>& S);

typedef HoleOverlap<complex<double>> HoleOverlapcd;

typedef HoleOverlap<double> HoleOverlapd;

#endif //MCTDH_HOLEOVERLAP_H
