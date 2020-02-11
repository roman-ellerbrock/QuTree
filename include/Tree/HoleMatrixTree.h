//
// Created by Roman on 3/29/2019.
//

#ifndef MCTDH_HOLEOVERLAP_H
#define MCTDH_HOLEOVERLAP_H
#include "FactorMatrixTree.h"

template<typename T>
class HoleMatrixTree: public TreeStructuredObject<FactorMatrix<T>>
/**
 * \class HoleMatrixTree
 * \ingroup Tree
 * \brief Calculate Hole-Overlaps (A, B)_(p_)
 *
 * This class calculates hole-overlaps and holds the resulting
 * Matrices. HoleMatrixTrees can be calculated from TensorTree
 * objects and corresponding FactorMatrixTree objects. Similar to
 * the FactorMatrixTree class, it works also for
 * non-orthogonal TensorTrees.
 */
{
public:
	using TreeStructuredObject<Matrix<T>>::attributes_;

	/// Default Constructor
	HoleMatrixTree() = default;

	/// Construct from stream
	explicit HoleMatrixTree(istream& is);

	/// Construct from file
	explicit HoleMatrixTree(const string& filename);

	/// Construct and allocate memory for every node
	explicit HoleMatrixTree(const TTBasis& basis);

	/// Construct, allocate and calculate
	HoleMatrixTree(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const FactorMatrixTree<T>& S, const TTBasis& basis);

	/// Construct, allocate and calculate non-othorgonal
	HoleMatrixTree(const TensorTree<T>& Psi, const TTBasis& basis);

	/// Default destructor
	~HoleMatrixTree() = default;

	/// (Re-)allocate memory for every node
	void Initialize(const TTBasis& basis);

	/// calculate hole-matrices locally, non-othorgonal
	void CalculateLayer(const TensorTree<T>& Psi,
		const TensorTree<T>& Chi, const FactorMatrixTree<T>& S,
		const Node& node);

	/// Calculate hole-matrices for the whole tree, non-othorgonal
	void Calculate(
		const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const FactorMatrixTree<T>& S, const TTBasis& basis);

	/// Calculate hole-matrices assuming orthogonal tensor tree
	void Calculate(
		const TensorTree<T>& Psi, const TTBasis& basis);

	/// Calculate hole-matrices locally
	void CalculateLayer(const TensorTree<T>& Psi, const Node& node);

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
ostream& operator<<(ostream& os, const HoleMatrixTree<T>& S);

template<typename T>
istream& operator>>(istream& is, HoleMatrixTree<T>& S);

typedef HoleMatrixTree<complex<double>> HoleMatrixTreecd;

typedef HoleMatrixTree<double> HoleMatrixTreed;

#endif //MCTDH_HOLEOVERLAP_H
