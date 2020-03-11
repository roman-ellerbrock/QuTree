//
// Created by Roman Ellerbrock on 3/10/20.
//

#ifndef SOPMATRIXTREES_H
#define SOPMATRIXTREES_H
#include "TreeClasses/SparseMatrixTree.h"
#include "TreeOperators/SumOfProductsOperator.h"

template <typename T>
class SOPMatrixTrees {
public:
	SOPMatrixTrees(const SOP<T>& H, const Tree& tree) {
		for (const auto& M : H) {
			matrices_.push_back(SparseMatrixTree<T>(M, tree));
			contractions_.push_back(SparseMatrixTree<T>(M, tree));
		}
	}

	void print() const {
		for (size_t k = 0; k < matrices_.size(); ++k) {
			cout << "MatrixTree[" << k << "]:\n";
			matrices_[k].print();
		}
		for (size_t k = 0; k < contractions_.size(); ++k) {
			cout << "ContractionTree[" << k << "]:\n";
			contractions_[k].print();
		}
	}

	size_t size() const {
		assert(matrices_.size() == contractions_.size());
		return matrices_.size();
	}

	vector<SparseMatrixTree<T>> matrices_;
	vector<SparseMatrixTree<T>> contractions_;
};

typedef SOPMatrixTrees<complex<double>> MatrixTreescd;
typedef SOPMatrixTrees<double> MatrixTreesd;


#endif //SOPMATRIXTREES_H
