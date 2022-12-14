//
// Created by Roman Ellerbrock on 2019-07-18.
//

#ifndef NUMBERPARTITIONING_MLO_H
#define NUMBERPARTITIONING_MLO_H
#include "multilayer_Operator.h"

class NumberPartitioning_MLO: public MLO {
public:
	explicit NumberPartitioning_MLO(const mctdhBasis& basis)
		: MLO(basis) { Initialize(basis); }

	~NumberPartitioning_MLO() = default;

private:
	void SpecialInitialize(const mctdhBasis& basis) override;

	lSOPlist InitializeUpper(const Node& node);

	lSOPlist InitializeBottom(const Node& node,
		const vector<size_t>& n, double A) const;

	void SpecialInitializeBottom(const mctdhBasis& basis) override;
};


#endif //NUMBERPARTITIONING_MLO_H
