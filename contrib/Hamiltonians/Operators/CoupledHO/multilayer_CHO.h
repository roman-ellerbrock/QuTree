//
// Created by Roman Ellerbrock on 2019-07-13.
//

#ifndef MULTILAYER_CHO_H
#define MULTILAYER_CHO_H
#include "multilayer_Operator.h"

class multilayer_CHO: public MLO {
public:
	explicit multilayer_CHO(const mctdhBasis& basis) : MLO(basis) { Initialize(basis); }

	~multilayer_CHO() = default;

private:
	void SpecialInitialize(const mctdhBasis& basis) override;
	void SpecialInitializeBottom(const mctdhBasis& basis) override;
	lSOPlist BuildHamiltonianUpper(const Node& node, bool coupling)const;

	lSOPlist BuildHamiltonianBottom(const Node& node, double v_coeff,
		double lambda, bool coupling)const;
};


#endif //MULTILAYER_CHO_H
