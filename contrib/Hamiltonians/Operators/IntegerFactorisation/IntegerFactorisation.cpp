//
// Created by Roman Ellerbrock on 2019-07-12.
//

#include "IntegerFactorisation.h"


IntegerFactorisation::IntegerFactorisation(const mctdhBasis& basis, size_t F) {
	SpecialInitialize(basis, F);
}

void IntegerFactorisation::SpecialInitialize(const mctdhBasis& basis, const size_t F) {
    size_t f  = basis.nLeaves();
    SOP S;
    for (size_t i = 0; i < f; ++i) {
        MultiParticleOperator M(PauliMatrices::bit_x, i);
        auto factor = pow(2., i);
        S.push_back(M, factor);
    }
    S = S * S;
    MPO Fop(PauliMatrices::Identity, 0);
    S.push_back(Fop, -F);
    S = S * S;
    for (size_t k = 0; k < S.size(); ++k) {
        push_back(S[k], S.Coeff(k));
    }
}