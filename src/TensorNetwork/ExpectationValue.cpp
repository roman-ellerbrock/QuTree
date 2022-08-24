#include "TensorNetwork/ExpectationValue.h"
#include "TensorNetwork/contractions.h"

Tensorcd expectationValue(const TensorTreecd& Psi, const SOPcd& H, const Tree& tree) {
    const Node& root = tree.root();
    Edge edge(&root, (Node*) nullptr);

    vector<TensorTreecd> Hs = contraction(Psi, Psi, H);

    Tensorcd hA = Psi[root];
    apply(hA, Hs, H, &edge);
    auto exp = contraction(Psi[root], hA, edge);

    return exp;
}
