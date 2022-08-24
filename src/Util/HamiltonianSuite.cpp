#include "Util/HamiltonianSuite.h"

SOPcd coupledHarmonicOssciallator(size_t nLeaves) {

    LeafFuncd p(&BasisAPI::applyP);
    LeafFuncd kin(&BasisAPI::applyKin);
    LeafFuncd x(&BasisAPI::applyX);
    LeafFuncd x2(&BasisAPI::applyX2);

    double cm = 219474.6313705;
    double lambda = 2000. / cm;
    double omega = 4000. / cm;
    double c = 0.5 * omega * omega;

    SOPcd H;
    for (size_t k = 0; k < nLeaves; ++k) {
        ProductOperatorcd T;
        T.push_back(kin, k);
        H.push_back(T, 1.);

        ProductOperatorcd V;
        V.push_back(x2, k);
        H.push_back(V, c);
    }

    for (size_t k = 0; k < nLeaves; ++k) {
        size_t kn = k + 1 % nLeaves;
        ProductOperatorcd M;
        M.push_back(x, k);
        M.push_back(x, kn);
        H.push_back(M, lambda * lambda);
    }

    return H;
}
