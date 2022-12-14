//
// Created by Roman on 3/31/2019.
//

#include "CoupledHO.h"

SOPcd CoupledHO(const Tree& tree) {
	LeafFuncd x = &LeafInterface::applyX;
	LeafFuncd x2 = &LeafInterface::applyX2;
	LeafFuncd kin = &LeafInterface::applyKin;
	LeafFuncd p = &LeafInterface::applyP;
    LeafFuncd id = &LeafInterface::identity;

	constexpr double cm = 219474.6313705;
//	constexpr double lambda = 1000. / cm;
	constexpr double lambda = 2000. / cm;
	constexpr double omega = 4000. / cm;

	constexpr double c = 0.5 * omega * omega;

	size_t f = tree.nLeaves();

	SOPcd H;
// Kinetic energy and uncoupled HO
	for (size_t k = 0; k < f; ++k) {
		{
			MLOcd M;
			M.push_back(kin, k);
			H.push_back(M, 1.);
		}
		{
			MLOcd M;
			M.push_back(x2, k);
			H.push_back(M, c);
		}
	}

// coupling
	for (size_t k = 0; k < f; ++k) {
		size_t kn = (k + 1) % f;
		MLOcd M;
		M.push_back(x, k);
		M.push_back(x, kn);
		H.push_back(M, lambda * lambda);
	}

// normalize ground state to zero
    {
        MLOcd M;
        M.push_back(id,0);
        // determine ground state energy
        double gse = 0.;
        for(int i = 0; i < f; ++i){
            gse += sqrt(omega * omega + 2. * lambda * lambda * cos(2. * M_PI * (1.*i)/(1.*f)));
        }
        gse /= (-2.);
        H.push_back(M,gse);
    }


	return H;

	/// H = a * (h1 * h2 *...*hd) + b(g1 * g2 * ...)
}



