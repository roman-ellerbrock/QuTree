//
// Created by Roman Ellerbrock on 10/20/22.
//

#include "FermiGas.h"
#include "JordanWigner.h"

void fermiPlus(MLOcd& M, size_t idx) {
	M.push_back(JordanWigner::sigmaPlus(), idx);
	for (size_t k = 0; k < idx; ++k) {
		M.push_back(JordanWigner::sigmaZ(), k);
	}
}

void fermiMinus(MLOcd& M, size_t idx) {
	M.push_back(JordanWigner::sigmaMinus(), idx);
	for (size_t k = 0; k < idx; ++k) {
		M.push_back(JordanWigner::sigmaZ(), k);
	}
}

size_t idx2D(size_t i, size_t j, size_t sigma, size_t N) {
	size_t idx = i + j * N;
	if (sigma) { idx += N * N; }
	return idx;
}

MLOcd fermiHopping(size_t idx1, size_t idx2) {
	MLOcd M;
	fermiPlus(M, idx1);
	fermiMinus(M, idx2);
	return M;
}

MLOcd fermiHoppingAdj(size_t idx1, size_t idx2) {
	MLOcd M;
	fermiMinus(M, idx1);
	fermiPlus(M, idx2);
	return M;
}

SOPcd fermiGas2D(size_t N) {
	/// Create Fermi Gas Hamiltonian in 2D with open boundary conditions
	double eps = 1e0 / (N + 1e0);
	double eps2 = eps * eps;
	cout << "eps = " << eps << endl;
	cout << "eps2 = " << eps2 << endl;

	SOPcd H;
	/// Hopping terms
	for (size_t sigma = 0; sigma < 2; ++sigma) {
		/// <ij> hopping
		for (size_t i = 0; i < N; ++i) {
			for (size_t j = 0; j < N; ++j) {

				size_t idx = idx2D(i, j, sigma, N);
				if (i < N - 1) {
					size_t idx2 = idx2D(i + 1, j, sigma, N);
					H.push_back(fermiHopping(idx, idx2), 1e0 / eps2);
					H.push_back(fermiHoppingAdj(idx, idx2), 1e0 / eps2);
					cout << "(" << i << ", " << j << ")=" << idx
						 << " --> (" << i + 1 << ", " << j << ")=" << idx2
						 << ", sigma = " << sigma << "\n";
				}

				if (j < N - 1) {
					size_t idx2 = idx2D(i, j + 1, sigma, N);
					H.push_back(fermiHopping(idx, idx2), 1e0 / eps2);
					H.push_back(fermiHoppingAdj(idx, idx2), 1e0 / eps2);
					cout << "(" << i << ", " << j << ")=" << idx
						 << " --> (" << i << ", " << j + 1 << ")=" << idx2
						 << ", sigma = " << sigma << "\n";
				}
			}
		}

		/// ii kin
		for (size_t i = 0; i < N; ++i) {
			for (size_t j = 0; j < N; ++j) {
				size_t idx = idx2D(i, j, sigma, N);
				/// Adj shouldn't matter
				H.push_back(fermiHoppingAdj(idx, idx), -4e0 / eps2);
				cout << "(" << i << ", " << j << ") = " << idx << ", sigma = " << sigma << "\n";
			}
		}
	}
	return H;
}