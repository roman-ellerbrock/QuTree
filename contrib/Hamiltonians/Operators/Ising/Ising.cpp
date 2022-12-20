#include "Ising.h"
#include <random>
#include "Pauli.h"


size_t idx(size_t lx, size_t ly, size_t Lx, size_t Ly) {
	lx = lx % Lx;
	cout << "(" << lx << ", " << ly << ")";
	if (ly >= Ly) {
		cerr << "wrong ly index.\n";
		exit(1);
	}
	return ly * Lx + lx;
}

void sigma_x_1024(Tensorcd& hPhi, const Tensorcd& phi, size_t mode) {
	/// mode \in [0, 9]
	const TensorShape shape({2, 2, 2, 2, 2, 2, 2, 2, 2, 2});
	size_t before = shape.before(mode);
	size_t after = shape.after(mode);
	size_t beaft = before * after;
	for (size_t n = 0; n < shape.after(mode); n++) {
		for (size_t bef = 0; bef < shape.before(mode); ++bef) {
			hPhi[bef + n * beaft] = phi[bef + before + n * beaft];
		}
		for (size_t bef = 0; bef < shape.before(mode); ++bef) {
			hPhi[bef + before + n * beaft] = phi[bef + n * beaft];
		}
	}
}



void subH(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& A) {
	size_t Ly = 4;
	hA.zero();
	for (size_t l = 0; l < Ly; ++l) {

	}
}

SOPcd ising2D(size_t Lx, size_t Ly, double h) {
	SOPcd S;
	for (size_t lx = 0; lx < Lx; ++lx) {
		/// - X_mn * X_(m+1)n

		for (size_t ly = 0; ly < Ly; ++ly) {
			{
				MLOcd N;
				cout << "-X";
				N.push_back(PauliMatrices::sigma_x, idx(lx, ly, Lx, Ly));
				cout << " * X";
				N.push_back(PauliMatrices::sigma_x, idx((lx + 1) % Lx, ly, Lx, Ly));
				cout << endl;
				S.push_back(N, -1.);
			}

			/// - X_mn * X_m(n+1)
			if (ly < Ly - 1) {
				MLOcd M;
				cout << "-X";
				M.push_back(PauliMatrices::sigma_x, idx(lx, ly, Lx, Ly));
				cout << " * X";
				M.push_back(PauliMatrices::sigma_x, idx(lx, ly + 1, Lx, Ly));
				cout << endl;
				S.push_back(M, -1.);
			}

			///
			{
				MLOcd Z;
				cout << "+hZ";
				Z.push_back(PauliMatrices::sigma_z, idx(lx, ly, Lx, Ly));
				S.push_back(Z, h);
				cout << endl;
			}
			cout << endl;
		}
	}
	getchar();
	return S;
}

