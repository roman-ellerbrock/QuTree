//
// Created by Roman Ellerbrock on 2019-07-18.
//

#include "NumberPartitioning_MLO.h"
#include "Pauli.h"


void NumberPartitioning_MLO::SpecialInitialize(const mctdhBasis& basis) {

	/*!
	 * The operator has the form
	 * H = (sum_i s_i*n_i)**2
	 * Following setup is used:
	 * Bottomlayer:
	 * {sigma_z, sigma_z**2}
	 * Upper layer:
	 * {z1+z2+..., z1**2+z2**2+z1*z2+ ...}
	 * Toplayer:
	 * {H_sep + v_coup + x_l*x_r}
	 */

	size_t f = basis.nLeaves();
	double A = 1. / (1. * f * f);
	vector<size_t> n;
	for (size_t l = 0; l < f; ++l) {
		n.push_back(l % 8 + 1);
//		n.push_back(l + 1);
	}
	for (const Node& node : basis) {
		if (node.IsBottomlayer()) {
			operator[](node) = InitializeBottom(node, n, A);
		} else {
			operator[](node) = InitializeUpper(node);
		}
	}
}

lSOPlist NumberPartitioning_MLO::InitializeUpper(const Node& node) {
	constexpr size_t SUM = 0;
	constexpr size_t sigma_z2 = 1;

	size_t d = node.nChildren();

	// create the sum of sigma_z's from underlying sums
	lSOP sum;
	for (size_t k = 0; k < d; ++k) {
		lSPO hk(SUM, k);
		sum.push_back({hk}, 1.);
	}

	// Local H=sum_ij z_i*z_j
	lSOP H_local;
	for (size_t k = 0; k < d; ++k) {
		for (size_t j = 0; j < d; ++j) {
			if (j != k) {
				lSPO hk(SUM, k);
				lSPO hj(SUM, j);
				lMPO M({hk, hj});
				H_local.push_back(M, 1.);
			}
		}
	}

	// Add H_local from children
	for (size_t k = 0; k < d; ++k) {
		lSPO hk(sigma_z2, k);
		H_local.push_back({hk}, 1.);
	}

	if (node.IsToplayer()) {
		return {H_local};
	} else {
		return {sum, H_local};
	}
}

lSOPlist NumberPartitioning_MLO::InitializeBottom(const Node& node,
	const vector<size_t>& n, double A) const {

	// Link to the base operators and set coefficients
	const Leaf& phy = node.PhysCoord();
	size_t mode = phy.Mode();
	constexpr size_t SUM = 0;
	constexpr size_t sigma_z2 = 3;

	// build two base operators
	lSOP sum;
	lSOP H_local;
	{
		lSPO z(SUM, mode);
		sum.push_back({z}, sqrt(A) * n[mode]);
	}
	{
		lSPO z2(sigma_z2, mode);
		H_local.push_back({z2}, A * n[mode] * n[mode]);
	}

	return {sum, H_local};
}

void NumberPartitioning_MLO::SpecialInitializeBottom(const mctdhBasis& basis) {
	function<void(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi)> sigma_x = PauliMatrices::sigma_x;
	function<void(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi)> sigma_y = PauliMatrices::sigma_y;
	function<void(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi)> sigma_z = PauliMatrices::sigma_z;
	function<void(const PrimitiveBasis& grid, Tensorcd& psi, const Tensorcd& phi)> Identity = PauliMatrices::Identity;

	for (size_t k = 0; k < basis.nLeaves(); ++k) {
		bottom_lib lib;
		vector<string> names;

		lib.emplace_back(make_shared<StandardRefSPO>(sigma_z));
		names.emplace_back("z");

		lib.emplace_back(make_shared<StandardRefSPO>(sigma_y));
		names.emplace_back("y");

		lib.emplace_back(make_shared<StandardRefSPO>(sigma_x));
		names.emplace_back("x");

		lib.emplace_back(make_shared<StandardRefSPO>(Identity));
		names.emplace_back("z2");

		push_back(lib, names);
	}
}
