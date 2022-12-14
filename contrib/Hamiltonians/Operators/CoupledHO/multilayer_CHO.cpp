//
// Created by Roman Ellerbrock on 2019-07-13.
//

#include "multilayer_CHO.h"

void multilayer_CHO::SpecialInitialize(const mctdhBasis& basis) {
	constexpr double cm = 219474.6313705;
	constexpr double lambda = 2000. / cm;
	constexpr double omega = 4000. / cm;
	constexpr double v = 0.5 * omega * omega;

	/*!
	 * Following setup is used:
	 * Bottomlayer:
	 * {H_sep, x}
	 * Upper layer:
	 * {H_sep, x_l, x_r, v_coup}
	 * Toplayer:
	 * {H_sep + v_coup + x_l*x_r}
	 */
	bool couple = true;

	for (const Node& node : basis) {
		if (node.IsBottomlayer()) {
			operator[](node) = BuildHamiltonianBottom(node, v, lambda, couple);
		} else {
			operator[](node) = BuildHamiltonianUpper(node, couple);
		}
	}
}

/*	for (const GetNode& node : basis) {
		size_t d = node.nChildren();
		if (node.IsBottomlayer()) {
			// Link to the base operators and set coefficients
			const Leaf& phy = node.PhysCoord();
			size_t mode = phy.Mode();
			lSPO T(0, mode);
			lSPO x2(1, mode);
			lSPO x(2, mode);
			// build two base operators
			lSOP H_sep;
			H_sep.push_back({T}, 1.);
			H_sep.push_back({x2}, v);

			if (couple) {
				lSOP x_left;
				x_left.push_back({x}, lambda);
				lSOP x_right;
				x_right.push_back({x}, lambda);
				// the list on this node is the seperabel and the coupling potential
				attributes.push_back({H_sep, x_left, x_right});
			} else {
				attributes.push_back({H_sep});
			}
		} else {
			// Add the seperable hamiltonians
			lSOP H_sep;
			for (size_t k = 0; k < d; ++k) {
				lSPO hk(Hs, k);
				H_sep.push_back({hk}, 1.);
			}

			if (couple) {
				// Add left- and right-most x-operators
				lSOP x_left;
				lSOP x_right;
				{
					lSPO xl(XL, 0);
					lSPO xr(XR, d - 1);
					x_left.push_back({xl}, 1.);
					x_right.push_back({xr}, 1.);
				}

				// Create the local coupling potential
				// 1.) the sublying coupling potentials
				lSOP local_v_coup;
				for (size_t k = 0; k < d; ++k) {
					const GetNode& child = node.Down(k);
					if (child.IsBottomlayer()) { continue; }
					lSPO vc_sub(Vc, k);
					local_v_coup.push_back({vc_sub}, 1.);
				}

				// 2.) the missing products of x's: \sum_k x_k * x_{k+1}
				for (size_t k = 0; k < d - 1; ++k) {
					// get the "right" x of the first child
					lSPO xk(XR, k);
					// get the "left" x of the next child
					lSPO xknext(XL, k + 1);
					lMPO xx = {xk, xknext};
					local_v_coup.push_back(xx, 1.);
				}

				// 3.) For toplayer: x_0 * x_{d-1}
				if (node.IsToplayer()) {
					lSPO x0(XL, 0);
					// get the "left" x of the next child
					lSPO xlast(XR, d - 1);
					lMPO xx = {x0, xlast};
					local_v_coup.push_back(xx, 1.);
				}

				if (node.IsToplayer()) {
					attributes.push_back({H_sep, local_v_coup});
				} else {
					attributes.push_back({H_sep, x_left, x_right, local_v_coup});
				}
			}
		}
		*/

lSOPlist multilayer_CHO::BuildHamiltonianUpper(const Node& node, bool coupling)const {
	constexpr size_t Hs = 0;
	constexpr size_t XL = 1;
	constexpr size_t XR = 2;
	constexpr size_t Vc = 3;
	// Add the seperable hamiltonians
	lSOP H_sep;
	size_t d = node.nChildren();
	for (size_t k = 0; k < d; ++k) {
		lSPO hk(Hs, k);
		H_sep.push_back({hk}, 1.);
	}

	if (coupling) {
		// Add left- and right-most x-operators
		lSOP x_left;
		lSOP x_right;
		{
			lSPO xl(XL, 0);
			lSPO xr(XR, d - 1);
			x_left.push_back({xl}, 1.);
			x_right.push_back({xr}, 1.);
		}

		// Create the local coupling potential
		// 1.) the sublying coupling potentials
		lSOP local_v_coup;
		for (size_t k = 0; k < d; ++k) {
			const Node& child = node.Down(k);
			if (child.IsBottomlayer()) { continue; }
			lSPO vc_sub(Vc, k);
			local_v_coup.push_back({vc_sub}, 1.);
		}

		// 2.) the missing products of x's: \sum_k x_k * x_{k+1}
		for (size_t k = 0; k < d - 1; ++k) {
			// get the "right" x of the first child
			lSPO xk(XR, k);
			// get the "left" x of the next child
			lSPO xknext(XL, k + 1);
			lMPO xx = {xk, xknext};
			local_v_coup.push_back(xx, 1.);
		}

		// 3.) For toplayer: x_0 * x_{d-1}
		if (node.IsToplayer()) {
			lSPO x0(XL, 0);
			// get the "left" x of the next child
			lSPO xlast(XR, d - 1);
			lMPO xx = {x0, xlast};
			local_v_coup.push_back(xx, 1.);
		}

		if (node.IsToplayer()) {
			return {H_sep, local_v_coup};
		} else {
			return {H_sep, x_left, x_right, local_v_coup};
		}
	} else {
		return {H_sep};
	}
}

lSOPlist multilayer_CHO::BuildHamiltonianBottom(const Node& node, double v_coeff,
	double lambda, bool coupling)const {

	// Link to the base operators and set coefficients
	const Leaf& phy = node.PhysCoord();
	size_t mode = phy.Mode();
	lSPO T(0, mode);
	lSPO x2(1, mode);
	lSPO x(2, mode);
	// build two base operators
	lSOP H_sep;
	H_sep.push_back({T}, 1.);
	H_sep.push_back({x2}, v_coeff);

	if (coupling) {
		lSOP x_left;
		x_left.push_back({x}, lambda);
		lSOP x_right;
		x_right.push_back({x}, lambda);
		// the list on this node is the seperabel and the coupling potential
		return {H_sep, x_left, x_right};
	} else {
		return {H_sep};
	}

}

void multilayer_CHO::SpecialInitializeBottom(const mctdhBasis& basis) {
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x2 = &PrimitiveBasis::ApplyX2;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> kin = &PrimitiveBasis::ApplyKin;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> p = &PrimitiveBasis::ApplyP;

	for (size_t k = 0; k < basis.nLeaves(); ++k) {
		bottom_lib lib;
		vector<string> names;
		// kinetic energy
		{
			lib.emplace_back(make_shared<StandardSPO>(kin));
			names.emplace_back("T");
		}

		{
			// V1 = x**2
			lib.emplace_back(make_shared<StandardSPO>(x2));
			names.emplace_back("x**2");
		}

		{
			// V2 = x
			lib.emplace_back(make_shared<StandardSPO>(x));
			names.emplace_back("x");
		}

		push_back(lib, names);
	}
}

