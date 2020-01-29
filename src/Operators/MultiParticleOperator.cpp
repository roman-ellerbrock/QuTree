#include "MultiParticleOperator.h"


template <typename T>
MultiParticleOperator<T>::MultiParticleOperator()
	: hasV(false) {
	mode_.clear();
}

template <typename T>
Tensor<T> MultiParticleOperator<T>::ApplyBottomLayer(Tensor<T> Phi,
	const Leaf& Phy) const {
	Tensor<T> hPhi(Phi.Dim());
	size_t mode_x = Phy.Mode();
	const PrimitiveBasis& grid = Phy.PrimitiveGrid();
	bool switchbool = true;

	// Applying the MPO uses switching of the result Tensor to increase performance.
	for (size_t l = 0; l < SingParOp.size(); ++l) {
		if (mode_x != mode_[l]) { continue; }

		shared_ptr<SPO<T>> spo = operator[](l);
		// apply it
		if (switchbool) {
			spo->Apply(grid, hPhi, Phi);
		} else {
			spo->Apply(grid, Phi, hPhi);
		}
		switchbool = !switchbool;
	}

	if (switchbool) {
		return Phi;
	} else {
		return hPhi;
	}
}

template <typename T>
Tensor<T> MultiParticleOperator<T>::ApplyBottomLayer(Tensor<T> Acoeffs,
	const vector<int>& list, const PrimitiveBasis& grid) const {
	Tensor<T> hAcoeff(Acoeffs.Dim());
	bool switchbool = true;
	// Applying the MPO uses switching of the result Tensor to increase performance.
	for (size_t l = 0; l < list.size(); l++) {
		// get the active part in the MPO
		int part = list[l];
		shared_ptr<SPO<T>> spo = operator[](part);

		// apply it
		if (switchbool) {
			spo->Apply(grid, hAcoeff, Acoeffs);
		} else {
			spo->Apply(grid, Acoeffs, hAcoeff);
		}
		switchbool = !switchbool;
	}

	if (switchbool) {
		return Acoeffs;
	} else {
		return hAcoeff;
	}
}

template<typename T>
TensorTree<T> MultiParticleOperator<T>::Apply(TensorTree<T> Psi,
	const TTBasis& basis) const {
	for (size_t i = 0; i < basis.nNodes(); i++) {
		const Node& node = basis.GetNode(i);
		if (node.IsBottomlayer()) {
			const Leaf& phy = node.PhysCoord();
			const PrimitiveBasis& grid = phy.PrimitiveGrid();
			int coord = phy.Mode();

			// build list with active parts
			vector<int> activelayerparts;
			for (size_t l = 0; l < SingParOp.size(); l++) {
				if (mode_[l] == coord) {
					activelayerparts.push_back(l);
				}
			}

			Tensor<T>& Acoeff = Psi[node];
			Acoeff = ApplyBottomLayer(Acoeff, activelayerparts, grid);
		}
	}
	return Psi;
}

template <typename T>
void MultiParticleOperator<T>::SetV(const PotentialOperator& V_) {
	v = V_;
	hasV = true;
}

template class MultiParticleOperator<complex<double>>;
template class MultiParticleOperator<double>;
