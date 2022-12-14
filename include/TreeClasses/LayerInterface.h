//
// Created by Roman Ellerbrock on 3/8/20.
//

#ifndef LAYERINTERFACE_H
#define LAYERINTERFACE_H
#include "TreeClasses/HamiltonianRepresentation.h"

class LayerInterface {
public:
	LayerInterface(const Hamiltonian& H, const HamiltonianRepresentation& hRep,
		const Node& node, complex<double> propagation_phase)
		: h_(&H), hRep_(&hRep), node_(&node), propagation_phase_(propagation_phase) {}

	~LayerInterface() = default;

	void Derivative(double time, Tensorcd& dPhi, const Tensorcd& Phi) {
		LayerDerivative(dPhi, time, Phi, *h_, *hRep_, *node_, propagation_phase_);
	}

	double Error(const Tensorcd& Phi, const Tensorcd& Chi) const;

	complex<double> propagation_phase_;
private:
	const Hamiltonian* h_;
	const HamiltonianRepresentation* hRep_;
	const Node* node_;
};

#endif //LAYERINTERFACE_H
