//
// Created by Roman on 3/19/2019.
//

#pragma once

#include "Hamiltonian.h"
#include "FermionNumberBasis.h"
#include "Graph.h"
#include <chrono>
#include "Pauli.h"

class Ising :
	public SOP
{
public:
	Ising(const mctdhBasis& basis) {}
	Ising(const mctdhBasis& basis, size_t N, size_t method);
    Ising(const mctdhBasis& basis, double J, double H) {
        SpinChain(basis, H, H);
    };
	~Ising() = default;

	void HamiltonianPaths(const mctdhBasis & basis, size_t N, size_t method);
    void SpinChain(const mctdhBasis& basis, double J, double H);

};



