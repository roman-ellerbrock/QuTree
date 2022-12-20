//
// Created by Roman on 3/19/2019.
//

#pragma once
#include "TreeOperators/Hamiltonian.h"
#include "Pauli.h"


/// create a rectangular ising model
SOPcd ising2D(size_t Lx, size_t Ly, double h);
