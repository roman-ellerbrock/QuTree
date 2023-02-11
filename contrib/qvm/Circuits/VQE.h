//
// Created by Shaun Datta on 7/27/21.
//

#ifndef QUTREEVM_VQE_H
#define QUTREEVM_VQE_H
#include "TreeOperators/SOPVector.h"

typedef SOPVectorcd Circuit;
/// makes circuit with staircase of CNOTs
Circuit exampleCirc(size_t num_qubits);

/**
 * IsingCirc is a benchmark 1D entanglement circuit with interleaved CNOTs and Hadamards
 * @param num_qubits
 * @param depth
 * @return
 */

//Circuit IsingCirc(size_t num_qubits, size_t depth);

/**
 * addLayer generates a single layer of IsingCirc.
 * @param num_qubits
 * @param offset specifies whether the layer is odd- or even-numbered,
 * i.e. whether the top gate should act on the 0th or 1st qubit
 * @return
 */
Circuit addLayer(size_t num_qubits, bool offset, double theta=0, bool adjoint = false);

Circuit paramCirc(size_t num_qubits, size_t depth, vector<Vectord> params, bool adjoint);

void vqe(TensorTreecd& psi, const string& filename, const Tree& tree, const size_t depth, size_t max_iter);

Vectord updateParameters(const Vectord& params, const Vectord& grad, double step);
Vectord getGradient(const Vectord& params, const Vectord& eps, double E);
double getH(Vectord params);
Vectord flatten(const vector<Vectord>& vec, int num_params, const size_t depth);
vector<Vectord> reshape(const Vectord& vec, const size_t depth, const size_t length);
size_t size(const vector<Vectord>& x);
vector<Vectord> initializeVQEparameters(mt19937& gen, size_t depth, const Tree& tree);

#endif //QUTREEVM_VQE_H
