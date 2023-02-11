//
// Created by Shaun Datta on 7/27/21.
//

#include "VQE.h"
#include "GateOperators.h"
#include "TreeClasses/MatrixTensorTreeFunctions.h"
#include "Applications/TreeApplyOperator.h"
#include "VQE/ElectronicStructure.h"
#include "VQE/JordanWigner.h"
#include "VQE/SOPexpectedvalue.h"

using namespace Circuits;

Circuit exampleCirc(size_t num_qubits) {
    Circuit circ;
    MLOcd H;
    H.push_back(Hadamard, 0);
    SOPcd SH;
    SH.push_back(H, 1.);
    circ.push_back(SH);

    for (int i = 0; i < num_qubits - 1; i++) {
        circ.push_back(CNot(i, i + 1));
    }
    return circ;
}

//overloaded addLayer for parameterized circuits
Circuit addLayer(size_t num_qubits, bool offset, Vectord params, bool adjoint) {
    Circuit circ;
    for (int i = offset; i < num_qubits - 2; i += 3) {
//        cout << "params index: " << (i - offset) / 3 << " | " << i << "  " << i + 1 << "  " << i + 2 << endl;
        circ.push_back(cGivens(params[(i - offset) / 3], i, i + 1, i + 2));
//        MLOcd H;
//        H.push_back(Hadamard, i);
//        H.push_back(Hadamard, i+1);
//        LeafMatrixcd lm(givens_rot(0.), adjoint);
//        H.push_back(lm,i);
//        SOPcd SH;
//        SH.push_back(H, 1.);
//        circ.push_back(SH);
    }
    return circ;
}

/*Circuit IsingCirc(size_t num_qubits, size_t depth){
    Circuit circ;

    for (int layer=0; layer<depth; layer++){
        bool offset = layer%2==0;
        circ.append(addLayer(num_qubits,offset));
    }

    return circ;
}*/

Circuit paramCirc(size_t num_qubits, size_t depth, vector<Vectord> params, bool adjoint) {
    Circuit circ;

    for (int layer = 0; layer < depth; layer++) {
        bool offset = layer % 3;
        circ.append(addLayer(num_qubits, offset, params[layer], adjoint));
    }

    return circ;
}

/**
 * update parameters based: params' = params - step * gradient of <H>
 */
Vectord updateParameters(const Vectord& params, const Vectord& grad, double step){
    Vectord newParams(params.dim());
    newParams = (params-grad*step);
    return newParams;
}

/**
 * k'th derivative of <H> = (<H>(params+eps_k) - <H>(params))/|eps_k|
 * eps_k is a vector of zeroes but for the k'th entry
 */
Vectord getGradient(const Vectord& params, const Vectord& eps, double E){
    int len = params.dim();
    Vectord grad(len);
    for (int i = 0; i < len; i++) {
        Vectord eps_k(len); // initialize to zero
        eps_k[i] = eps[i];
        double newE = getH(params+eps_k);
        grad[i] = ((newE-E)/abs(eps_k[i]));
    }
    return grad;
}

// scarecrow method for getH, which computes <H> given a set of parameters
double getH(Vectord params){
    return 0.;
}

/**
 * flattens an array of vectors into a single vector
 */
Vectord flatten(const vector<Vectord>& vec, int num_params, const size_t depth){
    Vectord flattened(num_params);
    int index = 0;
    for (int i=0; i<depth; ++i){
        for(int j=0; j<vec[i].dim(); ++j){
            flattened[index] = vec[i][j];
            ++index;
        }
    }
    return flattened;
}

/**
 * reshapes a single vector into vector<Vectord>, corresponding to
 * the circuit parameters organized layer-by-layer
 * vector<Vectord> has depth-many entries of length-long vectors
 */
vector<Vectord> reshape(const Vectord& vec, const size_t depth, const size_t length){
    vector<Vectord> params;
    int index = 0;
    for (int j = 0; j < depth; ++j) {
       Vectord layer_params(length);
       for (int k=0; k<length; ++k){
           layer_params[k]=vec[index];
           ++index;
       }
       params.emplace_back(layer_params);
    }
    return params;
}

size_t size(const vector<Vectord>& x) {
	int num_params = 0;
	for (const auto& y : x) {
		num_params += y.dim();
	}
	return num_params;
}

vector<Vectord> initializeVQEparameters(mt19937& gen, size_t depth, const Tree& tree) {
	vector<Vectord> params; // params is depth-number of vectors of length n/3

	uniform_real_distribution<double> unif(0., 1e-5);
	normal_distribution<double> gauss(0., 1e-5);

	for (int j = 0; j < depth; j++) {
		Vectord theta(tree.nLeaves() / 3);
		// for controlled Givens, need one angle per three qubits
		for (size_t i = 0; i < theta.dim(); ++i) {
			theta[i] = gauss(gen);
			cout << theta[i] << " ";
		}
		theta.print();
		params.emplace_back(theta);
	}
	return params;
}

void vqe(TensorTreecd &psi, const string &filename, const Tree &tree, const size_t depth,
         size_t max_iter) {
    size_t nqubits = tree.nLeaves();
    mt19937 gen;
    vector<Vectord> params = initializeVQEparameters(gen, depth, tree);
    Fidelity f;

    int num_params = size(params);

    Vectord gradient(num_params);
    for (int k=0; k<num_params; ++k){
        gradient[k] = 1.;
    }

	SOPcd H = electronicStructure("hamiltonian.dat");

    double Hexp = expectationVal(H, psi, tree);
    cout << Hexp << endl;

	double eps = 1e-3;
	double eps_convergence = 1e-7;
	auto psi0 = psi;
    for (int i = 0; i < max_iter || gradient.norm()<eps_convergence; i++) {
        Circuit circ = paramCirc(nqubits, depth, params, 0);
        Circuit circ_adj = paramCirc(nqubits, depth, params, 1);
        MatrixTreecd rho = TreeFunctions::contraction(psi, tree, true);

        // act on zero state with circ
        psi = psi0;
        TreeFunctions::applyOperator(psi, rho,
                                     circ, circ_adj,
                                     tree, f, gen);

        Vectord flattened = flatten(params, num_params, depth);
        //double energy = getH(flattened);
        double energy = expectationVal(H, psi, tree);
        cout <<"\n" << i << " " << energy << "\n";
        Vectord epsilons(num_params);
        for(int i=0; i<num_params; ++i){
            epsilons[i] = eps;
        }
        gradient = getGradient(flattened, epsilons, energy);
        flattened = updateParameters(flattened, gradient, 0.1);
        // step size 0.1
        params = reshape(flattened, depth, tree.nLeaves() / 3);
    }
}