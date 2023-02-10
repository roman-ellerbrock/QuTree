//
// Created by Grace Johnson on 4/21/20.
//

#include "Random.h"
#include "../FullRank.h"
#include "TreeShape/TreeFactory.h"
#include <vector>
#include <random>
#include "GateOperators.h"
#include "../RandomProjectorInstruction.h"

namespace Circuits {

	void print(const string& operations, size_t d) {
		cout << "d = " << d;
		if (d < 10) {
			cout << "  | ";
		} else {
			cout << " | ";
		}
		for (size_t j = 0; j < operations.size(); ++j) {
			if (operations[j] == 'C') {
				cout << operations[j] << "-" << operations[j + 1] << " ";
				++j;
			} else {
				cout << operations[j] << " ";
			}
		}
		cout << endl;
	}

	vector<vector<size_t>> red() {
		vector<vector<size_t>> circuit;
		circuit.push_back({1, 6});
		circuit.push_back({2, 7});
		circuit.push_back({3, 8});
		circuit.push_back({4, 9});

		circuit.push_back({10, 15});
		circuit.push_back({11, 16});
		circuit.push_back({12, 17});
		circuit.push_back({13, 18});

		circuit.push_back({19, 24});
		circuit.push_back({20, 25});
		circuit.push_back({21, 26});
		circuit.push_back({22, 27});

		circuit.push_back({28, 33});
		circuit.push_back({29, 34});
		circuit.push_back({30, 35});
		circuit.push_back({31, 36});

		circuit.push_back({37, 42});
		circuit.push_back({38, 43});
		circuit.push_back({39, 44});
		circuit.push_back({40, 45});

		circuit.push_back({46, 51});
		circuit.push_back({47, 52});
		circuit.push_back({48, 53});
		circuit.push_back({49, 54});

		return circuit;
	}

	vector<vector<size_t>> dark_green() {
		vector<vector<size_t>> circuit;
		circuit.push_back({2, 6});
		circuit.push_back({3, 7});
		circuit.push_back({4, 8});
		circuit.push_back({5, 9});

		circuit.push_back({11, 15});
		circuit.push_back({12, 16});
		circuit.push_back({13, 17});
		circuit.push_back({14, 18});

		circuit.push_back({20, 24});
		circuit.push_back({21, 25});
		circuit.push_back({22, 26});
		circuit.push_back({23, 27});

		circuit.push_back({29, 33});
		circuit.push_back({30, 34});
		circuit.push_back({31, 35});
		circuit.push_back({32, 36});

		circuit.push_back({38, 42});
		circuit.push_back({39, 43});
		circuit.push_back({40, 44});
		circuit.push_back({41, 45});

		circuit.push_back({47, 51});
		circuit.push_back({48, 52});
		circuit.push_back({49, 53});
		circuit.push_back({50, 54});

		return circuit;
	}

	vector<vector<size_t>> pink() {
		vector<vector<size_t>> circuit;
		circuit.push_back({6, 10});
		circuit.push_back({7, 11});
		circuit.push_back({8, 12});
		circuit.push_back({9, 13});

		circuit.push_back({15, 19});
		circuit.push_back({16, 20});
		circuit.push_back({17, 21});
		circuit.push_back({18, 22});

		circuit.push_back({24, 28});
		circuit.push_back({25, 29});
		circuit.push_back({26, 30});
		circuit.push_back({27, 31});

		circuit.push_back({33, 37});
		circuit.push_back({34, 38});
		circuit.push_back({35, 39});
		circuit.push_back({36, 40});

		circuit.push_back({42, 46});
		circuit.push_back({43, 47});
		circuit.push_back({44, 48});
		circuit.push_back({45, 49});

		return circuit;
	}

	vector<vector<size_t>> light_green() {
		vector<vector<size_t>> circuit;

		circuit.push_back({6, 11});
		circuit.push_back({7, 12});
		circuit.push_back({8, 13});
		circuit.push_back({9, 14});

		circuit.push_back({15, 20});
		circuit.push_back({16, 21});
		circuit.push_back({17, 22});
		circuit.push_back({18, 23});

		circuit.push_back({24, 29});
		circuit.push_back({25, 30});
		circuit.push_back({26, 31});
		circuit.push_back({27, 32});

		circuit.push_back({33, 38});
		circuit.push_back({34, 39});
		circuit.push_back({35, 40});
		circuit.push_back({36, 41});

		circuit.push_back({42, 47});
		circuit.push_back({43, 48});
		circuit.push_back({44, 49});
		circuit.push_back({45, 50});

		return circuit;
	}

	SOPVectorcd Translate2Gate(const vector<vector<size_t>>& targets, bool adj) {
		SOPVectorcd circuit;
		for (const auto& twoqubits : targets) {
			if (twoqubits.size() != 2) {
				cerr << "2 Qubit targets required.\n";
				exit(1);
			}
			size_t q1 = twoqubits[0] - 1;
			size_t q2 = twoqubits[1] - 1;
			if (q1 > 53) {
				cerr << "wrong qubit index on q1!\n";
				exit(1);
			}
			if (q2 > 53) {
				cerr << "wrong qubit index on q2!\n";
				exit(1);
			}

//			circuit.push_back(Circuits::CNot(q1, q2));
			auto Z = make_shared<LeafMatrixcd>(sigma_z(), false);
			circuit.push_back(Circuits::makeCGate(Z, q1, q2));
//			circuit.append(Circuits::iSWAPSOP(q1, q2, QM::pi / 6., adj));
		}

/*		SOPVectorcd circuit_contr; /// Do multiple operations in a single operator
		size_t n = 2;
		for (size_t i = 0; i < circuit.size(); i += n) {
			SOPcd x(circuit[i]);
			for (size_t j = 1; j < n; ++j) {
				if (i + j >= circuit.size()) { continue; }
				x = circuit[i + j] * x;
			}
			circuit_contr.append(x);
			cout << "x size: " << x.size() << endl;
		}
		circuit = circuit_contr;
		cout << "circuit size: " << circuit.size() << endl;
		*/
		return circuit;
	}

	SOPVectorcd ColoredGrid(const string& color, bool adj) {
		vector<vector<size_t>> targets;
		cout << "color: " << color << endl;
		if (color == "red") {
			targets = red();
		} else if (color == "darkgreen") {
			targets = dark_green();
		} else if (color == "pink") {
			targets = pink();
		} else if (color == "lightgreen") {
			targets = light_green();
		} else {
			cout << "Color unknown. Choose from red, darkgreen, lightgreen and pink.\n";
			exit(1);
		}
		return Translate2Gate(targets, adj);
	}

	Matrixcd SigmaMix(double alpha, double phi) {
		Matrixcd s = sin(alpha) * sin(phi) * sigma_x();
		s += sin(alpha) * cos(phi) * sigma_y();
		s += cos(alpha) * sigma_z();
		return s;
	}

	Matrixcd RandomRot(mt19937& gen, bool adjoint) {
		uniform_real_distribution<double> dist(0., QM::two_pi);
		double alpha = dist(gen);
		double phi = dist(gen);
		double theta = dist(gen) / 2.;
		auto A = SigmaMix(alpha, phi);
		auto x = diagonalize(A);
		const Matrixcd& U = x.first;
		Vectord& ew1 = x.second;
		Vectorcd ew(2);
		ew(0) = exp(QM::im * theta * ew1(0));
		ew(1) = exp(QM::im * theta * ew1(1));
		Matrixcd B(2, 2);
		for (size_t i = 0; i < B.dim1(); ++i) {
			for (size_t j = 0; j < B.dim2(); ++j) {
				for (size_t k = 0; k < B.dim2(); ++k) {
					B(i, j) += U(i, k) * ew(k) * conj(U(j, k));
				}
			}
		}
		if (adjoint) { B = B.adjoint(); }
		return B;
	}

	MLOcd Random1DGates(const Register& reg, mt19937& gen, bool adjoint) {
		MLOcd M;
		for (size_t i = 0; i < reg.size(); ++i) {
			M.push_back(RandomRot(gen, adjoint), i);
		}
		return M;
//		SOPcd S(M);
//		return S;
	}

	SOPVectorcd Random2D(const Register& reg, size_t depth, mt19937& gen, bool adjoint) {
		SOPVectorcd circuit;
		SOPVectorcd circuit_adj;
		vector<string> colors = {"red", "lightgreen", "pink", "darkgreen", "pink", "darkgreen", "red", "lightgreen"};
		size_t size = colors.size();
		for (size_t d = 0; d < depth; ++d) {
			string color = colors[d % size];
			cout << color << endl;
			circuit.append(Random1DGates(reg, gen, adjoint));
			circuit.append(ColoredGrid(color, adjoint));
		}
		return circuit;
	}

	QuantumCircuit random2D(const Register& reg, size_t depth, mt19937& gen, double p) {
		QuantumCircuit circuit;
		vector<string> colors = {"red", "lightgreen", "pink", "darkgreen", "pink", "darkgreen", "red", "lightgreen"};
		size_t size = colors.size();
		for (size_t d = 0; d < depth; ++d) {
			string color = colors[d % size];
			cout << color << endl;
			circuit.emplace_back(translate(Random1DGates, reg, gen));
			circuit.emplace_back(translate(ColoredGrid, color));
//			if (color == "pink" || color == "lightgreen" ) {
			circuit.emplace_back(randomMeasurement(gen, reg, p));
//			}
		}
		return circuit;
	}

	SOPVectorcd random2D(const Register& reg, size_t depth, size_t seed, bool adjoint) {
		mt19937 gen(seed);
		return Random2D(reg, depth, gen, adjoint);
	}

	MLOcd RandomClifford(mt19937& gen, size_t q) {
		uniform_real_distribution<double> dist(0., 1.);
		double r = dist(gen);
		if (r < 0.5) {
			return MLOcd(Hadamard, q);
		} else {
			return MLOcd(S, q);
		}
	}

	SOPVectorcd RandomClifford(mt19937& gen, size_t q1, size_t q2) {
		uniform_real_distribution<double> dist(0., 1.);
		double r = dist(gen);
		if (r < 0.25E0) {
			return SOPVectorcd(CNot(q1, q2));
		} else if (r < 0.75E0) {
			MLOcd M = RandomClifford(gen, q1);
			M = M * RandomClifford(gen, q2);
			return SOPVectorcd(M);
		} else {
			return SWAP(q1, q2);
		}
	}

	SOPVectorcd Random1DSlice(mt19937& gen, const Register& reg, size_t off) {
		SOPVectorcd circ;
		for (size_t i = off; i < reg.size(); i += 2) {
			size_t q1 = i;
			size_t q2 = (i + 1) % (reg.size() - 1);
			circ.append(RandomClifford(gen, q1, q2));
		}
		return circ;
	}

	SOPVectorcd CNot1D(const Register& reg, size_t off, bool adjoint) {
		SOPVectorcd circ;
		for (size_t i = 0; i < reg.size() - 1; i += 2) {
			size_t q1 = (i + off);
			size_t q2 = (i + off + 1) % (reg.size());
			circ.append(CNot(q1, q2));
		}
		return circ;
	}

	shared_ptr<MeasurementInstruction> randomMeasurement(mt19937& gen,
		const Register& reg, double p) {
		vector<size_t> targets;
		for (size_t i = 0; i < reg.size(); ++i) {
			uniform_real_distribution<double> dist(0., 1.);
			double r = dist(gen);
			if (r < p) { targets.push_back(i); }
		}
		string creg("creg");
		return make_shared<MeasurementInstruction>(targets, creg, true);
	}

	QuantumCircuit random1D(const Register& reg,
		size_t depth, double p, size_t seed) {
		QuantumCircuit circ;
		mt19937 gen(seed);
		for (size_t d = 0; d < depth; ++d) {
			size_t off = d % 2;
			circ.emplace_back(translate(Random1DGates, reg, gen));
			circ.emplace_back(translate(CNot1D, reg, off));
			if (p > 0) { circ.emplace_back(make_shared<RandomProjectorInstruction>()); }
			if (d == (depth / 2)) { circ.emplace_back(make_shared<RandomProjectorInstruction>()); }
//			if (p > 0) { circ.emplace_back(randomMeasurement(gen, reg, p)); }
		}
//		@TODO: Add Instruction for random natural orbital projection
//		circ.emplace_back(make_shared<RandomProjectorInstruction>());
		return circ;
	}

}

