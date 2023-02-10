//
// Created by Roman Ellerbrock on 4/21/20.
//

#include "QuantumCircuit.h"
#include "QuantumInstruction.h"
#include "Label.h"
#include "GoTo.h"
#include "IfElse.h"
#include "MeasurementInstruction.h"
#include "Measurements.h"
#include "ConditionalQuantumInstruction.h"
#include "InTimeQuantumInstruction.h"
#include "Circuits/QuantumCircuits.h"
#include "OutputInstruction.h"
#include "Circuits/Random.h"

QuantumCircuit::QuantumCircuit(const YAML::Node& node, const Register& reg) {
	if (auto par = node["name"]) {
		name_ = par.as<string>();
	} else {
		cerr << "Did not specify quantum program name\n";
		exit(1);
	}

	if (node["instructions"] && node["file"]) {
		cerr << "Select instructions or file, not both.\n";
		exit(1);
	}

	if (!node["instructions"] && !node["file"]) {
		if (node["type"]) {
			auto type = eval<string>(node, "type");
			if (type == "Shor" || type == "Shors") {
				auto a = (cpp_int) eval<size_t>(node, "a");
				auto N = (cpp_int) eval<size_t>(node, "N");
				size_t nbit = reg.size();
//				*this = Circuits::Shor(reg, long_integer(a, nbit), long_integer(N, nbit));
				*this = Circuits::Shor(reg, a, N);
            } else if (type == "ShorFull") {
                auto a = (cpp_int) eval<size_t>(node, "a");
                auto N = (cpp_int) eval<size_t>(node, "N");
                *this = Circuits::ShorFull(reg, a, N);
			} else if (type == "random1D") {
				auto p = (double) eval<double>(node, "p");
				auto depth = (size_t) eval<size_t>(node, "depth");
				auto seed = (size_t) eval<size_t>(node, "seed");
				*this = Circuits::random1D(reg, depth, p, seed);
			} else if (type == "ising10") {
				auto p = (double) eval<double>(node, "p");
				auto steps = (double) eval<size_t>(node, "steps");
				*this = Circuits::isingModel(p, steps);
			} else if (type == "isingModel") {
				auto p = (double) eval<double>(node, "p");
				auto depth = (double) eval<double>(node, "depth");
				auto steps = (double) eval<size_t>(node, "steps");
				*this = Circuits::isingIntegration(reg, depth, p, steps);

			} else if (type == "random2D") {
				auto p = (double) eval<double>(node, "p");
				auto depth = (size_t) eval<size_t>(node, "depth");
				mt19937 gen(time(nullptr));
				*this = Circuits::random2D(reg, depth, gen, p);
			} else {
				cerr << "Predefined quantum circuit name unknown.\n";
				exit(3);
			}
		} else {
			cerr << "Select instructions or file, not both.\n";
			exit(1);
		}
	}

	if (node["instructions"]) {
		for (const YAML::Node& instr : node["instructions"]) {
//			if (auto instruction_type = instr["type"]) {
			if (auto instruction_type = instr["instr"]) {
				auto type = instruction_type.as<string>();
				if (type == "label") {
					emplace_back(make_shared<Label>(instr));
				} else if (type == "measurement") {
					emplace_back(make_shared<MeasurementInstruction>(instr, reg));
				} else if (type == "goto") {
					emplace_back(make_shared<GoTo>(instr));
				} else if (type == "ifelse") {
					emplace_back(make_shared<IfElse>(instr));
				} else if (type == "quantum") {
					emplace_back(make_shared<QuantumInstruction>(instr, reg));
				} else if (type == "conditionalquantum") {
					emplace_back(make_shared<ConditionalQuantumInstruction>(instr, reg));
				} else if (type == "intimequantum") {
					emplace_back(make_shared<InTimeQuantumInstruction>(instr, reg));
				} else if (type == "output") {
					emplace_back(make_shared<OutputInstruction>());
				} else {
					cerr << "Error: Instruction type unknown.\n";
					exit(1);
				}
			} else {
				cerr << "Did not specify instruction type for:\n";
				cerr << instr << endl;
				exit(1);
			}
		}
	} else if (node["file"]) {
		cerr << "File selection not implemented, yet.\n";
		exit(1);
	}
}

void setLabels(QVMState& state, const QuantumCircuit& prog) {
	/// Gather labels
	state.label_.clear();
	size_t pos = 0.;
	for (auto& instruction : prog) {
		instruction->setLabels(pos, state.label_);
	}

	if (!state.label_.empty()) {
		cout << "Labels:\n";
		for (const auto& x : state.label_) {
			cout << "name: " << x.first << ", pos: " << x.second << endl;
		}
		cout << endl;
	}
}

void run(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, const QuantumCircuit& prog) {
/// Run
/// Loop until at end instruction
	cout << "Quantum simulation:" << endl;
	setLabels(state, prog);
	state.pos_ = 0;
	while (state.pos_ != prog.size()) {
		auto stack_pos = state.pos_; /// Position in stack
		auto& next_instruction = prog[stack_pos]; /// get next instruction
//		cout << "Next instruction: " << stack_pos << endl;
		next_instruction->apply(Psi, rho, state, cout);
	}
	cout << "Finished quantum simulation.\n";
}

void QuantumCircuit::execute(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state) const {
	run(Psi, rho, state, *this);
}

void run(FullRank::Wavefunction& Psi, QVMState& state, const QuantumCircuit& prog) {
	/// Loop until at end instruction
	cout << "Full-rank Quantum simulation:" << endl;
	setLabels(state, prog);
	state.pos_ = 0;
	while (state.pos_ != prog.size()) {
		auto stack_pos = state.pos_; /// Position in stack
		auto& next_instruction = prog[stack_pos]; /// get next instruction
//		cout << "Next instruction: " << stack_pos << endl;
		next_instruction->apply(Psi, state, cout);
	}
}

void QuantumCircuit::execute(FullRank::Wavefunction& Psi, QVMState& state) const {
	run(Psi, state, *this);
}
