//
// Created by Roman Ellerbrock on 4/15/20.
//

#include "QuantumInstruction.h"
//#include "Core/TreeApplyOperatorDynamic.h"
#include "Applications/ApplyFirstOrder.h"
#include "Applications/TreeApplyOperator.h"
#include "Circuits/GateOperators.h"
#include "Circuits/QuantumCircuits.h"
#include "Circuits/Arithmetic.h"
#include "Circuits/QFT.h"
#include "Circuits/Random.h"
#include "yaml-cpp/yaml.h"
#include "Instruction.h"

Register TargetFromScalar(const YAML::Node& node, const Register& reg) {
	auto node_str = node.as<string>();
	if (std::is_number(node_str)) {
		size_t idx = stoi(node_str);
		assert(idx < reg.size());
		return Register(reg.front() + idx, 1, reg.name() + "[" + to_string(idx) + "]");
	} else {
		return reg[node];
	}
}

Register SingleTargetFromSequence(const YAML::Node& node, const Register& reg) {
	assert(node.IsSequence());
	/// Return in-place subreg? sub[i]
	Register t = reg[node];
	string back;
	for (const auto& x : node) {
		back = x.as<string>();
	}
	if (std::is_number(back)) {
		size_t idx = stoi(back);
		assert(idx < t.size());
		return Register(t.front() + idx, 1, t.name() + "[" + to_string(idx) + "]");
	} else {
		return t;
	}
}

vector<Register> TargetSequence(const YAML::Node& node, const Register& reg) {
	vector<Register> targets;
	assert(node.IsSequence());
	/// Check whether a subnode is sequence of sequence
	bool multiple_targets = false;
	for (const auto& sub : node) {
		if (sub.IsSequence()) { multiple_targets = true; }
	}

	if (multiple_targets) {
		for (const auto& sub : node) {
			assert(sub.IsSequence()); //@TODO: Add non-sequences here?
			targets.emplace_back(SingleTargetFromSequence(sub, reg));
		}
	} else {
		targets.emplace_back(SingleTargetFromSequence(node, reg));
	}
	return targets;
}

vector<Register> TargetRegisters(const YAML::Node& node,
	const Register& reg) {
	if (node.IsMap()) {
		cerr << "Maps as targets not implemented.\n";
		exit(2);
	} else if (node.IsSequence()) {
		return TargetSequence(node, reg);
	} else if (node.IsScalar()) {
		return {TargetFromScalar(node, reg)};
	} else {
		cerr << "Unknown node type for target Register.\n";
		exit(2);
	}
}

void SelectQuantumInstruction(SOPVectorcd& stack,
	const string& circuit_name, const vector<Register>& targets,
	bool adjoint, const YAML::Node *node) {
	assert(!targets.empty());

	stack.clear();

	using namespace Circuits;
	if (circuit_name == "QFT") {
		stack = QFT(targets.front(), adjoint);
	} else if (circuit_name == "iQFT") {
		stack = iQFT(targets.front(), adjoint);
		cout << "iQFT: " << stack.size() << endl;
	} else if (circuit_name == "H") {
		stack = HadamardChain(targets.front());
	} else if (circuit_name == "X") {
		stack.append(distribute(X, targets.front(), adjoint));
	} else if (circuit_name == "CNot") {
		assert(targets.size() == 2);
		stack.append(CNot(targets[0], targets[1]));
	} else if (circuit_name == "CZ") {
		assert(targets.size() == 2);
		stack.append(CZ(targets[0], targets[1]));
	} else if (circuit_name == "Toffoli") {
		assert(targets.size() == 3);
		stack.append(Toffoli(targets[0], targets[1], targets[2]));
	} else if (circuit_name == "T") {
		stack.append(distibute(T(), targets.front(), adjoint));
	} else if (circuit_name == "sqrtX") {
		stack.append(distibute(sqrtX(), targets.front(), adjoint));
	} else if (circuit_name == "sqrtY") {
		stack.append(distibute(sqrtY(), targets.front(), adjoint));
	} else if (circuit_name == "CNotChain") {
		stack.append(CNotChain(targets.front()));
	} else if (circuit_name == "QuantumAddition") {
		assert(targets.size() == 2);
		stack.append(Circuits::add_qft(targets[0], targets[1], adjoint, 48));
	} else if (circuit_name == "SetNumber") {
		assert(targets.size() == 1);
		auto b = eval<size_t>(*node, "value");
		long_integer bb(b, targets[0].size());
		stack.append(Circuits::Set_Number(targets[0], bb));
	} else if (circuit_name == "cUa") {
		assert(targets.size() == 3);
		assert(node != nullptr);
		stack.append(Circuits::cUa(*node, targets[0], targets[1],
			targets[2], targets[3], adjoint, 48));
	} else if (circuit_name == "cMULTmod") {
		assert(targets.size() == 3);
		assert(node != nullptr);
		stack.append(Circuits::cMULTmod(*node, targets[0], targets[1],
			targets[2], targets[3], adjoint));
	} else if (circuit_name == "Shor") {
		if (targets.size() != 4) { cerr << "wrong number of targets.\n"; exit(1); }
		cpp_int a = 7;
		if (auto par = (*node)["a"]) {
			a = par.as<size_t>();
		}
		cpp_int N = 15;
		if (auto par = (*node)["N"]) {
			N = par.as<size_t>();
		}
//		size_t mode_0 = targets[3].SmallestBit();
		size_t approx = 8;
		if (auto par = (*node)["approx"]) {
			approx = par.as<size_t>();
		}
		stack.append(Circuits::Shors<cpp_int>(targets[0], targets[1],
			targets[2], targets[3], a, N, adjoint, approx));
		cout << "Shor's circuit size: " << stack.size() << endl;
	} else if (circuit_name == "ccADDmod") {
		assert(targets.size() == 3);
		stack.append(Circuits::ccADDmod(*node, targets[0], targets[1],
			targets[2], targets[3], adjoint));
	} else if (circuit_name == "ising10") {
		stack.append(Circuits::transversalIsing(100));
	} else if (circuit_name == "random1D") {
		assert(targets.size() == 1);
		size_t depth = 8;
		if (node != nullptr) {
			if ((*node)["depth"]) {
				depth = node->operator[]("depth").as<size_t>();
			}
		}
		double p = 0.;
		size_t seed = 1;
		cerr << "wrong circuit selection.\n";
		getchar();
//		stack.append(random1D(targets.front(), depth, adjoint));
//		stack.append(Circuits::random1D(targets, depth, p, seed));
	} else if (circuit_name == "random2D") {
		assert(targets.size() == 1);
		size_t depth = 8;
		if (auto par = (*node)["depth"]) {
			depth = par.as<size_t>();
		}
		size_t seed = time(nullptr);
		if (auto par = (*node)["seed"]) {
			seed = par.as<size_t>();
		}
		cout << "Random 2D Circuit with depth " << depth << " and seed " << seed << endl;
		stack.append(Circuits::random2D(targets[0], depth, seed, adjoint));
	} else {
		cerr << "Circuit name unknown: " << circuit_name << endl;
		exit(1);
	}
}

QuantumInstruction::QuantumInstruction(const YAML::Node& node,
	const Register& reg) {

	type_ = "QuantumInstruction";

	vector<Register> targets;
	if (auto par = node["reg"]) {
		/// Get target register if "reg" is specified
		targets = TargetRegisters(par, reg);
	} else {
		/// Default: target is everything
		targets.emplace_back(reg);
		cout << "Use default register.\n";
	}

	string circuit_name;
	if (auto par = node["type"]) {
		circuit_name = par.as<string>();
	} else {
		cerr << "Did not specify type of quantum instruction.\n";
	}

	for (const auto& target : targets) {
		cout << "Instruction: " << circuit_name << "[" << target.name() << "]\n";
	}
	SelectQuantumInstruction(stack_, circuit_name, targets, false, &node);
	SelectQuantumInstruction(stack_adj_, circuit_name, targets, true, &node);
}

QuantumInstruction::QuantumInstruction(const string& circuit_name,
	const vector<Register>& regs) {
	SelectQuantumInstruction(stack_, circuit_name, regs, false);
	SelectQuantumInstruction(stack_adj_, circuit_name, regs, true);
}

void QuantumInstruction::apply(Wavefunction& Psi, MatrixTreecd& rho,
	QVMState& state, ostream& os) {
	const string method = "2ndOrder";
//	const string method = "symmetric";
	if (method == "SCF") {
		/// 1st Order method
		TreeFunctions::applyOperatorSCF(Psi, rho, stack_, state.tree_);
	} else if (method == "2ndOrder") {
		/// 2nd Order method - Static bond dimension
		TreeFunctions::applyOperator(Psi, rho, stack_, stack_adj_, state.tree_, state.f_, state.gen_);
	} else if (method == "2ndOrderDynamic") {
		/// 2nd Order method - Dynamic bond dimensio
		TreeFunctions::applyOperator(Psi, state.tree_, stack_, stack_adj_,
			state.par_.eps, state.par_.max_spf, state.par_.plus_spf);
	} else {
		cerr << "Unknown method type.\n";
		exit(1);
	}
	++state.pos_;
}

void QuantumInstruction::apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) {
	Psi = FullRank::applyOperator(Psi, stack_, state.tree_);
	++state.pos_;
}

QuantumInstruction::QuantumInstruction(const YAML::Node& node,
	const QVMState& state) {

	type_ = "QuantumInstruction";

	string circuit_name;
	if (auto par = node["type"]) {
		circuit_name = par.as<string>();
	} else {
		cerr << "Did not specify type of quantum instruction.\n";
	}

	vector<Register> targets;
	if (auto par = node["reg"]) {
		/// Get target register if "reg" is specified
		targets = TargetRegisters(par, state.reg_);
	} else {
		/// Default: target is everything
		targets.emplace_back(state.reg_);
		cout << "Use default register.\n";
	}

	for (const auto& target : targets) {
		cout << "Instruction: " << circuit_name << "[" << target.name() << "]\n";
	}

	InTimeSelectQuantumInstruction(stack_, circuit_name, targets, state, false, &node);
	InTimeSelectQuantumInstruction(stack_adj_, circuit_name, targets, state, true, &node);
}

void InTimeSelectQuantumInstruction(SOPVectorcd& stack,
	const string& circuit_name, const vector<Register>& targets,
	const QVMState& state, bool adjoint, const YAML::Node *node) {
	assert(!targets.empty());

	stack.clear();

	using namespace Circuits;

	if (circuit_name == "Rspecial") {
		stack.append(Rspecial(*node, targets[0], state, adjoint));
	} else {
		cerr << "Unknown circuit name: " << circuit_name << "\n";
		exit(1);
	}

}



