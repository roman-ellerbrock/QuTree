//
// Created by Roman Ellerbrock on 4/23/20.
//

#include "MeasurementInstruction.h"
#include "Measurements.h"
#include "QuantumInstruction.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"

MeasurementInstruction::MeasurementInstruction(const vector<size_t>& targets,
	const string& creg, bool append)
	: targets_(targets), classical_reg_name_(creg), append_(append) {
}

void MeasurementInstruction::apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) {
	++state.pos_;
	if (targets_.empty()) { return; }

	const Tree& tree = state.tree_;
	auto Psilast = Psi;
	auto measurement = Measurements::measurement(Psi, state.gen_, targets_, tree);
	SparseTree stree(targets_, tree);
	TreeFunctions::contraction(rho, Psi, stree, true);

	/// Does measurement count as decreasing the fidelity?
	bool measure_fidelity = false;
	if (measure_fidelity) {
		auto S = TreeFunctions::dotProduct(Psi, Psilast, state.tree_);
		double f = real(S[tree.topNode()](0, 0));
		state.f_.append(f);
	}

	qrOrthogonal(Psi, tree);
	orthonormal(Psi, tree);

	/// either push back measurement or (over-)write
	if (state.classical_reg_.count(classical_reg_name_) && append_) {
		vector<size_t>& m = state.classical_reg_[classical_reg_name_];
		for (auto x : measurement) {
			m.push_back(x);
		}
	} else {
		state.classical_reg_[classical_reg_name_] = measurement;
	}
}

void MeasurementInstruction::apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) {
	++state.pos_;
	cerr << "Measurement instruction not implemented for full-rank wavefunction.\n";
	exit(1);
}

MeasurementInstruction::MeasurementInstruction(const Register& qreg,
	const string& creg, bool append) : append_(append)  {
	type_ = "measurement";

	targets_.clear();
	for (size_t i = 0; i < qreg.size(); ++i) {
		size_t idx = qreg.front() + i;
		targets_.emplace_back(idx);
	}
	if (targets_.empty()) {
		cerr << "Error: Measurement Instruction on empty register!\n";
		exit(1);
	}

	classical_reg_name_ = creg;
}

MeasurementInstruction::MeasurementInstruction(const YAML::Node &node,
	const Register& reg) {

	type_ = "measurement";
	if (auto type_par = node["instr"]) {
		assert(type_par.as<string>() == type_);
	} else {
		cerr << "Type flag missing\n";
		exit(2);
	}

	append_ = false;
	if (node["append"]) {
		append_ = eval<bool>(node, "append");
	}

	if (auto quantum_reg_node = node["quantum_reg"]) {
		const auto& quantum_regs = TargetRegisters(quantum_reg_node, reg);
//		const Register& quantum_reg = reg[quantum_reg_node];
		targets_.clear();
		for (const Register& quantum_reg : quantum_regs) {
			for (size_t i = 0; i < quantum_reg.size(); ++i) {
				size_t idx = quantum_reg.front() + i;
				targets_.emplace_back(idx);
			}
		}
	} else {
		cerr << "Specify quantum_reg for MeasurementInstruction\n";
		exit(2);
	}

	if (auto classical_reg_par = node["classical_reg"]) {
		classical_reg_name_ = classical_reg_par.as<string>();
	} else {
		cerr << "Specify classical_reg for MeasurementInstruction\n";
		exit(2);
	}

}

