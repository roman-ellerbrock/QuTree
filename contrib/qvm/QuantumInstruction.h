//
// Created by Roman Ellerbrock on 4/15/20.
//

#ifndef QUANTUMINSTRUCTION_H
#define QUANTUMINSTRUCTION_H
#include "TreeOperators/SOPVector.h"
#include "Applications/TreeApplyOperatorDynamic.h"
#include "Instruction.h"
#include "Circuits/Register.h"

vector<Register> TargetRegisters(const YAML::Node& node, const Register& reg);

class QuantumInstruction: public Instruction {
public:
	QuantumInstruction() = default;

	QuantumInstruction(const YAML::Node& node, const Register& reg);

	QuantumInstruction(const YAML::Node& node, const QVMState& state);

	~QuantumInstruction() = default;

	/// Quasi-move constructor
	QuantumInstruction(SOPVectorcd&& stack, SOPVectorcd&& stack_adj)
		: stack_(move(stack)), stack_adj_(move(stack_adj)) {}

	/// Quasi-Copy constructor
	QuantumInstruction(SOPVectorcd stack, SOPVectorcd stack_adj)
		: stack_(move(stack)), stack_adj_(move(stack_adj)) {}

	/// Generate Quantum Instruction from a operator-name
	QuantumInstruction(const string& circuit_name, const vector<Register>& regs);

	///
	void apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) override;

	void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) override;

private:
	SOPVectorcd stack_;
	SOPVectorcd stack_adj_;
};

void SelectQuantumInstruction(SOPVectorcd& stack,
	const string& circuit_name, const vector<Register>& targets, bool adjoint,
	const YAML::Node* node = nullptr);

void InTimeSelectQuantumInstruction(SOPVectorcd& stack,
	const string& circuit_name, const vector<Register>& targets,
	const QVMState& state, bool adjoint, const YAML::Node *node);

template <typename F, typename... Args>
auto translate(const F &f, Args &&...args) -> std::enable_if_t<
	std::is_convertible_v<decltype(f(std::forward<Args>(args)..., std::declval<bool>())), SOPVectorcd>,
	shared_ptr<QuantumInstruction>>
{
	SOPVectorcd Q = f(std::forward<Args>(args)..., false);
	SOPVectorcd Qadj = f(std::forward<Args>(args)..., true);
	return make_shared<QuantumInstruction>(Q, Qadj);
}


#endif //QUANTUMINSTRUCTION_H
