//
// Created by Roman Ellerbrock on 4/15/20.
//
#ifndef QVM_H
#define QVM_H
#include "Instruction.h"
#include "Applications/TreeApplyOperatorDynamic.h"
#include "QuantumCircuit.h"
#include "yaml-cpp/yaml.h"
#include "Circuits/Register.h"
#include "QVMState.h"
#include "Measurements.h"


class QVM {
public:
	QVM() = default;
	explicit QVM(const YAML::Node& node);
	~QVM() = default;

	void run(const YAML::Node& node);

	QVMState state_;
	map<string, QuantumCircuit> circuit_;
	Measurements::Sample sample_;
private:

	void setTree(const YAML::Node& node, const Register& reg);
	void setTree(Tree& tree, const YAML::Node& node, const Register& reg);
	void setWavefunction(const YAML::Node& node);
	void output(const YAML::Node& node, ostream& os);
	void crossEntropyBenchmark(const YAML::Node& job);
	void write(const YAML::Node& job);
    void runvqe(const YAML::Node& job);

	void sample(const YAML::Node& job);
	void runFullRank(const YAML::Node& job);
	void writeFidelity(const YAML::Node& job);

	void runUtility(const YAML::Node& job, const Tree& tree, const Tree& fr_tree);

	void overlapFidelity(); /// quick hack
};


#endif //QVM_H
