//
// Created by Roman Ellerbrock on 2/19/20.
//
#include "TreeClasses/TensorTree.h"
#include "TreeShape/TreeFactory.h"
#include "Core/stdafx.h"
#include <random>
#include "Circuits/Register.h"
#include "Circuits/QFT.h"
#include "TreeClasses/TreeIO.h"
#include "Measurements.h"
#include "FullRank.h"
#include "QVM.h"

void FullRankSimulation();
void TTSimulation();
void runQVM(const string& filename);

int main(int argc, char* argv[]) {

	if (argc != 2){
		cerr << "Please provide a name of a .yaml-file as argument." << endl;
		exit(1);
	}
	string filename(argv[1]);

	runQVM(filename);
}

void runQVM(const string& filename) {
	auto node = YAML::LoadFile(filename); // e.g. ../examples/example2.yaml
	QVM qvm(node);
}

void TTSimulation() {
	mt19937 gen(1993);
	Tree tree = TreeFactory::balancedTree(8, 2,  2);
	TensorTreecd Psi(gen, tree, true);

	Register reg(0, 8, "main");
	auto Hs = Circuits::HadamardChain(reg);
	auto QFT = Circuits::QFT(reg, false, 48);
	auto QFT_T = Circuits::QFT(reg, true, 48);

	MatrixTreecd rho = TreeFunctions::contraction(Psi, tree, true);
//	TreeFunctions::ApplyOperator(Psi, rho, Hs, Hs, tree);
//	TreeFunctions::ApplyOperator(Psi, rho, QFT, QFT_T, tree);

	TreeIO::output(Psi, tree);

	// @TODO: Check changing dimensions. / doesn't seem to work, yet.

	vector<size_t> targets;
	for (size_t k = 0; k < tree.nLeaves(); ++k) {
		targets.emplace_back(k);
	}

	auto M = Measurements::sample(Psi, gen, 200, targets, tree);
	Measurements::print(M);
}

void FullRankSimulation() {

	mt19937 gen(1993);
	Tree tree = TreeFactory::balancedTree(12, 2,  2);

	Register reg(0, 12, "main");
	auto Hs = Circuits::HadamardChain(reg);
	auto QFT = Circuits::QFT(reg, false, 48);
	auto QFT_T = Circuits::QFT(reg, true, 48);

	FullRank::Wavefunction Psi = FullRank::initialize(tree);
	Psi(0) = 1.;
	Psi = FullRank::applyOperator(Psi, Hs, tree);
	Psi = FullRank::applyOperator(Psi, QFT, tree);
	FullRank::print(Psi);
}

/**
 * =============================
 * Cross-entropy benchmarking
 * =============================
 *
 * 1.) Run full-rank simulation (save wavefunction)
 * 2.) Run tensortree simulation
 * 3.) Sample ttwf
 * 4.) Evaluate correct probablities for measurements
 * 5.) evaluate cross-entropy
 *
 */




