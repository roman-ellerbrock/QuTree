//
// Created by Roman Ellerbrock on 8/28/20.
//

#include "OutputInstruction.h"
#include "TreeClasses/TreeIO.h"
#include "Util/long_integer.h"

void OutputInstruction::apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) {

//	Orthogonal(Psi, state.tree_);
//	Orthonormal(Psi, state.tree_);

//	cout << "Number of entries:" << state.wavefunction_.size() << endl;
	TreeIO::output(Psi, state.tree_, os);

	/// Perform sampling of wavefunction
/*	const Register& reg = state.reg_["x"];
	vector<size_t> targets;
	for (size_t i = 0; i < reg.size(); ++i) {
		size_t idx = reg.front() + i;
		targets.push_back(idx);
	}
	auto sample = Measurements::sample(Psi, state.gen_, 50, targets, state.tree_);
	for (const auto& x : sample) {
		for (auto y : x.first) {
			cout << y << " ";
		}
		vector<size_t> xs = x.first;
		long_integer l(0, xs.size());
		cout << " | ";
		for (size_t j = 0; j < l.size(); ++j){
			l.set_bit(j, xs[l.size() - 1 - j]);
		}
		cout << l.convert() << " | (" << x.second << ")" << endl;
	}
*/
	state.pos_++;
//	getchar();
}

void OutputInstruction::apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) {
	state.pos_++;
	cerr << "Output instruction not employed.\n";
	exit(1);
}
