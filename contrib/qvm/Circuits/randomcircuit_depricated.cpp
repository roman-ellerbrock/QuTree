//
// Created by Roman Ellerbrock on 2/17/21.
//

namespace Circuit {
	SOPVectorcd random1D(const Register& reg, size_t depth, bool adjungate) {
		mt19937 gen(2020);
		uniform_int_distribution<size_t> distribution(0, 1);
		SOPVectorcd stack;

		cout << "=============================\n";
		cout << "1D Random circuit with depth " << depth << endl;
		cout << "=============================\n";
		cout << "Format:\n";
		cout << "X => sqrtX\n";
		cout << "Y => sqrtY\n";
		cout << "T => T\n";
		cout << "C-Z => controlled-Z gate\n";
		cout << "adjoint: " << adjungate << endl;

		// Apply Hadamard gates to all qubits
		stack.append(HadamardChain(reg));

		// Create a history vector for the random gates applied (rules on p. 596 of Characterizing quantum supremacy)
		//    0 = no gate applied yet, 1 = X^1/2 applied, 2 = Y^1/2 applied, 3 = T applied
		std::vector<size_t> history(reg.size(), 0);
		std::vector<bool> ready(reg.size(), true);
		//size_t map[4] = {3, 1, 4, 2};
		size_t map[4] = {2, 0, 3, 1};

		// Apply cycles of CZ and random gates until depth d is reached
		for (size_t d = 0; d < depth; d++) {
			size_t offset = map[d % 4];
			SOPVectorcd layer;
			MLOcd random_gates;
			string operations(reg.size(), '-');
			for (size_t i = 0; i < reg.size(); i++) {
				size_t target = reg.front() + i;
				size_t mod_i = i % 4;
				if (mod_i == offset) {    // Handle CZ gates here
					if (i != reg.size() - 1) {
						layer.append(makeCGate(Z, target, target + 1));
						// set ready to true for next layer
						ready[i] = true;
						ready[i + 1] = true;
						operations[target] = 'C';
						operations[target + 1] = 'Z';
					} else {
						// set ready to false for next layer
						ready[i] = false;
					}
				} else if ((mod_i != offset and mod_i != (offset + 1) % 4) or i == 0) {   // Handle single qubit gates
					if (ready[i]) {
						Matrixcd Op = identityMatrixcd(2);
						switch (history[i]) {
							case 0: {
								Op = T();
								operations[target] = 'T';
								history[i] = 3;
								break;
							}
							case 1: {
								size_t rand = distribution(gen); // append either Y^1/2 or T gate
								if (rand) {
									Op = sqrtY();
									operations[target] = 'Y';
									history[i] = 2;
								} else {
									Op = T();
									operations[target] = 'T';
									history[i] = 3;
								}
								break;
							}
							case 2: {
								size_t rand = distribution(gen); // append either X^1/2 or T gate
								if (rand) {
									Op = sqrtX();
									operations[target] = 'X';
									history[i] = 1;
								} else {
									Op = T();
									operations[target] = 'T';
									history[i] = 3;
								}
								break;
							}
							case 3: {
								size_t rand = distribution(gen); // append either X^1/2 or Y^1/2 gate
								if (rand) {
									Op = sqrtX();
									operations[target] = 'X';
									history[i] = 1;
								} else {
									Op = sqrtY();
									operations[target] = 'Y';
									history[i] = 2;
								}
								break;
							}
						}
						random_gates.push_back(Op, target, adjungate); // append T gate
					}
					// set ready to false for next layer
					ready[i] = false;
				}
			}
			print(operations, d);
			layer.append(random_gates);
			stack.append(layer);
		}

		return stack;
	}
}