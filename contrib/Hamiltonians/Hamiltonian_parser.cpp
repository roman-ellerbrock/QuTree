//
// Created by Roman Ellerbrock on 12/31/22.
//

#include "Hamiltonian_parser.h"
#include "yaml_extension.h"

shared_ptr<Hamiltonian> read_hamiltonian(const YAML::Node& node,
	const Tree& tree) {

	auto name = evaluate<string>(node, "name");
	shared_ptr<Hamiltonian> H_ptr(new Hamiltonian);
	Hamiltonian& H = *H_ptr;

	if (name == "coupled_ho") {
		H = CoupledHO(tree);
	} else if (name == "kinetic_energy") {
		H = Operator::KineticEnergy(tree);
	} else if (name == "nocl") {
		bool V = evaluate<bool>(node, "V", bool(true));
		H = Operator::NOCl_H(V);
	} else if (name == "ch3_meanfield") {
		H = Operator::CH3_meanfield();
	} else if (name == "exciton") {
		H = Operator::Exciton("matrix.out", tree);
	} else if (name == "ndi") {
		H = Operator::ndi();
	} else if (name == "ch3_quasiexact") {
		CH3_quasiexact Hch3(tree);
		H = Hch3;
	} else if (name == "electronicstructure") {
		auto hfile = evaluate<string>(node, string("hamiltonian"));
		H = electronicStructure(hfile);
	} else if (name == "fermigas") {
		H = fermiGas2D(2);
	}else if (name == "schaepers") {
		// find the masses supplied for this hamiltonian
		auto masses = evaluate<string>(node, string("masses"));
		stringstream masses_ss(masses);
		vector<double> massvec;
		while(masses_ss.good()){
			string substr;
			getline(masses_ss, substr, ',');
			massvec.push_back(stod(substr));
		}

		// find the coupling construction for this hamiltonian
		auto coupling = evaluate<string>(node, string("coupling"));
		stringstream coupling_ss(coupling);
		vector<int> couplingvec;
		if(massvec.size() > 2){
			while(coupling_ss.good()){
				string substr;
				getline(coupling_ss, substr, ',');
				couplingvec.push_back(stoi(substr));
			}
		}

		// init schaepers vector
		H = Operator::schaepers(tree, massvec, couplingvec);

	}else if (name == "portfoliooptimization") {
		auto Na = evaluate<size_t>(node, "nAssets", 25);
		auto NaTot = evaluate<size_t>(node, "nAssetsTotal", 99);
		auto Nt = evaluate<size_t>(node, "nTime", 1);
		auto NtTot = evaluate<size_t>(node, "nTimeTotal", 180);
		auto Nq = evaluate<size_t>(node, "nQubitsPerAsset", 1);
		auto alpha = evaluate<double>(node, "alpha", 1.);
		auto gamma = evaluate<double>(node, "gamma", 1.);
		auto rho = evaluate<double>(node, "rho", 1.);
		auto K = evaluate<double>(node, "investment", 1.);
		auto tickers = evaluate<string>(node, "tickers", "merged.csv");
		H = meanVarianceAnalysis(tickers, Na, Nt, NaTot, NtTot, Nq,
			alpha, gamma, rho, K);
	} else if (name == "portfoliooptimization2") {
		auto Na = evaluate<size_t>(node, "nAssets", 25);
		auto NaTot = evaluate<size_t>(node, "nAssetsTotal", 99);
		auto Nt = evaluate<size_t>(node, "nTime", 1);
		auto NtTot = evaluate<size_t>(node, "nTimeTotal", 180);
		auto Nq = evaluate<size_t>(node, "nQubitsPerAsset", 1);
		auto alpha = evaluate<double>(node, "alpha", 1.);
		auto gamma = evaluate<double>(node, "gamma", 1.);
		auto rho = evaluate<double>(node, "rho", 1.);
		auto K = evaluate<double>(node, "investment", 1.);
		auto tickers = evaluate<string>(node, "tickers", "merged.csv");
		meanVarianceAnalysisOptimization(tickers, Na, Nt, NaTot, NtTot, Nq,
			alpha, gamma, rho, K);
	} else if (name == "ising2D") {
		auto Lx = evaluate<size_t>(node, "Lx", 32);
		auto Ly = evaluate<size_t>(node, "Ly", 5);
		auto h = evaluate<double>(node, "h", 2.9);
		cout << Lx << " x " << Ly << endl;
		cout << "h = " << h << endl;
		H = ising2D(Lx, Ly, h);

	} else {
		cout << "No valid Hamiltonian name." << endl;
		cout << "Chosen name: " << name << endl;
		cout << "See Parser/yaml_parser.cpp function read_hamiltonian for a list of choices.\n" << endl;
		exit(1);
	}
	assert(H_ptr->size() > 0);
	return H_ptr;
}

PotentialOperator set_potential(const YAML::Node& node, const Tree& tree) {
	auto name = evaluate<string>(node, "name");
	if (name == "coupled_ho") {
		auto coupling = evaluate<bool>(node, "coupling", "true");
		auto V = make_shared<CDVRModelV>(tree.nLeaves(), coupling);
		return PotentialOperator(V, 0, 0);
	} else if (name == "nocl") {
		auto state = evaluate<string>(node, "state", "S1");
		bool S1 = (state == "S1");
		auto V = make_shared<NOClPotential>(S1);
		return PotentialOperator(V, 0, 0);
	} else if(name == "ch3") {
		auto V = make_shared<CH3Potential>();
		Vectord mass(4);
		mass(0) = 12.0;
		mass(1) = 1.007;
		mass(2) = 1.007;
		mass(3) = 1.007;
		PotentialOperator Vop(V, 0, 0);
		Vop.Q_ = make_shared<TrafoCH3Quasiexact>(mass);
		return Vop;
	} else if(name == "ch3_normal_modes") {
		auto V = make_shared<CH3Potential>();

		PotentialOperator Vop(V, 0, 0);
		Vop.Q_ = make_shared<NormalModes>("U.dat", "x0.dat", 6);
		return Vop;
	} else if(name == "liuch4cl") {
		// find the masses supplied for this potential
		auto masses = evaluate<string>(node, "masses");
		stringstream masses_ss(masses);
		vector<double> massvec;
		while(masses_ss.good()){
			string substr;
			getline(masses_ss, substr, ',');
			massvec.push_back(stod(substr));
		}
		// find the coupling construction for this potential
		auto coupling = evaluate<string>(node, "coupling");
		stringstream coupling_ss(coupling);
		vector<int> couplingvec;
		if(massvec.size() > 2){
			while(coupling_ss.good()){
				string substr;
				getline(coupling_ss, substr, ',');
				couplingvec.push_back(stoi(substr));
			}
		}
		auto V = make_shared<liuch4cl>(massvec,couplingvec);
		PotentialOperator Vop(V, 0, 0);
		/*
		Vectord mass(4);
		mass(0) = 12.0;
		mass(1) = 1.007;
		mass(2) = 1.007;
		mass(3) = 1.007;
		Vop.Q_ = std::make_shared<TrafoCH3Quasiexact>(mass);
		 */
		return Vop;

	} else {
		cerr << "Did not recognise potential energy operator name\n";
		exit(1);
	}
	return PotentialOperator();
}

