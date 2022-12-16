//
// Created by Roman Ellerbrock on 2/27/20.
//

#include <string>
#include "mctdh.h"
#include "TreeShape/TreeFactory.h"
#include "Hamiltonians.h"
#include "Util/Overlaps.h"
#include "Util/normal_modes.h"
#include "Applications/Eigenstates.h"
#include "Applications/SCF.h"
#include <algorithm>

namespace parser {

	template <typename T>
	T evaluate(const YAML::Node& node, const string& key) {
		T val;
		if (auto par = node[key]) {
			return par.as<T>();
		} else {
			cerr << "Did not specify key '" << key << "' in yaml node " << node << endl;
			exit(3);
		}
	}

	template <typename T>
	T evaluate(const YAML::Node& node, const string& key, T default_) {
		T val;
		if (auto par = node[key]) {
			return par.as<T>();
		} else {
			return default_;
		}
	}

	Leaf create_leaf(const YAML::Node& yaml_node) {
		auto dim = yaml_node["dim"].as<size_t>();
		auto mode = yaml_node["mode"].as<size_t>();
		auto type = yaml_node["leaftype"].as<size_t>();
		PhysPar par;
		return Leaf(dim, mode, type, 0, par);
	}

	Node create_node(const YAML::Node& yaml_node) {
		Node node;
		auto dim = yaml_node["dim"].as<size_t>();
		vector<size_t> dims;
		for (const auto& child : yaml_node["children"]) {
			auto type = child["type"].as<string>();
			dims.push_back(child["dim"].as<size_t>());
			if (type == "leaf") {
				Leaf leaf = create_leaf(child);
				node = Node(leaf, dim);
			} else if (type == "node") {
				node.push_back(create_node(child));
			} else {
				cerr << "Invalid node type." << endl;
				cerr << "Choices: (leaf, node)" << endl;
				exit(1);
			}
		}

		dims.push_back(dim);
		node.shape() = TensorShape(dims);
		return node;
	}

	void read_leaf_parameters(Tree& tree, const YAML::Node& node) {
		for (const auto& child : node["leaves"]) {
			auto mode = evaluate<size_t>(child, "mode");
			auto& leaf = tree.getLeaf(mode);
			auto r0 = evaluate<double>(child, "r0");
			auto wfr0 = evaluate<double>(child, "wfr0");
			auto omega = evaluate<double>(child, "omega");
			auto wfomega = evaluate<double>(child, "wfomega");
			auto& grid = leaf.interface();
			grid.initialize(omega, r0, wfr0, wfomega);
		}
	}

	Tree create_tree(const YAML::Node& node) {
		Tree tree;
		Node root = create_node(node["tree"]);
		tree.setRoot(root);
		tree.update();
		read_leaf_parameters(tree, node);
		return tree;
	}

	Tree read_tree(const YAML::Node& node) {
		auto type = evaluate<string>(node, "type");
		if (type == "balanced") {
			auto num_leaves = evaluate<size_t>(node, "number_leaves");
			auto dim_leaves = evaluate<size_t>(node, "dimension_leaves", 2);
			auto dim_nodes = evaluate<size_t>(node, "dimension_nodes", 2);
			auto dim_inc = evaluate<size_t>(node, "dimension_increment", 0);
			auto leaf_type = evaluate<size_t>(node, "leaf_type", 0);
			auto shuffle = evaluate<size_t>(node, "shuffle_indices", 0);
			///
			auto omega = evaluate<double>(node, "omega", 1.);
			auto r0 = evaluate<double>(node, "r0", 0.);
			auto wfr0 = evaluate<double>(node, "wfr0", 0.);
			auto wfomega = evaluate<double>(node, "wfomoega", 1.);
			Tree tree = TreeFactory::balancedTree(num_leaves, dim_leaves,
				dim_nodes, dim_inc, leaf_type,
				omega, r0, wfr0, wfomega);
			auto checkGrid = evaluate<double>(node, "checkGrid", 1.);
			if (checkGrid) {
				for (const Node& node: tree) {
					if (node.isBottomlayer()) {
						const Leaf& leaf = node.getLeaf();
						leaf.interface().getX().print();
					}
				}
			}
			if (shuffle) {
				/// shuffle leaf modes
				std::vector<size_t> v(tree.nLeaves());
				std::generate(v.begin(), v.end(), [n = 0]() mutable { return n++; });
				std::shuffle(v.begin(), v.end(), mt19937(time(nullptr)+23498));
				map<size_t, size_t> m;
				size_t i = 0;
				for (auto x : v) {
					m[i] = x;
					i++;
				}
				tree.reindexLeafModes(m);
			}
			return tree;
		} else if (type == "manual") {
			return create_tree(node);
		} else if (type == "compact") {
			auto tree_str = evaluate<string>(node, "tree");
			stringstream ss(tree_str);
			Tree tree(ss);
			if (!tree.isWorking()) {
				cerr << "Failed to read tree with .yaml parser in compact format.\n";
				cerr << "Error caused by tree input string reading:\n";
				cerr << tree_str << endl;
				exit(2);
			}
			return tree;
		} else {
			cerr << "No valid tree type." << endl;
			cerr << "Choices: (balanced, manual)" << endl;
			exit(1);
		}
		return Tree();
	}

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
			bool V = evaluate<bool>(node, "V", true);
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
			auto hfile = evaluate<string>(node, "hamiltonian");
			H = electronicStructure(hfile);
		} else if (name == "fermigas") {
			H = fermiGas2D(2);
		}else if (name == "schaepers") {
			// find the masses supplied for this hamiltonian
			auto masses = evaluate<string>(node, "masses");
			stringstream masses_ss(masses);
			vector<double> massvec;
			while(masses_ss.good()){
				string substr;
				getline(masses_ss, substr, ',');
				massvec.push_back(stod(substr));
			}

			// find the coupling construction for this hamiltonian
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
			auto rho  = evaluate<double>(node, "rho", 1.);
			auto K  = evaluate<double>(node, "investment", 1.);
			auto tickers = evaluate<string>(node, "tickers", "merged.csv");
			H = meanVarianceAnalysis(tickers, Na, Nt, NaTot, NtTot, Nq,
				alpha, gamma, rho, K);
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

	void new_wavefunction(mctdh_state& state, const YAML::Node& node) {
		auto name = evaluate<string>(node, "name");
		auto type = evaluate<string>(node, "type");
		if (type == "read") {
            if(evaluate<string>(node, "filename").empty()){
                cerr << "supply filename to save and read directive" << endl;
                exit(1);
            }
            auto filename = evaluate<string>(node, "filename");
			Wavefunction Psi(state.tree_);
			ifstream is(filename);
			while(is.peek() != EOF) { /// read last wavefunction
				Psi.read(is);
			}
			state.wavefunctions_[name] = Psi;
			is.close();
		} else if (type == "create") {
			bool Hartree = evaluate<bool>(node, "Hartree", true);
			Wavefunction Psi(state.rng_, state.tree_, Hartree);
//			occupyCIS(Psi, state.tree_);
			state.wavefunctions_[name] = Psi;
		} else if (type == "save") {
		    if(evaluate<string>(node, "filename").empty()){
		        cerr << "supply filename to save and read directive" << endl;
		        exit(1);
		    }
		    auto filename = evaluate<string>(node, "filename");
		    Wavefunction Psi(state.tree_);
		    Psi = state.wavefunctions_[name];
		    ofstream of(filename);
		    Psi.write(of);
		    of.close();
		} else {
			cerr << "No valid Wavefunction initialization type." << endl;
			cerr << "Choices: (read, create, save)" << endl;
			exit(1);
		}
	}

	IntegratorVariables new_ivar(const YAML::Node& node, mctdh_state& state) {
	    // TODO: this routine does not use the 'save'-directive
	    // TODO: also, it does not save the wavefunction and only works with a wavefunction called "Psi"
		auto t_end = evaluate<double>(node, "t_end", 100*41.362);
		auto t = evaluate<double>(node, "t", 0.);
		auto out = evaluate<double>(node, "out", 41.362);
		auto dt = evaluate<double>(node, "dt", 1.);
		auto cmf = evaluate<double>(node, "cmf", 1e-4);
		auto bs = evaluate<double>(node, "bs", 1e-5);
		auto file_in = evaluate<string>(node, "file_in", "in.dat");
		auto file_out = evaluate<string>(node, "file_out", "out.dat");
		auto append = evaluate<bool>(node, "append", true);
		IntegratorVariables ivar(t, t_end, dt, out, cmf, bs,
			state.wavefunctions_["Psi"], *state.hamiltonian_,
			state.tree_, state.cdvrtree_, file_out, file_in, append);
		return ivar;
	}

	SCF_parameters scf_parameters(const YAML::Node& node, mctdh_state& state) {
		SCF_parameters par;
		par.nIter = evaluate<size_t>(node, "nIter", 20);
		par.nKrylov = evaluate<size_t>(node, "nKrylov", 20);
		par.nITP = evaluate<size_t>(node, "nITP", 0);
		par.beta = evaluate<double>(node, "beta", 1);
		par.output = evaluate<bool>(node, "output", 1);
		par.psi = &state.wavefunctions_["Psi"];
		par.h = state.hamiltonian_.get();
		par.tree = &state.tree_;
		return par;
	}

	mctdh_state read(const string& yaml_filename) {
		YAML::Node config = YAML::LoadFile(yaml_filename);
		mctdh_state job;
		job.tree_ = read_tree(config["tree"]);
		job.hamiltonian_ = read_hamiltonian(config["hamiltonian"], job.tree_);
		new_wavefunction(job, config["wavefunction"]);
		return job;
	}

	mctdh_state run(const string& yaml_filename) {
		YAML::Node config = YAML::LoadFile(yaml_filename);
		mctdh_state state;

		for (const auto& node : config["run"]) {
			const auto& job = node["job"].as<string>();
			if (job == "tree") {
				state.tree_ = read_tree(node);
				state.cdvrtree_ = state.tree_;
//				if (state.cdvrtree_.nNodes() == 0) { state.cdvrtree_ = state.tree_; }

			} else if (job == "hamiltonian") {
				state.hamiltonian_ = read_hamiltonian(node, state.tree_);
			} else if (job == "potential") {
				PotentialOperator V = set_potential(node, state.tree_);
				state.hamiltonian_->V_ = V;
				state.hamiltonian_->hasV = true;
			} else if (job == "wavefunction") {
				new_wavefunction(state, node);
			} else if (job == "eigenstates") {
				auto ivar = new_ivar(node, state);
				Eigenstates(ivar);
			} else if (job == "scf") {
				auto par = scf_parameters(node, state);
				scf(par);
			} else if (job == "cmf") {
				auto ivar = new_ivar(node, state);
				const Hamiltonian& H = *ivar.h;
				const Tree& tree = *ivar.tree;
				const Tree& cdvrtree = *ivar.cdvrtree;
				CMFIntegrator cmf(H, tree, cdvrtree, 1.);
				cmf.Integrate(ivar);
			} else if (job == "cdvrtree") {
				state.cdvrtree_ = read_tree(node);
			} else if (job == "overlaps") {
				auto file1 = node["bra"].as<string>();
				auto file2 = node["ket"].as<string>();
				wavefunctionOverlap(file1, file2, state.tree_);
			} else if (job == "find_minimum") {
				auto ivar = new_ivar(node, state);
				const auto H = state.hamiltonian_;
				const Tree& tree = state.tree_;
				find_minimum(*H, tree);
			}
		}
		return state;
	}
}


