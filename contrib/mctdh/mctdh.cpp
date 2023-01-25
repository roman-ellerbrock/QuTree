//
// Created by Roman Ellerbrock on 2/27/20.
//

#include <string>
#include "mctdh.h"
#include "TreeShape/TreeFactory.h"
#include "Hamiltonian_parser.h"
#include "Util/Overlaps.h"
#include "Util/normal_modes.h"
#include "Applications/Eigenstates.h"
#include "Applications/SCF.h"
#include <algorithm>
#include "yaml_extension.h"

namespace parser {

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
		return tree;}

	Tree read_tree(const YAML::Node& node) {
		auto type = evaluate<string>(node, "type");
		if (type == "balanced") {
			auto num_leaves = evaluate<size_t>(node, "number_leaves");
			auto dim_leaves = evaluate<size_t>(node, "dimension_leaves", 2);
			auto dim_nodes = evaluate<size_t>(node, "dimension_nodes", 2);
			auto dim_inc = evaluate<size_t>(node, "dimension_increment", 0);
			auto leaf_type = evaluate<size_t>(node, "leaf_type", 0);
			auto shuffle = evaluate<size_t>(node, "shuffle_indices", 0);
			auto integers = evaluate<size_t>(node, "integers", 0);
			auto bits = evaluate<size_t>(node, "bits", 0);
			///
			auto omega = evaluate<double>(node, "omega", 1.);
			auto r0 = evaluate<double>(node, "r0", 0.);
			auto wfr0 = evaluate<double>(node, "wfr0", 0.);
			auto wfomega = evaluate<double>(node, "wfomoega", 1.);
			Tree tree = TreeFactory::balancedTree(num_leaves, dim_leaves,
				dim_nodes, dim_inc, leaf_type,
				omega, r0, wfr0, wfomega);
			auto checkGrid = evaluate<double>(node, "checkGrid", 0);
			if (checkGrid) {
				cout << "checking grid:\n";
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
			if (integers && bits) {
				std::map<size_t, size_t> idxs = TreeFactory::leaves_staggered_integers(integers, bits);
				tree.reindexLeafModes(idxs);
/*				for (const Node& node : tree) {
					if (node.isBottomlayer()) {
						const Leaf& leaf = node.getLeaf();
						node.info();
						cout << "mode = " << leaf.mode() << endl;
					}
				}
				getchar();*/
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
		par.conversion = evaluate<double>(node, "conversion", 219474.6313705e0);
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
		state.rng_ = mt19937(time(NULL));

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
			} else if (job == "seed") {
				auto seed = evaluate<size_t>(node, "seed");
				state.rng_ = mt19937(seed);
			}
		}
		return state;
	}
}


