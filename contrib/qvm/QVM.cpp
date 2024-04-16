//
// Created by Roman Ellerbrock on 4/15/20.
//

#include "QVM.h"
#include "TreeShape/TreeFactory.h"
#include "TreeClasses/TreeIO.h"
#include "FullRank.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "Util/OverlapUtilities.h"
#include "Util/filenames.h"
#include "Circuits/VQE.h"
#include "Circuits/ShorSparsetensor.h"
#include "NumberTheory.h"

QVM::QVM(const YAML::Node& node) {

	if (node["Register"]) {
		state_.reg_ = Register(node["Register"], 0);
		cout << "Register:\n";
		state_.reg_.print();
		cout << endl;
	} else {
		cerr << "No registers specified.\n";
		exit(1);
	}

	if (!node["Circuit"]) {
		cerr << "No circuit specified.\n";
		exit(1);
	}
	cout << "Circuits:" << endl;
	auto circ = node["Circuit"];
	auto name = eval<string>(circ, "name");
	circuit_[name] = QuantumCircuit(circ, state_.reg_);
	cout << "size of circuit: " << circuit_.at(name).size() << endl;
	cout << endl;

	/// Execute job script
	cout << "QVM:" << endl;
	if (const auto& jobs = node["QVM"]) {
		for (const YAML::Node& job: jobs) {
			if (!job["job"]) {
				cerr << "specify job.\n";
				exit(1);
			}
			const auto type = job["job"].as<string>();
			if (type == "tree") {
				setTree(job, state_.reg_);
			} else if (type == "wavefunction") {
				setWavefunction(job);
			} else if (type == "run") {
				run(job);
			} else if (type == "output") {
				output(job, cout);
			} else if (type == "sample") {
				sample(job);
			} else if (type == "runfullrank") {
				runFullRank(job);
			} else if (type == "crossentropy") {
				crossEntropyBenchmark(job);
//			} else if (type == "vqe") {
//				runvqe(job);
			} else if (type == "write_fidelity") {
				writeFidelity(job);
			} else if (type == "write") {
				write(job);
			} else if (type == "seed") {
				state_.gen_.seed(eval<size_t>(job, "seed"));
			} else if (type == "utility") {
				runUtility(job, state_.tree_, state_.fr_tree_);
			} else {
				cerr << "Unknown command.\n";
				exit(1);
			}
		}
	} else {
		cerr << "No jobs specified.\n";
		exit(2);
	}
}

void QVM::sample(const YAML::Node& job) {
	auto name = eval<string>(job, "name");
	auto nsample = eval<size_t>(job, "number_samples");
	const Register& reg = state_.reg_;
	vector<size_t> targets;
	for (size_t i = 0; i < reg.size(); ++i) {
		size_t idx = reg.front() + i;
		targets.push_back(idx);
	}
	sample_ = Measurements::sample(state_.wavefunction_[name], state_.gen_,
		nsample, targets, state_.tree_);
	for (const auto& x: sample_) {
		for (auto y: x.first) {
			cout << y << " ";
		}
		cout << " | " << x.second;
		cout << endl;
//		/// print out measurement in integer format
//		size_t val = 0;
//		for (size_t i = 0; i < reg["M"].size(); ++i) {
//            /// big endian, read largest first
//            size_t i_inv = (reg["M"].size()-1) - i;
//		    val += x.first[i] * pow(2, i_inv);
//		}
//        size_t x_val = 0;
//        for (size_t i = 0; i < reg["x"].size(); ++i) {
//            /// big endian, read largest first
//            size_t i_inv = (reg["x"].size()-1) - i;
//            x_val += x.first[i + reg["M"].size()] * pow(2, i_inv);
//        }
//
//        cout << "\t" << val << "\t" << x_val << endl;
	}
}

void QVM::output(const YAML::Node& node, ostream& os) {
	auto name = eval<string>(node, "name");

	TreeIO::output(state_.wavefunction_[name], state_.tree_, os);

	cout << "Classical register(s):\n";
	for (const auto& x: state_.classical_reg_) {
		cout << x.first << " | ";
		long_integer m(0, 2 * x.second.size());
		size_t k = 0;
		for (auto bit: x.second) {
			cout << bit << " ";
			m.set_bit(k++, bit);
		}

		auto number = m.convert();
		cout << " | " << number << endl;
	}
}

void QVM::runFullRank(const YAML::Node& job) {
	/// Get circuit
	cout << "- run full-rank:" << endl;
	string prog_name = eval<string>(job, "circuit");
	const QuantumCircuit& prog = circuit_.at(prog_name);

	/// initialize wavefunction
	FullRank::Wavefunction& Psi = state_.frpsi_;
	Psi = FullRank::initialize(state_.tree_);

	prog.execute(Psi, state_);

//	modify(Psi);
//	Psi_approx = modify(Psi);

//	prog2.execute(Psi, state_);
//	prog2.execute(Psi_approx, state_);
}

void QVM::run(const YAML::Node& node) {
	/// Get Circuit
	cout << "- run:" << endl;
	string prog_name = eval<string>(node, "circuit");
	cout << "Circuits:\n";
	for (const auto& x: circuit_) {
		cout << x.first << ",\t number of instructions: " << x.second.size() << endl;
	}
	cout << "Run circuit: " << prog_name << endl << endl;
	QuantumCircuit& prog = circuit_.at(prog_name);

	/// Get Wavefunction
	string wf_name;
	if (auto wf_par = node["wavefunction"]) {
		wf_name = wf_par.as<string>();
	} else {
		cerr << "Specify a wavefunction name.\n";
		exit(2);
	}
	Wavefunction& Psi = state_.wavefunction_.at(wf_name);
	MatrixTreecd& rho = state_.rho_.at(wf_name);

	/// Initialize timer
	auto t1 = chrono::high_resolution_clock::now();

	//run(state_.wavefunction_.at(wf_name), state_, prog);
	prog.execute(Psi, rho, state_);

	auto t2 = chrono::high_resolution_clock::now();;
	microseconds global_time = chrono::duration_cast<chrono::microseconds>(t2 - t1);
	cout << "Total time after running Quantum Circuit: " << to_string(global_time.count() / 1000000.) << "s\n";
}

void QVM::setTree(const YAML::Node& node, const Register& reg) {
	string treename = "tree";
	if (node["name"]) {
		treename = eval<string>(node, "name");
	}
	if (treename == "tree") {
		setTree(state_.tree_, node, reg);
	} else if (treename == "fr") {
		setTree(state_.fr_tree_, node, reg);
	} else {
		cerr << "Error: wrong tree name [tree, fr].\n";
		exit(2);
	}
}

void QVM::setTree(Tree& tree, const YAML::Node& node, const Register& reg) {
	size_t n_qubit = reg.size();
	if (auto nqubits_par = node["number_qubits"]) {
		n_qubit = nqubits_par.as<size_t>();
		assert(n_qubit >= reg.size());
	}

	string treetype = "binary";
	if (auto type_par = node["type"]) {
		treetype = type_par.as<string>();
	} else {
		cerr << "Error: tree type not set.\n";
		exit(3);
	}

	size_t dim = 2;
	if (auto dim_par = node["dim"]) {
		dim = dim_par.as<size_t>();
		assert(dim > 0);
	}

	if (treetype == "binary") {
		cout << "- Generate binary tree with " << n_qubit << " qubits.\n";
		tree = TreeFactory::balancedTree(n_qubit, 2, dim, 0, 6);
	} else if (treetype == "train") {
		tree = TreeFactory::unbalancedTree(n_qubit, 2, dim, 6);
	} else if (treetype == "compact") {
		auto tree_str = eval<string>(node, "tree");
		stringstream ss(tree_str);
		tree = Tree(ss);
	} else {
		cerr << "Tree type unknown. Chosen type: " << treetype << endl;
		exit(3);
	}

	string indextype = "default";
	if (auto index_par = node["indexing"]) {
		indextype = index_par.as<string>();
		cout << "Leaf indexing: " << indextype << endl;
	}
	if (indextype == "staggered") {
		map<size_t, size_t> Map;
		size_t half = n_qubit / 2;
		for (size_t i = 0; i < half; ++i) {
			Map[2 * i] = i;
			Map[2 * i + 1] = i + half;
		}
		assert((n_qubit % 2) == 0); /// @TODO: CHECK EDGE CASE
		if (n_qubit % 2) {
			Map[2 * half] = 2 * half;
		}
		tree.reindexLeafModes(Map);
	}
}

void QVM::setWavefunction(const YAML::Node& node) {
	string name = "Psi";
	if (auto name_par = node["name"]) {
		name = name_par.as<string>();
		cout << "- Generate wavefunction named: " << name << endl;
	} else {
		cout << "- Generate wavefunction with default name: " << name << endl;
	}

	if (node["from_file"]) {
		auto filename_par = eval<string>(node, "from_file");
		//filename = filename_par.as<string>();
		cout << "- Loading wavefunction from file: " << filename_par << endl;
		Wavefunction Psi(filename_par);
		//canonicalTransformation(Psi, state_.tree_, true);
		state_.wavefunction_[name] = Psi;
		//state_.rho_[name] = TreeFunctions::contraction(Psi, state_.tree_, true);
	} else {
		Wavefunction Psi(state_.gen_, state_.tree_);
		canonicalTransformation(Psi, state_.tree_, true);
		state_.wavefunction_[name] = Psi;
		state_.rho_[name] = TreeFunctions::contraction(Psi, state_.tree_, true);
	}
	const Tree& tree = state_.tree_;
	size_t dim = 0;
	for (const Node& node: tree) {
		dim += node.shape().totalDimension();
	}
	cout << "Number of wavefunction parameters: " << dim << endl;
}

void QVM::crossEntropyBenchmark(const YAML::Node& job) {
	cout << "Evaluating cross-entropy of sample and full-rank wavefunction:\n";
	bool verbose = false;
	if (job["verbose"]) {
		verbose = job["verbose"].as<bool>();
	}
	auto H = Measurements::crossEntropy(sample_, state_.frpsi_, verbose);
	auto S = Measurements::entropy(state_.frpsi_);
	auto Spt = Measurements::sPorterThomas(state_.tree_.nLeaves());
	cout << "H = " << H << "\tS = " << S << endl;
	cout << "Delta = " << S - H << endl;
	cout << "sPorterThomas = " << Spt << endl;
	cout << "alpha = " << (S - H) / S << endl;
}

/*void QVM::runvqe(const YAML::Node& job) {
	auto filename = eval<string>(job, "filename");
	size_t depth = eval<size_t>(job, "depth");
	size_t max_iter = eval<size_t>(job, "max_iter");
	string wavefunction = eval<string>(job, "wavefunction");
	TensorTreecd& psi = state_.wavefunction_[wavefunction];
	const Tree& tree = state_.tree_;
	vqe(psi, filename, tree, depth, max_iter);
}*/

void QVM::writeFidelity(const YAML::Node& job) {
	string name = filename("fidelity", ".dat");
	string header = "# f - f_avg - logf_total";
	if (job["filename"]) {
		name = eval<string>(job, "filename");
		name = filename(name, ".dat");
	}
	if (job["header"]) { header = eval<string>(job, "header"); }
	state_.f_.write(name, header);
}

void QVM::write(const YAML::Node& job) {
	string psiname = eval<string>(job, "wavefunction");
	const TensorTreecd& Psi = state_.wavefunction_.at(psiname);

	string name = eval<string>(job, "name");
	name = filename(name, ".dat");
	Psi.write(name);
}

void QVM::overlapFidelity() {
	vector<string> names({"psi/b.l6-2.dat", "psi/b.l6-3.dat", "psi/b.l6-4.dat", "psi/b.l6-5.dat",
						  "psi/b.l6-6.dat", "psi/b.l6-8.dat", "psi/b.l6-16.dat"});
	vector<TensorTreecd> Psis;
	for (const auto& x: names) {
		Psis.emplace_back(TensorTreecd(x));
	}
	string name = "psi/b.l6-16.dat";
	TensorTreecd Psi(name);
	vector<TensorTreecd> xs({Psi});
	Utility::wavefunctionOverlap(state_.tree_, xs, Psis);
}

#include <chrono>

void QVM::runUtility(const YAML::Node& job, const Tree& tree, const Tree& fr_tree) {

	auto type = eval<string>(job, "type");
	if (type == "stat") {
		Utility::statisticalWavefunctionOverlap(tree);
	} else if (type == "xeb") {
		string psi_fr = eval<string>(job, "file_fr");
		string psi_xeb = eval<string>(job, "file_xeb");
		Utility::xeb(fr_tree, tree, psi_fr, psi_xeb);
	} else if (type == "xeb_stat") {
		auto psiname = eval<string>(job, "file_bra");
		auto statname = eval<string>(job, "file_ket");
		auto xebname = eval<string>(job, "file_xeb");
		auto nsample = eval<size_t>(job, "nsample");
		xebname = filename(xebname, ".dat");
		Wavefunction Psi(psiname);
		statistical::Wavefunctions statPsis(statname, ".dat");
		ofstream os(xebname);
		Utility::xeb_stat(os, state_.gen_, Psi, statPsis, state_.tree_, nsample);
	} else if (type == "overlap") {
		overlapFidelity();
	} else if (type == "shors_random") {
		auto n = eval<size_t>(job, "n");
		auto nsample = eval<size_t>(job, "nsample");
		mt19937 gen(time(nullptr));
		runShorOnCoPrimes(n, gen, nsample);
	} else if (type == "shors_primes") {
		auto n = eval<size_t>(job, "n");
		auto m = eval<size_t>(job, "m");
		auto nsample = eval<size_t>(job, "nsample");
		mt19937 gen(time(nullptr));
		runShorOnCoPrimes(n, m, gen, nsample);
	} else if (type == "shors") {
		auto a = eval<size_t>(job, "a");
		auto N = eval<size_t>(job, "N");
		auto n = eval<size_t>(job, "n");
		auto nsample = eval<size_t>(job, "nsample");
		cout << "Run Sparsetensor version of Shor's...\n";
		cout << "a = " << a << endl;
		cout << "N = " << N << endl;
		cout << "n = " << n << endl;
		cout << "nsample = " << nsample << endl;
		mt19937 gen(time(nullptr));

		ShorSparsetensor state(a, N, n, gen);
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		for (size_t i = 0; i < nsample; ++i) {
			cout << i << "\t";
			cout.flush();
			state.sample();
			state.evaluateStatistics();
			state.printResult();
			state.clear();
		}
		end = std::chrono::system_clock::now();
		auto time = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000000.;
		auto timePerBitstring = time / (double) nsample;
		auto bitstringsPerSecond = 1. / timePerBitstring;
		cout << "s(r) = " << state.s_r_ << ", s(factor(r)) = " << state.s_factor_ << ", s_total = " << state.s_tot_
			 << endl;
		cout << "nSample: " << nsample << ", Time: " << time << "s, time/bitstring: " << timePerBitstring
			 << "s, bitstrings/s: " << bitstringsPerSecond << endl;
	} else if (type == "screen_a") {
		auto N = eval<size_t>(job, "N");
		auto nsample = eval<size_t>(job, "nsample");
		cout << "Run Sparsetensor version of Shor's...\n";
		cout << "N = " << N << endl;
		cout << "M = " << nsample << endl;
		mt19937 gen(time(nullptr));

		auto f = number_theory::factorize(N);
		cout << "factorization of N = " << N << "\n";
		cout << "prod_k p_k^e_k\n";
		cout << "p_k,\te_k\n";
		for (auto p: f) {
			cout << p.first << "\t" << p.second << endl;
		}
		cout << "Square free (N) = " << number_theory::is_square_free(f);

//		auto r = screen_a(N, gen, nsample);
/*		for (size_t i = 0; i < r.size(); ++i) {
			cout << r[i] << "\t" << N << endl;
		}*/

	} else if (type == "modpow") {
		cpp_int a = eval<size_t>(job, "a");
		cpp_int N = eval<size_t>(job, "N");
		cpp_int r = eval<size_t>(job, "r");
		cout << "a = " << a << endl;
		cout << "N = " << N << endl;
		cout << "r = " << r << endl;
		cpp_int aa = 1;
		for (size_t i = 0; i < r / 2; ++i) {
			aa *= a;
			aa %= N;
		}
		cout << "a^(r/2) mod N = " << aa << endl;
		cout << "gcd(x^(r/2) + 1) = " << gcd(aa + 1, N) << endl;
		cout << "gcd(x^(r/2) - 1) = " << gcd(aa - 1, N) << endl;
		aa = 1;
		for (size_t i = 0; i < r; ++i) {
			aa *= a;
			aa %= N;
		}
		cout << "a^r mod N = " << aa << endl;
	} else if (type == "factorize") {
		using namespace number_theory;
		auto N = eval<size_t>(job, "N");
		auto a = eval<size_t>(job, "a");

		auto f = number_theory::factorize(N);

		cout << "factorization of N = " << N << "\n";
		bool beg = true;
		for (auto p: f) {
			if (!beg) {
				cout << " x ";
			}
			beg = false;
			cout << "U(" << p.first << "^" << p.second << ")";
		}
		cout << " ~= ";
		beg = true;
		for (auto p: f) {
			if (!beg) {
				cout << " x ";
			}
			beg = false;
			cout << "C(" << p.first - 1 << ")";
		}
		cout << endl;

		cout << "Periods factorization:\n";
		vector<PrimeFactorization> rfs = factorizePeriods(f);

		cout << "gcd(r_1, r_2, ...) = " << gcd(rfs) << endl;

		/// Get Z_N^ *
//		vector<cpp_int> Z = number_theory::coprimes(N);

		/// phi(N) = |Z_N^*|
//		cpp_int phi = Z.size();
//		cout << "phi(N) = " << phi << endl;

		/// lambda(N)
		auto lambda = lcm(rfs);
		cout << "lambda(N) = " << product(lambda) << "\t(" << lambda << ")" << endl;
//		cout << "zeta(N) = phi(N) / lambda(N) = " << phi / product(lambda) << endl;

		/// get number of elements with lambda(N)
		cout << "order - number of elements with that order\n";
		cpp_int sum = 0;
		auto D = number_of_divisors(lambda);
		map<cpp_int, cpp_int> order_count;
		for (size_t I = 0; I < D; ++I) {
			auto d = divisor(lambda, I);
//			if ((double) (N / product(d)) > 0.0001 * N) { continue; }
			auto res = number_of_elements_with_order(rfs, d);
			sum += res;
			order_count[product(d)] = res;
		}
		double psum = 0.;
		double pfail = 0.;
		cpp_int phi;
		for (const auto& Pair: order_count) {
			phi += Pair.second;
		}
		for (const auto& Pair: order_count) {
			double p = (double) Pair.second / (double) phi;
			if (Pair.first % 2 == 0) {
				psum += p;
				cout << Pair.first << "\t" << Pair.second << "\t" << p << "\t" << psum << endl;
			} else {
//				phi += Pair.second;
				pfail += p;
			}
		}
		cout << "p_fail = " << pfail << endl;
		cout << "sum of all orders: " << sum << ", phi = " << phi << endl;

//		auto lambda = number_theory::lcm(rfs);

//		cout << "factorize with period finding:\n";
//		Shor::factorize(a, N);
		for (size_t i = 9000; i < 20000; ++i) {
			if ((i - 1000) % 20 == 0) { cout << endl; }
			cout << boost::math::prime(i) << "   ";
		}
		cout << endl;
	}

}

