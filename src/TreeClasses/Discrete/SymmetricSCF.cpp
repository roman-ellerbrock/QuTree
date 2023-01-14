//
// Created by Roman Ellerbrock on 12/20/22.
//

#include "TreeClasses/Discrete/SymmetricSCF.h"

template
class Configuration<size_t>;

template<typename T>
ostream& operator<<(ostream& os, const Configuration<T>& c) {
	if (!c.empty()) {
		for (size_t i = 0; i < c.size() - 1; ++i) {
			os << c[i] << " ";
		}
		cout << c.back();
	} else {
		cout << "-";
	}
	return os;
}

template ostream& operator<<(ostream& os, const Configuration<size_t>&);

template<typename T>
ostream& operator<<(ostream& os, const ConfigurationTensor<T>& c) {
	if (!c.empty()) {
		for (size_t i = 0; i < c.size() - 1; ++i) {
			os << c[i] << "\n";
		}
		cout << c.back();
	} else {
		cout << "-";
	}
	return os;
}

template ostream& operator<<(ostream& os, const ConfigurationTensor<size_t>& c);

template<typename T>
ostream& operator<<(ostream& os, const ConfigurationTree<T>& a) {
	cout << "up:\n";
	for (const auto& x: a.up_) {
		cout << x << endl;
	}
	cout << "down:\n";
	for (const auto& x: a.down_) {
		cout << x << endl;
	}
	return os;
}

template ostream& operator<<(ostream& os, const ConfigurationTree<size_t>& c);

ConfigurationTensor<> randomConfigurationTensor(const TensorShape& shape, mt19937& gen) {
	/// Create a range of integers range(0, shape.lastDimension)
	vector<size_t> x(shape.lastBefore());
	iota(x.begin(), x.end(), 0); //0 is the starting number

	/// Shuffle it and take the first shape.lastDimension elements
	shuffle(x.begin(), x.end(), gen);
	vector<size_t> y(x.begin(), x.begin() + shape.lastDimension());
	std::sort(y.begin(), y.end());

	/// convert to ConfigurationTensor
	ConfigurationTensor<> T;
	for (const auto& number: y) {
		Configuration<> c({number});
		T.emplace_back(c);
	}
	return T;
}

ConfigurationTree<> randomConfigurationTree(const Tree& tree, mt19937& gen) {
	ConfigurationTree<> Psi(tree);
	for (const Node& node: tree) {
		if (node.isToplayer()) { continue; }
		ConfigurationTensor<>& A = Psi.up_[node];
		const TensorShape& shape = node.shape();
		if (node.isBottomlayer()) {
			A = randomConfigurationTensor(shape, gen);
		} else {
			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child = node.child(k);
				A = A * Psi.up_[child];
			}
			/// shuffle and truncate
			std::shuffle(A.begin(), A.end(), gen);
			A.resize(shape.lastDimension());
		}
	}

	/// down
	for (int n = tree.nNodes() - 1; n >= 0; --n) {
		const Node& hole_child = tree.getNode(n);
		if (hole_child.isToplayer()) { continue; }

		ConfigurationTensor<>& A = Psi.down_[hole_child];
		const Node& node = hole_child.parent();
		size_t hole = hole_child.childIdx();

		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			if (child.address() == hole_child.address()) { continue; }
			A = A * Psi.up_[child];
		}

		A = A * Psi.down_[node];

		std::shuffle(A.begin(), A.end(), gen);
		A.resize(hole_child.shape().lastDimension());
	}

	return Psi;
}

ConfigurationTensor<> bottomTensor(const TensorShape& shape) {
	vector<size_t> x(shape.lastBefore());
	iota(x.begin(), x.end(), 0); //0 is the starting number
	ConfigurationTensor<> y;
	for (auto z: x) {
		y.push_back({z});
	}
	return y;
}

ConfigurationTensor<> cartesianProduct(Configuration<>& idx, const ConfigurationTree<>& Psi,
	const Node& node) {

	ConfigurationTensor<> A;
	if (node.isBottomlayer()) {
		A = Psi.up_[node];
		idx = Psi.idx_up_[node];
	} else {
		idx.clear();
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			A = A * Psi.up_[child];
			idx = idx * Psi.idx_up_[child];
		}
	}

	if (!node.isToplayer()) {
		A = A * Psi.down_[node];
		idx = idx * Psi.idx_down_[node];
	}
	return A;
}

Configuration<> resort(const Configuration<>& c, const Configuration<>& idx) {
	Configuration<> d(c);
	for (size_t i = 0; i < c.size(); ++i) { d[i] = c[idx[i]]; }
	return d;
}

vector<size_t> findIndices(const vector<size_t>& idx, const vector<size_t>& all) {
	/**
	 * example:
	 * all: {2 3 0 1}
	 * idx: {3 0 1}
	 * return: {1 2 3}
	 */
	vector<size_t> res;
	for (auto i: idx) {
		for (size_t I = 0; I < all.size(); ++I) {
			if (i == all[I]) { res.push_back(I); }
		}
	}
	return res;
}

ConfigurationTensor<> sliceDown(const ConfigurationTensor<>& B, const vector<size_t>& idx) {
	/// get relevant slice
	ConfigurationTensor<> sl;
	for (size_t i = 0; i < B.size(); ++i) {
		const Configuration<>& c = B[i];
		Configuration<> x;
		for (auto j: idx) {
			x.push_back(c[j]);
		}
		sl.push_back(x);
	}
	return sl;
}

ConfigurationTensor<> slice(const ConfigurationTensor<>& B, size_t start, size_t n) {
	/// get relevant slice
	ConfigurationTensor<> sl;
	for (size_t i = 0; i < B.size(); ++i) {
		vector<size_t> x_range(B[i].begin() + start, B[i].begin() + start + n);
		sl.push_back(x_range);
	}
	return sl;
}

ConfigurationTensor<> select_unique(const ConfigurationTensor<>& A, size_t n) {
	ConfigurationTensor<> B;
	map<Configuration<>, size_t> m;
	size_t n_saved = 0;
	for (size_t i = 0; i < A.size(); ++i) {
		const auto& a = A[i];
		// if not in map
		if (m.find(a) == m.end()) {
			m[a] = 1;
			B.push_back(a);
			n_saved++;
		}
		if (n_saved == n) { break; }
	}
	return B;
}

ConfigurationTensor<> sortForEnergy(const ConfigurationTensor<>& A, const vector<double>& E) {

	/// resort best to worst
	ConfigurationTensor<> B = A;
	std::vector<std::size_t> p(A.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](size_t i, size_t j) { return (E[i] < E[j]); });
	std::transform(p.begin(), p.end(), B.begin(),
		[&](std::size_t i) { return A[i]; });
	return B;
}

vector<double> integrateEnergy(const ConfigurationTensor<>& A,
	const vector<double>& E) {

	double beta = 1.;
	map<Configuration<>, double> m;
	for (size_t i = 0; i < A.size(); ++i) {
		if (m.find(A[i]) == m.end()) {
			m[A[i]] = -exp(-beta*E[i]);
		} else {
			m[A[i]] += -exp(-beta*E[i]);
		};
	}
	vector<double> E2;
	for (const auto& a : A) {
		E2.push_back(m[a]);
	}
	return E2;
}

vector<double> evaluateGrid(const ConfigurationTensor<>& A,
	const Configuration<>& idx,
	function<double(const Configuration<>& c)> f,
	pair<Configuration<>, double>& optimal,
	map<Configuration<>, double>& results,
	size_t verboseness) {

	vector<double> E(A.size());
	for (size_t i = 0; i < A.size(); ++i) {
		Configuration<> c = resort(A[i], idx);
//		if (results.find(c) == results.end()) {
			E[i] = f(c);
			if (E[i] < optimal.second) {
				optimal.first = c;
				optimal.second = E[i];
				if (verboseness > 1) { cout << "New record: " << c << " | " << E[i] << " | f_evals: " << results.size() << endl; }
//				cout << "New record: " << E[i] << " | f_evals: " << results.size() << endl;
			}
			results[c] = E[i];
/*		} else {
			E[i] = results[c];
		}*/
//		cout << c << " | " << E[i] << "\n";
	}
//	cout << endl;
	return E;
}

Configuration<> optimize(ConfigurationTree<>& Psi,
	function<double(const Configuration<>& c)> f,
	const Tree& tree, size_t n_sweep, size_t verboseness) {

	map<Configuration<>, double> results;
	pair<Configuration<>, double> optimal;
	optimal.second = 1e99;

	for (size_t sweep = 0; sweep < n_sweep; ++sweep) {
		/// bottom-up
		if (verboseness >= 1) { cout << "Sweep: " << sweep << endl; }
		for (const Node& node: tree) {
			if (node.isToplayer()) { continue; }
			if (node.shape().lastDimension() == node.shape().lastBefore()) { continue; } /// nothing to optimize here
			/// build tensor of all possible configurations
			Configuration<> idx;
			ConfigurationTensor<> A = cartesianProduct(idx, Psi, node);

			/// calculate energies for everything
			auto E = evaluateGrid(A, idx, f, optimal, results, verboseness);

			/// Sort for energy, cut off relevant part, select n best
			size_t start = 0;
			size_t number = Psi.idx_up_[node].size();
			ConfigurationTensor<> B = slice(A, start, number);

			E = integrateEnergy(B, E);

			B = sortForEnergy(B, E);

			/// select the first n unique vectors
			Psi.up_[node] = select_unique(B, node.shape().lastDimension());
		}

		/// top-down
		for (int n = tree.nNodes() - 1; n >= 0; --n) {
			const Node& hole = tree.getNode(n);
			if (hole.isToplayer()) { continue; }
			if (hole.isBottomlayer()) { continue; }
			const Node& node = hole.parent();

			Configuration<> idx;
			ConfigurationTensor<> A = cartesianProduct(idx, Psi, node);

			/// calculate energies for everything
			auto E = evaluateGrid(A, idx, f, optimal, results, verboseness);

			/// Sort for energy, cut off relevant part, select n best
			auto idx_map = findIndices(Psi.idx_down_[hole], idx);
			ConfigurationTensor<> B = sliceDown(A, idx_map);

			E = integrateEnergy(A, E);

			B = sortForEnergy(B, E);

			Psi.down_[hole] = select_unique(B, hole.shape().lastDimension());
		}
	}

	if (verboseness > 0) { cout << "Final result: " << optimal.first << " | " << optimal.second << endl; }
	if (verboseness > 0) { cout << "Number of unique function evaluations: " << results.size() << endl; }
	return optimal.first;
}

size_t to_integer(const Configuration<>& c) {
	size_t r{0};
	size_t factor = 1;
	for (const auto& x: c) {
		r += factor * x;
		factor *= 2;
	}
	return r;
}

void split_integers(vector<size_t>& vec, vector<size_t>& tmp, const Configuration<>& c, size_t n) {
	size_t N = c.size() / n; /// bits / integer
	for (size_t l = 0; l < n; ++l) {
		copy(c.begin() + l * N, c.begin() + (l + 1) * N, tmp.begin());
		vec[l] = to_integer(tmp);
	}
}

vector<size_t> split_integers(const Configuration<>& c, size_t n) {
	/**
	 * n: number of integers in c
	 */
	vector<size_t> vec(n);
	size_t N = c.size() / n; /// bits / integer
	vector<size_t> x(N);
	split_integers(vec, x, c, n);
	return vec;
}

double to_double(size_t i, size_t max_val) {
	return ((double) i / (double) max_val);
}

void split_doubles(vector<double>& xs, vector<size_t>& x_int,
	vector<size_t>& tmp, const Configuration<>& c,
	size_t n) {
	/**
	 * n: number of integers in c
	 */
	size_t N = c.size() / n; /// bits / integer
	split_integers(x_int, tmp, c, n);
	size_t max_val = pow(2, N) - 1;
	for (size_t l = 0; l < xs.size(); ++l) {
		xs[l] = to_double(x_int[l], max_val);
	}
}

vector<double> split_doubles(const Configuration<>& c, size_t n) {
	/**
	 * n: number of integers in c
	 */
	size_t N = c.size() / n; /// bits per integer
	vector<double> xs(n);
	vector<size_t> x_int(n);
	vector<size_t> single_int(N);
	split_doubles(xs, x_int, single_int, c, n);
	return xs;
}
