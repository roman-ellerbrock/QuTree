//
// Created by Roman Ellerbrock on 12/20/22.
//

#ifndef SYMMETRICSCF_H
#define SYMMETRICSCF_H
#include "NodeAttribute.h"
#include <iostream>
#include <random>

template <typename T = size_t>
class Configuration : public vector<T> {
public:
	using Base = vector<T>;

	Configuration() = default;
	~Configuration() = default;

	Configuration(const vector<T>& a) : Base(a) {}

	Configuration(const initializer_list<T>& list)
	: Base(list) {}

	Configuration operator*(const Configuration& r) const {
		Configuration<T> a(*this);
		a.insert(a.end(), r.begin(), r.end());
		return a;
	}

	bool operator==(const Configuration& r) const {
		if (Base::size() != r.size()) { return false; }
		for (size_t i = 0; i < r.size(); ++i) {
			if (this->operator[](i) != r[i]) { return false; }
		}
		return true;
	}

};

template <typename T>
ostream& operator<<(ostream& os, const Configuration<T>& c);

template <typename T = size_t>
class ConfigurationTensor : public vector<Configuration<T>> {
public:
	using Base = vector<Configuration<T>>;
	ConfigurationTensor() = default;
	~ConfigurationTensor() = default;

	explicit ConfigurationTensor(size_t len)
		: Base(len, Configuration<T>()) {
	}

	ConfigurationTensor(const initializer_list<Configuration<T>>& list)
		: Base(list) {
	}

	ConfigurationTensor(const initializer_list<initializer_list<T>>& list) {
		for (const auto& x: list) {
			Base::emplace_back(x);
		}
	}

	bool operator==(const ConfigurationTensor& r) const {
		if (Base::size() != r.size()) { return false; }
		for (size_t i = 0; i < r.size(); ++i) {
			if (this->operator[](i) != r[i]) { return false; }
		}
		return true;
	}

	ConfigurationTensor operator+(const ConfigurationTensor& r) const {
		ConfigurationTensor a(*this);
		a.insert(a.end(), r.begin(), r.end());
		return a;
	}

	ConfigurationTensor operator*(const ConfigurationTensor& r) const {
		ConfigurationTensor a;
		a.resize(this->size() * r.size());
		size_t idx = 0;
		for (const Configuration<>& x : (*this)) {
			for (size_t i = 0; i < r.size(); ++i) {
				a[idx++] = x * r[i];
			}
		}
		if (this->empty()) {
			a = r;
		}
		if (r.empty()) {
			a = *this;
		}
		return a;
	}

};

template <typename T>
ostream& operator<<(ostream& os, const ConfigurationTensor<T>& c);

template <typename T = size_t>
class ConfigurationTree {
public:
	ConfigurationTree() = default;
	~ConfigurationTree() = default;

	explicit ConfigurationTree(const Tree& tree) {
		initialize(tree);
	}

	void initialize(const Tree& tree) {
		up_.clear();
		down_.clear();
		for (const Node& node : tree) {
			up_.emplace_back(ConfigurationTensor<T>());
			down_.emplace_back(ConfigurationTensor<T>());
			idx_up_.emplace_back(Configuration<T>());
			idx_down_.emplace_back(Configuration<T>());
		}

		/// build index up
		for (const Node& node : tree) {
			if (node.isToplayer()) { continue; }
			Configuration<>& A = idx_up_[node];
			if (node.isBottomlayer()) {
				size_t mode = node.getLeaf().mode();
				A = Configuration<>({mode});
			} else {
				for (size_t k = 0; k < node.nChildren(); ++k) {
					const Node& child = node.child(k);
					A = A * idx_up_[child];
				}
			}
		}

		/// build index down
		for (int n = tree.nNodes() - 1; n >= 0; --n) {
			const Node& hole = tree.getNode(n);
			if (hole.isToplayer()) { continue; }
			Configuration<>& A = idx_down_[hole];
			const Node& node = hole.parent();

			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child = node.child(k);
				if (child.address() == hole.address()) { continue; }
				A = A * idx_up_[child];
			}

			A = A * idx_down_[node];
		}
	}

	void print(const Tree& tree) const {
		cout << "up:\n";
		for (const Node& node : tree) {
			node.info();
			cout << "idx = " << idx_up_[node] << endl;
			cout << up_[node] << endl;
		}
		cout << "down:\n";
		for (const Node& node : tree) {
			node.info();
			cout << "idx = " << idx_down_[node] << endl;
			cout << down_[node] << endl;
		}
	}

	NodeAttribute<ConfigurationTensor<T>> up_;
	NodeAttribute<ConfigurationTensor<T>> down_;
	NodeAttribute<Configuration<T>> idx_up_;
	NodeAttribute<Configuration<T>> idx_down_;
};

template <typename T = size_t>
ostream& operator<<(ostream& os, const ConfigurationTree<T>& a);

ConfigurationTensor<> randomConfigurationTensor(const TensorShape& shape, mt19937& gen);
ConfigurationTree<> randomConfigurationTree(const Tree& tree, mt19937& gen);

ConfigurationTensor<> bottomTensor(const TensorShape& shape);

Configuration<> optimize(ConfigurationTree<>& Psi,
	function<double(const Configuration<>& c)> f,
	const Tree& tree, size_t n_sweep = 3, size_t verboseness = 1);

size_t to_integer(const Configuration<>& c);
vector<size_t> split_integers(const Configuration<>& c, size_t n);
double to_double(size_t i, size_t max_val);
void split_doubles(vector<double>& xs, vector<size_t>& x_int,
	vector<size_t>& tmp, const Configuration<>& c,
	size_t n);
vector<double> split_doubles(const Configuration<>& c, size_t n);

#endif //SYMMETRICSCF_H
