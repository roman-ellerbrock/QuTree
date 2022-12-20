//
// Created by Roman Ellerbrock on 12/20/22.
//

#ifndef SYMMETRICSCF_H
#define SYMMETRICSCF_H
#include "NodeAttribute.h"
#include <iostream>

template <typename T = size_t>
class Configuration : public vector<T> {
public:
	using Base = vector<T>;

	Configuration() = default;
	~Configuration() = default;

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
		}
	}

	NodeAttribute<ConfigurationTensor<T>> up_;
	NodeAttribute<ConfigurationTensor<T>> down_;
};

template <typename T = size_t>
ostream& operator<<(ostream& os, ConfigurationTree<T>& a);

class SymmetricSCF {
public:
};


#endif //SYMMETRICSCF_H
