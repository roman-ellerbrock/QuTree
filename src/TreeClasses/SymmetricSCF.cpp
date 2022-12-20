//
// Created by Roman Ellerbrock on 12/20/22.
//

#include "TreeClasses/SymmetricSCF.h"

template class Configuration<size_t>;

template <typename T>
ostream& operator<<(ostream& os, const Configuration<T>& c) {
	for (size_t i = 0; i < c.size() - 1; ++i) {
		os << c[i] << " ";
	}
	cout << c.back();
	return os;
}

template ostream& operator<<(ostream& os, const Configuration<size_t>&);

template <typename T>
ostream& operator<<(ostream& os, const ConfigurationTensor<T>& c) {
	for (size_t i = 0; i < c.size() - 1; ++i) {
		os << c[i] << "\n";
	}
	cout << c.back();
	return os;
}

template ostream& operator<<(ostream& os, const ConfigurationTensor<size_t>& c);

template <typename T>
ostream& operator<<(ostream& os, ConfigurationTree<T>& a) {
	cout << "up:\n";
	for (const auto x : a.up_) {
		cout << x << endl;
	}
	cout << "down:\n";
	for (const auto x : a.down_) {
		cout << x << endl;
	}
	return os;
}
