//
// Created by Roman Ellerbrock on 2/25/22.
//

#ifndef SPARSE_VECTOR_H
#define SPARSE_VECTOR_H
#include <vector>
#include <map>
#include <iostream>

/**
 * \class sparse_vector
 *
 * sparse_vector is an ordered map that provides sequential access.
 */

template<class A>
class sparse_vector {
public:
	sparse_vector() = default;
	~sparse_vector() = default;

	void push_back(size_t dense_address, const A a) {
		objects_.push_back(a);
		map_[dense_address] = objects_.size() - 1;
	}

	bool contains(size_t dense_address) {
		if (map_.find(dense_address) == map_.end()) return false;
		return true;
	}

	const A& operator[](size_t dense_address) const {
		size_t sparse_idx = map_.at(dense_address);
		return objects_[sparse_idx];
	}

	auto begin() { return objects_.begin(); }
	auto end() { return objects_.end(); }

	auto begin() const { return objects_.begin(); }
	auto end() const { return objects_.end(); }

	void clear() {
		map_.clear();
		objects_.clear();
	}

	/// pointers to objects
	std::vector<A> objects_;
	/// dense address -> sparse address
	std::map<size_t, size_t> map_;
};

template <class T>
std::ostream& operator<<(std::ostream& os, const sparse_vector<T>& vec) {
	for (auto p : vec.map_) {
		os << p.first << " - " << vec.objects_.at(p.second) << std::endl;
	}
	return os;
}


#endif //SPARSE_VECTOR_H
