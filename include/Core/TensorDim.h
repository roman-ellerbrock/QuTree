#pragma once
#include "stdafx.h"


class TensorDim : public vector<size_t>
/**
 * \class TensorDim
 * \ingroup Core
 * \brief This class manages the dimensions of a Tensor.
 *
 * This class holds the tensor dimensions n=(n_1,...,n_order).
 * It is a decorated vector of (unsigned) integers (the dimensions).
 * The class, furthermore, manages contracted dimensions for quick
 * reshapings.
 * \image html TensorDim_1.png
 * \image latex TensorDim_1.eps_
 *
 * Usage:
 * TensorDim tdim({2, 3, 4});
 * */
{
public:
	TensorDim()
		: totalDimension_(0) {}

	TensorDim(const initializer_list<size_t>& dims);

	explicit TensorDim(const vector<size_t>& dim);

	explicit TensorDim(istream& is);

	explicit TensorDim(const string& file);

	~TensorDim() = default;

	void Initialize(const vector<size_t>& dim);

	void Write(ostream& os) const;

	void Write(const string& filename) const;

	void ReadDim(istream& is);

	inline size_t order() const { return size(); }

	inline size_t lastIdx() const { return size() - 1; }

	inline size_t totalDimension() const { return totalDimension_; }

	inline size_t lastBefore() const { return before_.back(); }

	inline size_t lastDimension() const { return back(); }

	void setDimension(size_t act, size_t k);

	vector<size_t> dimensions() const;

	void print(ostream& os = cout) const;

	size_t before(size_t k) const;
	size_t after(size_t k) const;

protected:
	size_t totalDimension_;
	vector<size_t> before_;
	vector<size_t> after_;
};

ostream& operator<<(ostream& os, const TensorDim& tdim);
istream& operator>>(istream& is, TensorDim& tdim);
bool operator==(const TensorDim& tdima, const TensorDim& tdimb);
bool operator!=(const TensorDim& tdima, const TensorDim& tdimb);

