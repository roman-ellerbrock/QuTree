#pragma once
#include "stdafx.h"


class TensorShape : public vector<size_t>
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
	TensorShape()
		: totalDimension_(0) {}

	TensorShape(const initializer_list<size_t>& dims);

	explicit TensorShape(const vector<size_t>& dim);

	explicit TensorShape(istream& is);

	explicit TensorShape(const string& file);

	~TensorShape() = default;

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

TensorShape replaceDimension(TensorShape shape, size_t target, size_t new_dimension);

ostream& operator<<(ostream& os, const TensorShape& tdim);
istream& operator>>(istream& is, TensorShape& tdim);
bool operator==(const TensorShape& tdima, const TensorShape& tdimb);
bool operator!=(const TensorShape& tdima, const TensorShape& tdimb);
