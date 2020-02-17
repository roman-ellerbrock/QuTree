#pragma once
#include "stdafx.h"

class TensorABC {
/**
 * \class TensorABC
 * \ingroup Core
 * \brief This class holds contracted (multi-shape) dimensions of a tensor.
 */
public:
	TensorABC()
		: before_(0), after_(0), active_(0) {}

	TensorABC(size_t k, vector<size_t> dim);

	inline size_t GetBefore() const { return before_; }

	inline size_t GetActive() const { return active_; }

	inline size_t GetAfter() const { return after_; }

private:
	size_t before_;
	size_t after_;
	size_t active_;
};

class TensorDim
/**
 * \class TensorDim
 * \ingroup Core
 * \brief This class manages the dimensions of a Tensor.
 *
 * The class manages dimensions and super-index mappings in a n-th
 * order tensor.
 * \image html TensorDim_1.png
 * \image latex TensorDim_1.eps_
 *
 * Usage:
 * TensorDim tdim({2, 3, 4}, 5);
 * */
{
public:
	TensorDim()
		: dimTot_(0) {}

	explicit TensorDim(const initializer_list<size_t>& dims);

	explicit TensorDim(const vector<size_t>& dim);

	explicit TensorDim(istream& is);

	explicit TensorDim(const string& file);

	~TensorDim() = default;

	void Initialize(const vector<size_t>& dim);

	void Write(ostream& os) const;

	void Write(const string& filename) const;

	void ReadDim(istream& is);

	inline size_t GetOrder() const { return abc_.size(); }

	inline size_t GetLastIdx() const { return abc_.size() - 1; }

	inline size_t GetDimTot() const { return dimTot_; }

	inline size_t LastBefore() const { return abc_.back().GetBefore(); }

	inline size_t LastActive() const { return abc_.back().GetActive(); }

	void SetActive(size_t act, size_t k);

	vector<size_t> GetDimList() const;

	void print(ostream& os = cout) const;

	const TensorABC& getabc(size_t k);
	size_t Before(size_t k) const;
	size_t Active(size_t k) const;
	size_t After(size_t k) const;

protected:
	size_t dimTot_;
	vector<TensorABC> abc_;
};

ostream& operator<<(ostream& os, const TensorDim& tdim);
istream& operator>>(istream& is, TensorDim& tdim);
bool operator==(const TensorDim& tdima, const TensorDim& tdimb);
bool operator!=(const TensorDim& tdima, const TensorDim& tdimb);

