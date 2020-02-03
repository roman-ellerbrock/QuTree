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
		: before_(0), after_(0), active_(0), total_(0) {}

	TensorABC(size_t k, vector<size_t> dim);

	inline size_t GetBefore() const { return before_; }

	inline size_t GetActive() const { return active_; }

	inline size_t GetAfter() const { return after_; }

	inline size_t GetTotal() const { return total_; }

private:
	size_t before_;
	size_t after_;
	size_t active_;
	size_t total_;
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
		: nTensor_(0), dimPart_(0), dimTot_(0), order_(0) {}

	explicit TensorDim(const vector<size_t>& dim, size_t ntensor_);

	explicit TensorDim(istream& is);

	explicit TensorDim(const string& file);

	~TensorDim() = default;

	void Initialize(const vector<size_t>& dim, size_t ntensor);

	void Write(ostream& os) const;

	void Write(const string& filename) const;

	void ReadDim(istream& is);

	inline size_t GetOrder() const { return order_; }

	inline size_t GetDimTot() const { return dimTot_; }

	inline size_t GetDimPart() const { return dimPart_; }

	inline size_t GetNumTensor() const { return nTensor_; }

	void SetNumTensor(size_t newntensor);
	void SetActive(size_t act, size_t k);

	vector<size_t> GetDimList() const;

	void print(ostream& os = cout) const;

	const TensorABC& getabc(size_t k);
	size_t Active(size_t k) const;
	size_t After(size_t k) const;
	size_t Before(size_t k) const;

protected:
	size_t order_;
	size_t dimTot_;
	size_t dimPart_;
	size_t nTensor_;
	vector<TensorABC> abc_;
};

ostream& operator<<(ostream& os, const TensorDim& tdim);
istream& operator>>(istream& is, TensorDim& tdim);
bool operator==(const TensorDim& tdima, const TensorDim& tdimb);
bool operator!=(const TensorDim& tdima, const TensorDim& tdimb);


