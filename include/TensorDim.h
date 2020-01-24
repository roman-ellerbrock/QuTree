#pragma once
#include "stdafx.h"

class TensorABC {
// class TensorABC holds contraction of a tensor of n-th order
public:
	TensorABC()
		: before_(0), after_(0), active_(0), total_(0) {}

	TensorABC(size_t k, vector<size_t> dim);

	inline size_t getbefore() const { return before_; }

	inline size_t getactive() const { return active_; }

	inline size_t getafter() const { return after_; }

	inline size_t gettotal() const { return total_; }

private:
	size_t before_;
	size_t after_;
	size_t active_;
	size_t total_;
};

/* *
 * \class TensorDim
 * \ingroup QD-lib
 * \brief This class manages the dimensions of a Tensor.
 * 
 * The class manages dimensions and super-index mappings in a n-th
 * order tensor.
 * \image html mctdh++_TensorDim_1.png
 * \image latex mctdh++_TensorDim_1.eps
 *
 * */

class TensorDim
	// TensorDim (short: TDim) is the dimension-class for a n-th order Tensor
{
public:
	TensorDim()
		: ntensor_(0), dimpart_(0), dimtot_(0), f_(0) {}

	explicit TensorDim(const vector<size_t>& dim, size_t ntensor_);

	explicit TensorDim(istream& is);

	explicit TensorDim(const string& file);

	~TensorDim() = default;

	void Initialize(const vector<size_t>& dim, size_t ntensor);

	void Write(ostream& os) const;

	void Write(const string& filename) const;

	void ReadDim(istream& is);

	inline size_t getf() const { return f_; }

	inline size_t F() const { return f_; }

	inline size_t getdimtot() const { return dimtot_; }

	inline size_t getdimpart() const { return dimpart_; }

	inline size_t getntensor() const { return ntensor_; }

	void setntensor(size_t newntensor);
	void setactive(size_t act, size_t k);

	vector<size_t> getdimlist() const;

	void print(ostream& os = cout) const;

	const TensorABC& getabc(size_t k);
	size_t Active(size_t k) const;
	size_t After(size_t k) const;
	size_t Before(size_t k) const;

protected:
	size_t f_;
	size_t dimtot_;
	size_t dimpart_;
	size_t ntensor_;
	vector<TensorABC> abc_;
};

ostream& operator<<(ostream& os, const TensorDim& tdim);
istream& operator>>(istream& is, TensorDim& tdim);
bool operator==(const TensorDim& tdima, const TensorDim& tdimb);
bool operator!=(const TensorDim& tdima, const TensorDim& tdimb);


