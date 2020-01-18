#pragma once
#include "stdafx.h"

class TensorABC {
// class TensorABC holds contraction of a tensor of n-th order
public:
	TensorABC()
		: before(0), after(0), active(0), total(0) {}

	TensorABC(size_t k, vector<size_t> dim);

	inline size_t getbefore() const { return before; }

	inline size_t getactive() const { return active; }

	inline size_t getafter() const { return after; }

	inline size_t gettotal() const { return total; }

private:
	size_t before;
	size_t after;
	size_t active;
	size_t total;
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
		: ntensor(0), dimpart(0), dimtot(0), f(0) {}

	explicit TensorDim(const vector<size_t>& dim, size_t ntensor_);

	explicit TensorDim(istream& is);

	explicit TensorDim(const string& file);

	~TensorDim() = default;

	void Initialize(const vector<size_t>& dim, size_t ntensor);

	void Write(ofstream& os) const;

	void Write(const string& filename) const;

	void ReadDim(istream& is);

	inline size_t getf() const { return f; }

	inline size_t F() const { return f; }

	inline size_t getdimtot() const { return dimtot; }

	inline size_t getdimpart() const { return dimpart; }

	inline size_t getntensor() const { return ntensor; }

	void setntensor(size_t newntensor);
	void setactive(size_t act, size_t k);

	vector<size_t> getdimlist() const;

	friend bool operator==(const TensorDim& tdima, const TensorDim& tdimb) {
		if (tdima.F() != tdimb.F()) { return false; }
		for (size_t k = 0; k < tdima.F(); k++) {
			if (tdima.Active(k) != tdimb.Active(k)) { return false; }
		}
		return (tdima.getntensor() == tdimb.getntensor());
	}

	friend bool operator!=(TensorDim& tdima, TensorDim& tdimb) {
		return !(tdima == tdimb);
	}

	void print(ostream& os = cout) const;

	const TensorABC& getabc(size_t k);
	size_t Active(size_t k) const;
	size_t After(size_t k) const;
	size_t Before(size_t k) const;

protected:
	size_t f;
	size_t dimtot;
	size_t dimpart;
	size_t ntensor;
	vector<TensorABC> abc;
};

