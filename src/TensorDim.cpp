#include"TensorDim.h"

TensorABC::TensorABC(size_t k, vector<size_t> dim)
	: before_(1), active_(1), after_(1), total_(1) {
	assert(k < dim.size());
	assert(k >= 0);

	before_ = 1;
	active_ = dim[k];
	after_ = 1;

	for (size_t i = 0; i < k; i++) {
		before_ *= dim[i];
	}
	for (size_t i = k + 1; i < dim.size(); i++) {
		after_ *= dim[i];
	}

	total_ = before_ * active_ * after_;
}

TensorDim::TensorDim(const vector<size_t>& dim, size_t ntensor)
	: TensorDim() {
	Initialize(dim, ntensor);
}

TensorDim::TensorDim(istream& is)
	: TensorDim() {
	ReadDim(is);
}

TensorDim::TensorDim(const string& file)
	: TensorDim() {
	ifstream is(file);
	ReadDim(is);
}

void TensorDim::Initialize(const vector<size_t>& dim, size_t ntensor) {
	assert(!dim.empty());
	assert(ntensor > 0);
	f_ = dim.size();
	abc_.clear();
	for (size_t i = 0; i < dim.size(); i++) {
		abc_.emplace_back(TensorABC(i, dim));
	}
	dimpart_ = abc_[0].gettotal();
	ntensor_ = ntensor;
	dimtot_ = dimpart_ * ntensor_;
}

void TensorDim::Write(ostream& os) const {
	// Write marker
	os.write("TDIM", 4);

	// Write dof
	int32_t n_write = getntensor();
	os.write((char *) &n_write, sizeof(int32_t));

	// Write dof
	int32_t f_write = f_;
	os.write((char *) &f_write, sizeof(int32_t));

	// Write active-dims
	for (size_t k = 0; k < f_; k++) {
		int32_t act = Active(k);
		os.write((char *) &act, sizeof(int32_t));
	}
}

void TensorDim::Write(const string& filename) const {
	ofstream os(filename);
	Write(os);
}

void TensorDim::ReadDim(istream& is) {
	// Check if binary string contains a TDim
	char check[5];
	is.read(check, 4);
	string s_check(check, 4);
	string s_key("TDIM");
	assert(s_check == s_key);

	// Read n-tensor
	int32_t n_read;
	is.read((char *) &n_read, sizeof(n_read));
	size_t ntens = n_read;

	// Read dof
	int32_t f_read;
	is.read((char *) &f_read, sizeof(f_read));

	// Read vector with dimensions of Tensor
	vector<size_t> dim_read;
	int32_t dim_now;
	f_ = f_read;
	for (size_t k = 0; k < f_; k++) {
		is.read((char *) &dim_now, sizeof(int32_t));
		dim_read.push_back(dim_now);
	}

	// Create this Tensor
	(*this) = TensorDim(dim_read, ntens);
}

vector<size_t> TensorDim::getdimlist() const {
	vector<size_t> dimlist;
	for (size_t i = 0; i < f_; i++) {
		dimlist.push_back(Active(i));
	}
	return dimlist;
}

const TensorABC& TensorDim::getabc(size_t k) {
	assert(k < f_);
	return abc_[k];
}

size_t TensorDim::Active(size_t k) const {
	assert(k < f_);
	return abc_[k].getactive();
}

size_t TensorDim::After(size_t k) const {
	assert(k < f_);
	return abc_[k].getafter();
}

size_t TensorDim::Before(size_t k) const {
	assert(k < f_);
	return abc_[k].getbefore();
}

void TensorDim::setntensor(size_t newntensor) {
	dimtot_ /= ntensor_;
	ntensor_ = newntensor;
	dimtot_ *= ntensor_;
}

void TensorDim::setactive(size_t act, size_t k) {
	assert(k < f_);

	vector<size_t> dim(F());
	for (int l = 0; l < F(); l++) {
		dim.emplace_back(Active(l));
	}

	dim[k] = act;
	Initialize(dim, ntensor_);
}

void TensorDim::print(ostream& os) const {
	if (f_ > 0) {
		os << "(";
		for (size_t k = 0; k < f_ - 1; ++k) {
			os << Active(k) << ", ";
		}
		os << Active(f_ - 1) << "); ";
		os << ntensor_ << endl;
	} else {
		os << "( ); " << ntensor_ << endl;
	}
}

ostream& operator<<(ostream& os, const TensorDim& tdim) {
	tdim.print(os);
	return os;
}

istream& operator>>(istream& is, TensorDim& tdim) {
	tdim.ReadDim(is);
	return is;
}

bool operator==(const TensorDim& tdima, const TensorDim& tdimb) {
	if (tdima.F() != tdimb.F()) { return false; }
	for (size_t k = 0; k < tdima.F(); k++) {
		if (tdima.Active(k) != tdimb.Active(k)) { return false; }
	}
	return (tdima.getntensor() == tdimb.getntensor());
}

bool operator!=(const TensorDim& tdima, const TensorDim& tdimb) {
	return !(tdima == tdimb);
}

