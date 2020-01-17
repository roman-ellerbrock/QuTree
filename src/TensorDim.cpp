#include"TensorDim.h"

TensorABC::TensorABC(size_t k, vector<size_t> dim) {
	assert(k < dim.size());
	assert(k >= 0);

	before = 1;
	active = dim[k];
	behind = 1;

	for (size_t i = 0; i < k; i++) {
		before *= dim[i];
	}
	for (size_t i = k + 1; i < dim.size(); i++) {
		behind *= dim[i];
	}

	total = before * active * behind;
}

TensorDim::TensorDim(const vector<size_t>& dim, size_t ntensor_) {
	Initialize(dim, ntensor_);
}

TensorDim::TensorDim(istream& is) {
	ReadDim(is);
}

TensorDim::TensorDim(const string& file) {
	ifstream is(file);
	ReadDim(is);
}

void TensorDim::Initialize(const vector<size_t>& dim, size_t ntensor_) {
	assert(dim.size() > 0);
	assert(ntensor_ > 0);
	f = dim.size();
	abc.clear();
	for (size_t i = 0; i < dim.size(); i++) {
		abc.push_back(TensorABC(i, dim));
	}
	dimpart = abc[0].gettotal();
	ntensor = ntensor_;
	dimtot = dimpart * ntensor_;
}

void TensorDim::WriteBin(ofstream& os) const {
	// Write marker
	os.write("TDIM", 4);

	// Write dof
	int32_t n_write = getntensor();
	os.write((char *) &n_write, sizeof(int32_t));

	// Write dof
	int32_t f_write = f;
	os.write((char *) &f_write, sizeof(int32_t));

	// Write active-dims
	for (size_t k = 0; k < f; k++) {
		int32_t act = Active(k);
		os.write((char *) &act, sizeof(int32_t));
	}
}

void TensorDim::info(ostream& os) const {
	cout << "ntensor = " << getntensor() << "\n";
	cout << "F = " << F() << "\n";
	for (size_t k = 0; k < F(); k++)
		cout << "active = " << Active(k) << "\n";
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
	f = f_read;
	for (size_t k = 0; k < f; k++) {
		is.read((char *) &dim_now, sizeof(int32_t));
		dim_read.push_back(dim_now);
	}

	// Create this Tensor
	(*this) = TensorDim(dim_read, ntens);
}

vector<size_t> TensorDim::getdimlist() const {
	vector<size_t> dimlist;
	for (size_t i = 0; i < f; i++) {
		dimlist.push_back(Active(i));
	}
	return dimlist;
}

const TensorABC& TensorDim::getabc(size_t k) {
	assert(k < f);
	return abc[k];
}

size_t TensorDim::Active(size_t k) const {
	assert(k < f);
	return abc[k].getactive();
}

size_t TensorDim::Behind(size_t k) const {
	assert(k < f);
	return abc[k].getbehind();
}

size_t TensorDim::Before(size_t k) const {
	assert(k < f);
	return abc[k].getbefore();
}

void TensorDim::setntensor(size_t newntensor) {
	dimtot /= ntensor;
	ntensor = newntensor;
	dimtot *= ntensor;
}

void TensorDim::setactive(size_t act, size_t k) {
	assert(k < f);

	vector<size_t> dim;
	for (int k = 0; k < F(); k++) {
		dim.push_back(Active(k));
	}

	dim[k] = act;
	Initialize(dim, ntensor);
}

void TensorDim::print(ostream& os)const {
	os << "d = " << f;
	for (size_t k = 0; k < f; ++k) {
		os << "\t" << k << "\t" << Active(k) << endl;
	}
}
