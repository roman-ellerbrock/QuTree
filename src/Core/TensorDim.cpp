#include"Core/TensorDim.h"

TensorABC::TensorABC(size_t k, vector<size_t> dim)
	: before_(1), active_(1), after_(1) {
	assert(k < dim.size());

	before_ = 1;
	active_ = dim[k];
	after_ = 1;

	for (size_t i = 0; i < k; i++) {
		before_ *= dim[i];
	}
	for (size_t i = k + 1; i < dim.size(); i++) {
		after_ *= dim[i];
	}
}

TensorDim::TensorDim(const vector<size_t>& dim)
	: TensorDim() {
	Initialize(dim);
}

TensorDim::TensorDim(const initializer_list<size_t>& dims)
	: TensorDim(vector<size_t>(dims)){ }

TensorDim::TensorDim(istream& is)
	: TensorDim() {
	ReadDim(is);
}

TensorDim::TensorDim(const string& file)
	: TensorDim() {
	ifstream is(file);
	ReadDim(is);
}

void TensorDim::Initialize(const vector<size_t>& dim) {
	assert(!dim.empty());
	abc_.clear();
	for (size_t i = 0; i < dim.size(); i++) {
		abc_.emplace_back(TensorABC(i, dim));
	}
	const TensorABC& cdim = abc_.back();
	dimTot_ = cdim.GetActive() * cdim.GetBefore();
}

void TensorDim::Write(ostream& os) const {
	// Write marker
	os.write("TDIM", 4);

	// Write dof
	int32_t f_write = GetOrder();
	os.write((char *) &f_write, sizeof(int32_t));

	// Write active_-dims
	for (size_t k = 0; k < GetOrder(); k++) {
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

	// Read dof
	int32_t f_read;
	is.read((char *) &f_read, sizeof(f_read));

	// Read vector with dimensions of Tensor
	size_t order = f_read;
	vector<size_t> dim_read(order);
	int32_t dim_now;
	for (size_t k = 0; k < order; k++) {
		is.read((char *) &dim_now, sizeof(int32_t));
		dim_read[k]=dim_now;
	}

	// Create this Tensor
	(*this) = TensorDim(dim_read);
}

vector<size_t> TensorDim::GetDimList() const {
	vector<size_t> dimlist;
	for (size_t i = 0; i < GetOrder(); i++) {
		dimlist.push_back(Active(i));
	}
	return dimlist;
}

const TensorABC& TensorDim::getabc(size_t k) {
	assert(k < GetOrder());
	return abc_[k];
}

size_t TensorDim::Active(size_t k) const {
	assert(k < GetOrder());
	return abc_[k].GetActive();
}

size_t TensorDim::Before(size_t k) const {
	assert(k < GetOrder());
	return abc_[k].GetBefore();
}

size_t TensorDim::After(size_t k) const {
	assert(k < GetOrder());
	return abc_[k].GetAfter();
}

void TensorDim::SetActive(size_t act, size_t k) {
	vector<size_t> dim = GetDimList();
	assert(k < dim.size());
	dim[k] = act;
	Initialize(dim);
}

void TensorDim::print(ostream& os) const {
	if (GetOrder() > 0) {
		os << "(";
		for (size_t k = 0; k < GetOrder(); ++k) {
			os << Active(k) << ", ";
		}
		os <<  ")" << endl;
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
	if (tdima.GetOrder() != tdimb.GetOrder()) { return false; }
	for (size_t k = 0; k < tdima.GetOrder(); k++) {
		if (tdima.Active(k) != tdimb.Active(k)) { return false; }
	}
	return true;
}

bool operator!=(const TensorDim& tdima, const TensorDim& tdimb) {
	return !(tdima == tdimb);
}

