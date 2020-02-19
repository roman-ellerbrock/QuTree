#include"Core/TensorDim.h"

size_t ContractDimensionsBefore(const vector<size_t>& dim, size_t k) {
	size_t bef = 1;
	for (size_t i = 0; i < k; i++) {
		bef *= dim[i];
	}
	return bef;
}

vector<size_t> ContractDimensionsBefore(const vector<size_t>& dim) {
	vector<size_t> befores(dim.size());
	for (size_t k = 0; k < dim.size(); ++k) {
		befores[k] = ContractDimensionsBefore(dim, k);
	}
	return befores;
}

size_t ContractDimensionsAfter(const vector<size_t>& dim, size_t k) {
	size_t aft = 1;
	for (size_t i = k + 1; i < dim.size(); i++) {
		 aft *= dim[i];
	}
	return aft;
}

vector<size_t> ContractDimensionsAfter(const vector<size_t>& dim) {
	vector<size_t> afters(dim.size());
	for (size_t k = 0; k < dim.size(); ++k) {
		afters[k] = ContractDimensionsAfter(dim, k);
	}
	return afters;
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
	resize(dim.size());
	for (size_t k = 0; k < dim.size(); ++k) {
		this->operator[](k) = dim[k];
	}
	after_ = ContractDimensionsAfter(dim);
	before_ = ContractDimensionsBefore(dim);
	totalDimension_ = after_.front() * front();
}

void TensorDim::Write(ostream& os) const {
	// Write marker
	os.write("TDIM", 4);

	// Write dof
	int32_t f_write = order();
	os.write((char *) &f_write, sizeof(int32_t));

	// Write active_-dims
	for (size_t k = 0; k < order(); k++) {
		int32_t act = dimension(k);
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

vector<size_t> TensorDim::dimensions() const {
	vector<size_t> dimlist;
	for (size_t i = 0; i < order(); i++) {
		dimlist.push_back(dimension(i));
	}
	return dimlist;
}

size_t TensorDim::dimension(size_t k) const {
	return operator[](k);
}

size_t TensorDim::before(size_t k) const {
	assert(k < order());
	return before_[k];
}

size_t TensorDim::after(size_t k) const {
	assert(k < order());
	return after_[k];
}

void TensorDim::setDimension(size_t act, size_t k) {
	vector<size_t> dim = dimensions();
	assert(k < dim.size());
	dim[k] = act;
	Initialize(dim);
}

void TensorDim::print(ostream& os) const {
	if (order() > 0) {
		os << "(";
		for (size_t k = 0; k < order() - 1; ++k) {
			os << dimension(k) << ", ";
		}
		os << dimension(order() - 1);
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
	if (tdima.order() != tdimb.order()) { return false; }
	for (size_t k = 0; k < tdima.order(); k++) {
		if (tdima.dimension(k) != tdimb.dimension(k)) { return false; }
	}
	return true;
}

bool operator!=(const TensorDim& tdima, const TensorDim& tdimb) {
	return !(tdima == tdimb);
}

