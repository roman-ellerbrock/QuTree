#include "Leaf.h"

Leaf::Leaf()
	: dim_(-1), type_(0), mode_(-1), subType_(0), up_(nullptr), nodeType_(0) {}

Leaf::Leaf(istream& file, AbstractNode *up_, NodePosition position_)
	: up_(up_), position_(move(position_)), type_(0), subType_(0), nodeType_(0), dim_(0), mode_(0) {
	// read basis size information
	file >> dim_;
	file >> type_;
	file >> mode_;
	assert(dim_ > 0);
	assert(type_ >= 0);
//	cout << "Leaf: " << dim << " " << type << " " << mode << endl;
	CreatePrimitiveBasis(type_, subType_, dim_);
}

Leaf::Leaf(const Leaf& old)
	: dim_(old.dim_), type_(old.type_), mode_(old.mode_), subType_(old.subType_),
	  nodeType_(old.nodeType_), up_(old.up_), position_(old.position_) {
	CreatePrimitiveBasis(type_, subType_, dim_);
	par_ = old.par_;
	primitiveBasis_->Initialize(par_.Omega(), par_.R0(), par_.WFR0(), par_.WFOmega());
}

void Leaf::CreatePrimitiveBasis(size_t type, size_t subtype, size_t dim) {
	// Construct Fundamental Operator class
	if (type == 0) {
		primitiveBasis_ = make_unique<HO_Basis>(dim);
	} else if (type == 1) {
		primitiveBasis_ = make_unique<FFTGrid>(dim);
	} else if (type == 2) {
		primitiveBasis_ = make_unique<LegendrePolynomials>(dim);
	} else if (type == 3) {
		primitiveBasis_ = make_unique<BosonNumberBasis>(dim);
	} else if (type == 4) {
		primitiveBasis_ = make_unique<FermionNumberBasis>(dim);
	} else if (type == 5) {
		primitiveBasis_ = make_unique<LogicalBasis>();
	} else if (type == 6) {
		primitiveBasis_ = make_unique<SpinGroup>(dim);
	} else {
		cout << "Error: This Basis Type is not in the known list of "
			 << "Typs. The iplemented ones are: \n"
			 << "0 = Hermite-DVR\n"
			 << "1 = FFT-Grid\n"
			 << "2 = Legendre-DVR\n"
			 << "3 = bosonic occupation numbers\n"
			 << "4 = fermionic occupation numbers\n"
			 << "5 = Logical Basis\n"
			 << "6 = Spin Group\n";
		assert(false);
	}
}

Leaf& Leaf::operator=(const Leaf& old) {
	Leaf phys(old);
	*this = move(phys);
	return *this;
}

Leaf::Leaf(size_t dim, size_t mode, size_t type, size_t subtype,
	PhysPar par)
	: type_(type), subType_(subtype), dim_(dim), par_(par),
	  mode_(mode), nodeType_(0), up_(nullptr) {
	CreatePrimitiveBasis(type, subtype, dim);
}

void Leaf::info(ostream& os) const {
	os << "Leaf" << endl;
	position_.info(os);
	os << "mode=" << Mode() << endl;
	par_.info(os);
}

void Leaf::Write(ostream& os) const {
	for (size_t l = 0; l < position_.Layer(); l++) { os << "\t"; }
	os << dim_ << "\t" << type_ << "\t" << mode_ << "\n";
}

void Leaf::Update(const NodePosition& p) {
	UpdatePosition(p);
}

void Leaf::UpdatePosition(const NodePosition& p) {
	position_ = p;
}

