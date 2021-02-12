#include "TreeShape/Leaf.h"
#include "TreeShape/LeafTypes/DVRBasis.h"
#include "TreeShape/LeafTypes/HO_Basis.h"
#include "TreeShape/LeafTypes/FFTGrid.h"
#include "TreeShape/LeafTypes/LegendrePolynomials.h"
#include "TreeShape/LeafTypes/SpinGroup.h"

Leaf::Leaf()
	: dim_(-1), type_(0), mode_(-1), subType_(0), parent_(nullptr), nodeType_(0) {}

Leaf::Leaf(istream& file, AbstractNode *up_, NodePosition position_)
	: parent_(up_), position_(move(position_)), type_(0), subType_(0), nodeType_(0), dim_(0), mode_(0) {
	// read basis size information
	file >> dim_;
	file >> type_;
	file >> mode_;
	assert(dim_ > 0);
	assert(type_ >= 0);
//	cout << "Leaf: " << dim_ << " " << type << " " << mode << endl;
	CreatePrimitiveBasis(type_, subType_, dim_);
}

Leaf::Leaf(const Leaf& old)
	: dim_(old.dim_), type_(old.type_), mode_(old.mode_), subType_(old.subType_),
	  nodeType_(old.nodeType_), parent_(old.parent_), position_(old.position_) {
	CreatePrimitiveBasis(type_, subType_, dim_);
	par_ = old.par_;
	interface_->initialize(par_.omega(), par_.r0(), par_.wfr0(), par_.wfOmega());
}

void Leaf::CreatePrimitiveBasis(size_t type, size_t subtype, size_t dim) {
	// Construct Fundamental Operator class
	if (type == 0) {
		interface_ = make_unique<HO_Basis>(dim);
	} else if (type == 1) {
		interface_ = make_unique<FFTGrid>(dim);
	} else if (type == 2) {
		interface_ = make_unique<LegendrePolynomials>(dim);
	} else if (type == 6) {
		interface_ = make_unique<SpinGroup>(dim);
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
	  mode_(mode), nodeType_(0), parent_(nullptr) {
	CreatePrimitiveBasis(type, subtype, dim);
}

void Leaf::info(ostream& os) const {
	os << "Leaf" << endl;
	position_.info(os);
	os << "mode=" << mode() << endl;
	par_.info(os);
}

void Leaf::write(ostream& os) const {
	for (size_t l = 0; l < position_.layer(); l++) { os << "\t"; }
	os << dim_ << "\t" << type_ << "\t" << mode_ << "\n";
}

void Leaf::update(const NodePosition& p) {
	updatePosition(p);
}

void Leaf::updatePosition(const NodePosition& p) {
	position_ = p;
}

