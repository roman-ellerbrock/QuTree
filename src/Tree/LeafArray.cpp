#include "Tree/LeafArray.h"

LeafArray::LeafArray(Node& root) {
	address_.clear();
	leaves_.clear();
	address_.resize(root.nLeaves());

	Node* node = &root;
	while(node) {
		for (Leaf& leaf : node->leaves_) {
			leaves_.emplace_back(reference_wrapper<Leaf>(leaf));
			size_t mode = leaf.api_.basis()->par_.mode_;
			address_[mode] = leaves_.size() - 1;
		}
		node = sweep(node);
	}
}

void LeafArray::push_back(Leaf& leaf) {
	leaves_.emplace_back(reference_wrapper<Leaf>(leaf));

	// Save address_ of this leafical mode
	size_t mode = leaf.api_.basis()->par_.mode_;
	if (mode >= address_.size()) {
		address_.resize(mode + 1);
	}
	int addr = leaves_.size() - 1;
	assert(addr < address_.size());
	address_[mode] = addr;
}


const Leaf& LeafArray::operator[](size_t mode) const {
	// The address_ of mode is saved in the vector "address_"
	assert(mode < address_.size());
	int addr = address_[mode];
	return leaves_[addr].get();
}

Leaf& LeafArray::operator[](size_t mode) {
	// The address_ of mode is saved in the vector "address_"
	assert(mode < address_.size());
	int addr = address_[mode];
	return leaves_[addr].get();
}

void LeafArray::readPars(istream& is) {
	for (size_t i = 0; i < size(); ++i) {
		double par0, par1, par2, par3;
		is >> par0 >> par1 >> par2 >> par3;
		BasisParameters& def = operator[](i).par();
		def.par0_ = par0;
		def.par1_ = par1;
		def.par2_ = par2;
		def.par3_ = par3;
	}
}
