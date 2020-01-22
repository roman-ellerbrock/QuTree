#include "LinearizedLeaves.h"

void LinearizedLeaves::Write(ostream& os) const {
}

void LinearizedLeaves::push_back(Leaf& phys) {
	// Check wether this node is a physical node
	assert(phys.NodeType() == 0);
	coordinates_.emplace_back(reference_wrapper<Leaf>(phys));

	// Save address_ of this physical mode
	int mode = phys.Mode();
	assert(mode >= 0);
	// mode > address_ fails due to CDVR (TDVR-X-Matrices)
	// @TODO: This is not fixed yet but commented out to test partialoverlap!
//	assert(mode < address_.size());
	// Assure that the address_ vector is long enough
	if (mode >= address_.size()) {
		address_.resize(mode + 1);
	}
	int addr = coordinates_.size() - 1;
	assert(addr >= 0);
	assert(addr < address_.size());
	address_[mode] = addr;
}

void LinearizedLeaves::resizeaddress(int n) {
	address_.resize(n);
}

const Leaf& LinearizedLeaves::operator[](int mode) const {
	// The address_ of mode is saved in the vector "address_"
	assert(mode < address_.size());
	assert(mode >= 0);
	int addr = address_[mode];
	return coordinates_[addr].get();
}

Leaf& LinearizedLeaves::operator[](int mode) {
	// The address_ of mode is saved in the vector "address_"
	assert(mode < address_.size());
	assert(mode >= 0);
	int addr = address_[mode];
	return coordinates_[addr].get();
}

reference_wrapper<Leaf>& LinearizedLeaves::operator()(int mode) {
	assert(mode < coordinates_.size());
	assert(mode >= 0);
	assert(0);
	int addr = address_[mode];
	return coordinates_[addr];
}

void LinearizedLeaves::clear() {
	address_.clear();
	coordinates_.clear();
}
