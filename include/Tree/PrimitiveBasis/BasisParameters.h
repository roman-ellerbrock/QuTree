//
// Created by Roman Ellerbrock on 11/27/21.
//
#ifndef BASISPARAMETERS_H
#define BASISPARAMETERS_H
#include <iostream>
#include "stdafx.h"

struct BasisParameters {
	/**
	 * \class Parameter container of a leaf.
	 */
	~BasisParameters() = default;

	void readTail(istream& file) {
		file >> par0_;
		file >> par1_;
		file >> par2_;
		file >> par3_;
	}
	void readDimline(istream& file) {
		file >> dim_;
		file >> type_;
		file >> mode_;
		subtype_ = 0;
	}

	bool operator==(const BasisParameters& par)const {
		if (par0_ != par.par0_) { return false; }
		if (par1_ != par.par1_) { return false; }
		if (par2_ != par.par2_) { return false; }
		if (par3_ != par.par3_) { return false; }
		if (dim_ != par.dim_) { return false; }
		if (type_ != par.type_) { return false; }
		if (mode_ != par.mode_) { return false; }
		if (subtype_ != par.subtype_) { return false; }
		return true;
	}

	[[nodiscard]] double omega() const { return par0_; }
	[[nodiscard]] double r0() const { return par1_; }
	[[nodiscard]] double wfr0() const { return par2_; }
	[[nodiscard]] double wfomega() const { return par3_; }

	bool operator!=(const BasisParameters& b) const {
		return !(*this == b);
	}

	double par0_ = 1.;
	double par1_ = 0.;
	double par2_ = 0.;
	double par3_ = 1.;

	size_t dim_ = 3;
	size_t type_ = 0;
	size_t mode_ = 0;
	size_t subtype_ = 0;
};

#endif //BASISPARAMETERS_H
