//
// Created by Roman Ellerbrock on 11/22/21.
//

#ifndef LEAFAPI_H
#define LEAFAPI_H
#include "PrimitiveBasis.h"
#include "DVR.h"
#include "FFTGrid.h"
#include "HarmonicOscillator.h"
#include "LegendrePolynomials.h"
#include "SpinGroup.h"

class LeafAPI {
public:
	LeafAPI() : basis_(nullptr) {}
	~LeafAPI() {
		delete basis_;
	}

	LeafAPI(const LeafAPI& prev) : LeafAPI() {
		initialize(prev.basis_->par_);
	}

	LeafAPI& operator=(const LeafAPI& prev) {
		if (&prev == this) {
			return *this;
		}
		initialize(prev.basis_->par_);
		return *this;
	}

	LeafAPI(LeafAPI&&) = delete;
	LeafAPI& operator=(LeafAPI&&) = delete;

	void initialize(const BasisParameters& par) {

		delete basis_;
		size_t type = par.type_;
		if (type == 0) {
			basis_ = new HarmonicOscillator;
		} else if (type == 1) {
			basis_ = new FFTGrid;
		} else if (type == 2) {
			basis_ = new LegendrePolynomials;
		} else if (type == 6) {
			basis_ = new SpinGroup;
		} else {
			cout << "Error: This Basis Type is not in the known list of "
			<< "Typs. The iplemented ones are: \n"
			<< "0 = Hermite-DVR\n"
			<< "1 = FFT-Grid\n"
			<< "2 = Legendre-DVR\n"
			<< "6 = Spin Group\n";
			exit(1);
		}
		basis_->initialize(par);
	}

	void write(ostream& os) const {
		basis_->write(os);
	}

	[[nodiscard]] const PrimitiveBasis* basis() const { return basis_; };

protected:
	PrimitiveBasis* basis_;
};


#endif //LEAFAPI_H
