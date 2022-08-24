//
// Created by Roman Ellerbrock on 11/22/21.
//

#ifndef BASISAPI_H
#define BASISAPI_H
#include "PrimitiveBasis.h"
#include "DVR.h"
#include "FFTGrid.h"
#include "HarmonicOscillator.h"
#include "LegendrePolynomials.h"
#include "SpinGroup.h"

class BasisAPI {
public:
	BasisAPI() : ptr_(nullptr) {}
	~BasisAPI() {
		delete ptr_;
	}

	BasisAPI(const BasisAPI& prev) : BasisAPI() {
		initialize(prev.ptr_->par_);
	}

	BasisAPI& operator=(const BasisAPI& prev) {
		if (&prev == this) {
			return *this;
		}
		initialize(prev.ptr_->par_);
		return *this;
	}

	BasisAPI(BasisAPI&&) = delete;
	BasisAPI& operator=(BasisAPI&&) = delete;

	void initialize(const BasisParameters& par) {

		delete ptr_;
		size_t type = par.type_;
		if (type == 0) {
			ptr_ = new HarmonicOscillator;
		} else if (type == 1) {
			ptr_ = new FFTGrid;
		} else if (type == 2) {
			ptr_ = new LegendrePolynomials;
		} else if (type == 6) {
			ptr_ = new SpinGroup;
		} else {
			cout << "Error: This Basis Type is not in the known list of "
			<< "Typs. The iplemented ones are: \n"
			<< "0 = Hermite-DVR\n"
			<< "1 = FFT-Grid\n"
			<< "2 = Legendre-DVR\n"
			<< "6 = Spin Group\n";
			exit(1);
		}
		ptr_->initialize(par);
	}

	void write(ostream& os) const {
		ptr_->write(os);
	}

	void applyX(Tensorcd& xA, const Tensorcd& A) const {
		if (!ptr_) {
			cerr << "Primitive Basis not initialized in BasisAPI.\n";
			exit(1);
		}
		ptr_->applyX(xA, A);
	}

	void applyX2(Tensorcd& xA, const Tensorcd& A) const {
		if (!ptr_) {
			cerr << "Primitive Basis not initialized in BasisAPI.\n";
			exit(1);
		}
		ptr_->applyX2(xA, A);
	}

	void applyP(Tensorcd& xA, const Tensorcd& A) const {
		if (!ptr_) {
			cerr << "Primitive Basis not initialized in BasisAPI.\n";
			exit(1);
		}
		ptr_->applyP(xA, A);
	}

	void applyKin(Tensorcd& xA, const Tensorcd& A) const {
		if (!ptr_) {
			cerr << "Primitive Basis not initialized in BasisAPI.\n";
			exit(1);
		}
		ptr_->applyKin(xA, A);
	}

	[[nodiscard]] const PrimitiveBasis* ptr() const { return ptr_; };
	[[nodiscard]] PrimitiveBasis* ptr() { return ptr_; };

protected:
	PrimitiveBasis* ptr_;
};


#endif //BASISAPI_H
