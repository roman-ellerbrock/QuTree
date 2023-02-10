//
// Created by Roman Ellerbrock on 4/14/20.
//
#include "FullRank.h"

namespace FullRank {

	void StatusTime(size_t it, size_t max, size_t freq, size_t length,
		chrono::high_resolution_clock::time_point& t1, chrono::high_resolution_clock::time_point& t2,
		chrono::microseconds& time) {

		double ptg = it * 1. / ((double) max);
		auto num = (size_t) (ptg * length);
		if (!(it % freq)) {

			string clear
				("                                                                                                     ");
			cout << "\r" << clear;
			cout.flush();
			t2 = chrono::high_resolution_clock::now();
			time += chrono::duration_cast<chrono::microseconds>(t2 - t1);
			t1 = t2;
			string msg = "Total time: " + to_string(time.count() / 1000000.) + "s\t";
			msg.append("Time per operation: " + to_string(time.count() / (it * 1000.)) + "ms\t");
			msg.append("Progress : [");
			for (size_t i = 0; i < num; ++i) {
				msg.append("#");
			}
			for (size_t i = num; i < length - 1; ++i) {
				msg.append(" ");
			}
			msg.append("]");
			cout << "\r" << msg;
			cout.flush();
		}
	}

	void Status(size_t it, size_t max, size_t freq, size_t length) {
		double ptg = it * 1. / ((double) max);
		auto num = (size_t) (ptg * length);
		if (!(it % freq)) {
			string msg = "Progress : [";
			for (size_t i = 0; i < num; ++i) {
				msg.append("#");
			}
			for (size_t i = num; i < length - 1; ++i) {
				msg.append(" ");
			}
			msg.append("]");
			cout << "\r" << msg;
			cout.flush();
		}
	}

	Wavefunction to_Wavefunction(const Wavefunction& Psi, const Tree& tree) {
		assert(0);
		return Psi;
	}

	Matrixcd toMatrix(const LeafOperatorcd& h, size_t mode, const Tree& tree) {
		TensorShape shape{2, 2};
		Tensorcd v(shape);
		Tensorcd hv(shape);
		v(0, 0) = 1.;
		v(1, 1) = 1.;
		const Leaf& leaf = tree.getLeaf(mode);
		h.apply(leaf.interface(), hv, v);
		Matrixcd mat = v.dotProduct(hv);
		return mat;
		//return mat.Adjoint().Transpose(); // TODO: this is a hack
	}

	Wavefunction ApplyOperator(Wavefunction Psi, const MLOcd& M, const Tree& tree) {
		for (size_t k = 0; k < M.size(); ++k) {
			const LeafOperatorcd& h = *M[k];
			const size_t mode = M.mode(k);
			auto mat = toMatrix(h, mode, tree);
			Psi = matrixTensor(mat, Psi, mode);
		}
		return Psi;
	}

	Wavefunction ApplyOperator(const Wavefunction& Psi, const SOPcd& sop, const Tree& tree) {
		Wavefunction sopPsi(Psi.shape());
		for (size_t i = 0; i < sop.size(); ++i) {
			const MLOcd& M = sop[i];
			const complex<double> coeff = sop.coeff(i);
			auto MPsi = ApplyOperator(Psi, M, tree);
			MPsi *= coeff;
			sopPsi += MPsi;
		}
		return sopPsi;
	}

	Wavefunction applyOperator(Wavefunction Psi, const SOPVectorcd& S, const Tree& tree) {
		size_t num = S.size();
		size_t idx = 0;
//		cout << "depth: " << S.size() << endl;
		for (const SOPcd& sop : S) {
			Status(idx++, num, 1, 20);
			Psi = ApplyOperator(Psi, sop, tree);
		}
		return Psi;
	}

	Wavefunction initialize(const Tree& tree) {
		vector<size_t> shape;
		for (size_t k = 0; k < tree.nLeaves(); ++k) {
			const Leaf& leaf = tree.getLeaf(k);
			shape.emplace_back(leaf.dim());
		}
		TensorShape x(shape);
		Wavefunction Psi(x);
		Psi(0) = 1.;
		return Psi;
	}

	void print(const Wavefunction& Psi) {
		const TensorShape& shape = Psi.shape();
		for (size_t k = 0; k < shape.order(); ++k) {
			auto rho = contraction(Psi, Psi, k);
			cout << "Qubit: " << k << endl;
			rho.print();
		}
	}

	double probability(const Wavefunction& Psi, size_t idx) {
		assert(idx <= Psi.shape().totalDimension());
		return pow(abs(Psi(idx)), 2);
	}

	void normalize(Wavefunction& Psi) {
		double norm = 0;
		for (size_t i = 0; i < Psi.shape().totalDimension(); ++i) {
			norm += pow(abs(Psi(i)), 2);
		}
		norm = sqrt(norm);
		Psi /= norm;
	}
}