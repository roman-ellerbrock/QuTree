//
// Created by thomas on 04.08.18.
//

#ifndef MCTDH_LANZCOS_H
#define MCTDH_LANZCOS_H

#include <glob.h>
#pragma once
#include "Core/TensorShape.h"
#include "Core/Tensor.h"
#include "Core/Matrix.h"
#include "Core/Vector.h"

template<class A, typename k>
class Lanzcos {
public:
	Lanzcos()
		: diagonal(false) {}

	//Creates a Krylov space of Order N
	void buildKrylovSpace(const Tensor<k>& v, size_t N,
		double parameter, A Operator,
		function<void(A&, const double, Tensor<k>&,
			const Tensor<k>&)> apply);

	//Diagonalizes the krylov matrix
	void diagonalizeKrylovSpace();

	//applys functions of the form f(parameter * Oerator) on
	//the states the krylov space is created with.
	Tensorcd applyfunction(double parameter,
		function<double(double)> fun) const;

	Tensorcd getEigenstates();
	vector<double> getEigenvalues();

private:

	void checkBuild() const;
	void checkDiagonal() const;

	vector<Matrixcd> a_; // diagonal
	vector<Matrixcd> b_; // off diagonal
	vector<Tensor<k>> kV_; // Krylov vector

	//Information after diagonalisation
	Matrixcd transformation_;
	Vectord eigenvalues_;

	bool diagonal_;
};

template<class A, typename k>
void Lanzcos<A, k>::buildKrylovSpace(const Tensor<k>& v, size_t N,
	double parameter, A Operator,
	function<void(A&, const double, Tensor<k>&, const Tensor<k>&)> apply) {
	//Initialize arrays
	a_.clear();
	b_.clear();
	kV_.clear();
	diagonal_ = false;

	//Define Storage for krylov Vectors
	for (size_t i = 0; i < N; ++i)
		kV_.emplace_back(v);

	//b[0] is zero
	auto states = v.Dim().getntensor();
	b_.emplace_back(Matrixcd(states, states));

	for (auto s = 0; s < N - 1; s++) {
		//calculate H Psi
		Operator.ApplyH(parameter, kV_[s + 1], kV_[s]);

		//calculate <H>^(s)
		a_.push_back(kV_[s].DotProduct(kV_[s + 1]));

		//calculate next state space
		kV_[s + 1] -= multStateAB(a_[s], kV_[s]);
		if (s > 0) {
			kV_[s + 1] -= multStateAB(b_[s], kV_[s - 1]);
		}

		//save for calculating B later
		Tensor<k> copy(kV_[s + 1]);

		//orthonormalize the state space
		GramSchmidt(kV_[s + 1]);

		//calculate normalisation matrix
		//b.push_back(copy.DotProduct(kV[s+1]));
		b_.push_back(kV_[s + 1].DotProduct(copy));
	}
}

template<class A, typename k>
void Lanzcos<A, k>::diagonalizeKrylovSpace() {
	checkBuild();

	size_t N = a.size();
	size_t states = a_[0].Dim1();

	//Create Krylovmatrix
	Matrixcd krylovmatrix(N * states, N * states);

	//Map blocks on big matrix
	for (size_t s = 0; s < N; s++)
		for (size_t i = 0; i < states; i++)
			for (size_t j = 0; j < states; j++) {
				krylovmatrix(s * states + i, s * states + j) = a_[s](i, j);
				if (s < N - 1) {
					krylovmatrix(s * states + i, (s + 1) * states + j)
						= b_[s + 1](j, i);
					krylovmatrix((s + 1) * states + i, s * states + j)
						= conj(b_[s + 1](i, j));
				}
			}

	//Diagonalize Krylov matrix
	Matrixcd transformation_(N * states, N * states);
	Vectord eigenvalues_(N * states);
	krylovmatrix.cDiag(transformation_, eigenvalues_);
	transformation = transformation_;
	eigenvalues = eigenvalues_;

	diagonal = true;
}

template<class A, typename k>
Tensorcd Lanzcos<A, k>::applyfunction(double parameter,
	function<double(
		double)> fun) const {
	checkBuild();

	return Tensorcd();
}

template<class A, typename k>
Tensorcd Lanzcos<A, k>::getEigenstates() {
	checkBuild();
	checkDiagonal();

	size_t N = a_.size();
	size_t dimpart = kV_[0].Dim().getdimpart();
	size_t states = kV_[0].Dim().getntensor();

	//create copy for result
	Tensorcd copy(kV_[0]);
	copy.Zero();

	//transform eigenvectors out of krylovspace
	for (size_t s = 0; s < N; ++s)
		for (size_t j = 0; j < states; ++j)
			for (size_t i = 0; i < states; ++i)
				for (size_t p = 0; p < dimpart; ++p) {
					copy(p, i) += kV_[s](p, j)
						* (transformation(s * states + j, i));
				}

	return copy;
}

template<class A, typename k>
vector<double> Lanzcos<A, k>::getEigenvalues() {
	checkBuild();
	checkDiagonal();

	auto states = a_[0].Dim1();

	vector<double> out;
	for (auto i = 0; i < states; ++i) {
		out.push_back(eigenvalues(i));
	}

	return out;
}

template<class A, typename k>
void Lanzcos<A, k>::checkBuild() const {
	if (a.size() == 0) {
		cout << "Error: Lanzcos Space was not build." << endl;
		exit(0);
	}
}

template<class A, typename k>
void Lanzcos<A, k>::checkDiagonal() const {
	if (!diagonal) {
		cout << "Error: Lanzcos matrix was not diagonalized." << endl;
		exit(0);
	}
}

#endif //MCTDH_LANZCOS_H
