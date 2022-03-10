//
// Created by Roman Ellerbrock on 2/19/20.
//

#ifndef SUMOFPRODUCTSOPERATOR_IMPLEMENTATION_H
#define SUMOFPRODUCTSOPERATOR_IMPLEMENTATION_H
#include "Operator/SumOfProductsOperator.h"

template<typename T>
SumOfProductsOperator<T>::SumOfProductsOperator(const PO<T>& M, T c) {
	push_back(M, c);
}

template<typename T>
SumOfProductsOperator<T> multAB(const SOP<T>& A, const SOP<T>& B) {
	SOP<T> C;

	for (size_t i = 0; i < A.size(); i++) {
		const PO<T> Ai = A(i);
		auto acoeff = A.coeff(i);
		for (size_t j = 0; j < B.size(); j++) {
			const PO<T> Bj = B(j);
			auto bcoeff = B.coeff(j);
			auto ccoeff = acoeff * bcoeff;
			C.push_back(Ai * Bj, ccoeff);
		}
	}

	return C;
}

template<typename T>
SumOfProductsOperator<T> operator*(T c, const SOP<T>& A) {
	SOP<T> C;
	for (size_t i = 0; i < A.size(); i++) {
		C.push_back(A(i), c * A.coeff(i));
	}
	return C;
}

template<typename T>
SumOfProductsOperator<T> operator*(const SOP<T>& A, T c) {
	SOP<T> C;
	for (size_t i = 0; i < A.size(); i++) {
		C.push_back(A(i), A.coeff(i) * c);
	}
	return C;
}

template<typename T>
SumOfProductsOperator<T> operator*(const PO<T>& M, const SOP<T>& A) {
	SOP<T> C;
	for (size_t i = 0; i < A.size(); i++) {
		PO<T> MA = M * A(i);
		C.push_back(MA, A.coeff(i));
	}
	return C;
}

template<typename T>
SumOfProductsOperator<T> operator*(const SOP<T>& A, const PO<T>& M) {
	SOP<T> C;
	for (size_t i = 0; i < A.size(); i++) {
		PO<T> MA = A(i) * M;
		C.push_back(MA, A.coeff(i));
	}
	return C;
}

template<typename T>
SumOfProductsOperator<T> operator*(const SOP<T>& A,
	const SOP<T>& B) {
	SOP<T> C;
	for (size_t i = 0; i < A.size(); i++) {
		const PO<T> Ai = A(i);
		auto acoeff = A.coeff(i);
		for (size_t j = 0; j < B.size(); j++) {
			const PO<T> Bj = B(j);
			auto bcoeff = B.coeff(j);
			auto ccoeff = acoeff * bcoeff;
			C.push_back(Ai * Bj, ccoeff);
		}
	}
	return C;
}

// add SoP-Operator
template<typename T>
SumOfProductsOperator<T> operator+(const SumOfProductsOperator<T>& A,
		const SumOfProductsOperator<T>& B) {
	SOP<T> C = B;

	for (size_t i = 0; i < A.size(); i++) {
		C.push_back(A(i), A.coeff(i));
	}

	return C;
}

#endif //SUMOFPRODUCTSOPERATOR_IMPLEMENTATION_H
