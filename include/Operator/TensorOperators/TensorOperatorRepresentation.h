//
// Created by Roman Ellerbrock on 8/4/21.
//

#ifndef TENSOROPERATORREPRESENTATION_H
#define TENSOROPERATORREPRESENTATION_H

template <typename T>
class TensorOperatorRepresentation : public TensorTree<T> {
public:
	TensorOperatorRepresentation() = default;
	TensorOperatorRepresentation(const Tree& tree, const SOP<T>& S) {
		size_t nparts = S.size();
	}
	~TensorOperatorRepresentation() = default;
};

typedef TensorOperatorRepresentation<double> TensorOperatorRepresentationd;
typedef TensorOperatorRepresentation<complex<double>> TensorOperatorRepresentationcd;

#endif //TENSOROPERATORREPRESENTATION_H
