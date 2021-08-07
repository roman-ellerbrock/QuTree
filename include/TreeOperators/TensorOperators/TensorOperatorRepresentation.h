//
// Created by Roman Ellerbrock on 8/4/21.
//

#ifndef TENSOROPERATORREPRESENTATION_H
#define TENSOROPERATORREPRESENTATION_H

class TensorOperatorRepresentation : public TensorTreed {
public:
	TensorOperatorRepresentation() = default;
	TensorOperatorRepresentation(const Tree& tree, const SOPd& S) {
		size_t nparts = S.size();
	}
	~TensorOperatorRepresentation() = default;
};

#endif //TENSOROPERATORREPRESENTATION_H
