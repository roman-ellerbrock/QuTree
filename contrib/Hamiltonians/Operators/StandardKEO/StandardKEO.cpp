#include "StandardKEO.h"


StandardKEO::StandardKEO(const mctdhBasis& basis, size_t dim)
{
	SpecialInitialize(basis, dim);
}

StandardKEO::~StandardKEO()
{
}

void StandardKEO::SpecialInitialize(const mctdhBasis& basis, size_t dim)
{
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> kin = &PrimitiveBasis::ApplyKin;
	cout << "dof for Standard KEO: " << dim << endl;

	for (size_t l = 0; l < dim; l++) {
		MultiParticleOperator M;
		M.push_back(kin, l);
		push_back(M, 1.);
	}
}

