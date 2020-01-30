#include "Core/TensorDim_Extension.h"

namespace TensorDim_Extension {
// Adjust dimensions of TensorDim
TensorDim ReplaceActive(const TensorDim& tdim, size_t mode, size_t new_dim) {

	// Replace the active dim in mode
	vector<size_t> dimlist = tdim.getdimlist();
	assert(mode < dimlist.size());
	dimlist[mode] = new_dim;

	return TensorDim(dimlist, tdim.getntensor());
}
TensorDim ReplaceNtensor(const TensorDim& tdim, size_t ntensor) {
	// Replace ntensor in tensordim
	vector<size_t> dimlist = tdim.getdimlist();
	return TensorDim(dimlist, ntensor);
}

}
