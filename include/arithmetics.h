#include "GraphFactory.h"
#include "TensorNetwork.h"

namespace qutree {
void dotProduct(TensorNetwork &mt, const TensorNetwork &Bra,
                const TensorNetwork &Ket);

void contract(TensorNetwork &mt, const TensorNetwork &Bra,
                const TensorNetwork &Ket);

TensorNetwork qr(TensorNetwork tn);

} // namespace qutree
