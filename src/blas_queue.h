
#ifndef BLASQUEUE_H
#define BLASQUEUE_H
#include <blas.hh>

namespace qutree {
	static int device{0};
	static int batch{1000};
	static blas::Queue queue(device, batch);
}

#endif // BLASQUEUE_H
