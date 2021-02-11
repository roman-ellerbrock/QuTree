//
// Created by Grace Johnson on 1/30/20.
//
// Tensor creation, population, basic operations, dot products, hole products
//

#include "Core/Tensor.h"
#include "Core/Matrix.h"

// Demonstrate various ways to create a Tensor object
Tensorcd create_tensor() {
    cout << "\ncreate_tensor:\n" << endl;
    // 1. Create from a TensorDim object
    TensorShape tdim({2, 3, 4, 2});
    Tensorcd A(tdim); // create tensor, all entries set to Zero
    cout << "A=" << endl;
    A.print();
	A.shape().print();

    // 2. Create from another Tensor
    Tensorcd B(A);

    // 3. Copy/move from another Tensor
    Tensorcd C = A;

    // 4. Read from an istream or file
    A.Write("tensor_A.dat");
    Tensorcd D("tensor_A.dat");

    return A;
}

Tensorcd fill_tensor(Tensorcd A){
    cout << "\nfill_tensor:\n" << endl;
    // Loop through complete super-index list
    cout << A.shape().totalDimension() << endl;
    for (size_t i = 0; i < A.shape().totalDimension(); i++) {
        A(i) = i * 0.01;
    }
    cout << "A=" << endl;
    A.print();

    // Loop through ntensor, then lower-index list
    Tensorcd B(A.shape());
    for (size_t n = 0; n < A.shape().lastDimension(); n++) {
        for (size_t j = 0; j < A.shape().lastBefore(); j++) {
            B(j, n) = (n * j + j * 2) * 0.002;
        }
    }
    cout << "B=" << endl;
    B.print();

    // Note: arithmetic operations work element-wise
    Tensorcd C = A - B;
    cout << "C=" << endl;
    C.print();

    return A;
}

Matrixcd dot_product(Tensorcd A, Tensorcd B) {
    cout << "\ndot_product:\n" << endl;
    Matrixcd w = A.DotProduct(B);
    cout << "w=" << endl;
    w.print();
    cout << "(" << w.dim1() << "," << w.dim2() << ")\n";
    return w;
}

void hole_product(Tensorcd A, const Tensorcd& B) {
    cout << "\nhole_product:\n" << endl;
    for (size_t k = 0; k < A.shape().order(); k++) {
        Matrixcd h = Contraction(A, B, k);
        cout << "\nk = " << k << ":\n h =" << endl;
        h.print();
        cout << "trace: " << h.trace() << endl;
        // Multiply
		Tensorcd C = MatrixTensor(h, A, k);
	}
}

void reshape (Tensorcd A) {
    cout << "\nreshape:\n" << endl;
    cout << "A.shape() = " << endl;
	A.shape().print();
    TensorShape tdim({4, 3, 2, 2});
    cout << "d = " << tdim.lastBefore() << endl;
    A.Reshape(tdim);
    cout << "A.shape() = " << endl;
	A.shape().print();
	assert(tdim == A.shape());
}

int main() {
    Tensorcd A = create_tensor();
    A = fill_tensor(A);
    Tensorcd B = A;
    Matrixcd w = dot_product(A, B);
    hole_product(A, B);
    reshape(A);
}
