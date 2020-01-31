//
// Created by Grace Johnson on 1/30/20.
//
// Tensor creation, population, basic operations, dot products, hole products
//

#include "Core/Tensor.h"

// Demonstrate various ways to create a Tensor object
Tensorcd create_tensor() {
    cout << "\ncreate_tensor:\n" << endl;
    // 1. Create from a TensorDim object
    vector<size_t> dims = {2, 3, 4}; // third order tensor
    size_t ntensor = 2; // group of two tensors
    TensorDim tdim(dims, ntensor);
    Tensorcd A(tdim); // create tensor, all entries set to Zero
    cout << "A=" << endl;
    A.print();
    A.Dim().print();

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
    cout << A.Dim().GetDimTot() << endl;
    for (size_t i = 0; i < A.Dim().GetDimTot(); i++) {
        A(i) = i * 0.01;
    }
    cout << "A=" << endl;
    A.print();

    // Loop through ntensor, then lower-index list
    Tensorcd B(A.Dim());
    for (size_t n = 0; n < A.Dim().GetNumTensor(); n++) {
        for (size_t j = 0; j < A.Dim().GetDimPart(); j++) {
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
    cout << "(" << w.Dim1() << "," << w.Dim2() << ")\n";
    return w;
}

Matrixcd hole_product(Tensorcd A, Tensorcd B) {
    cout << "\nhole_product:\n" << endl;
    Matrixcd h;
    for (size_t k = 0; k < A.Dim().GetOrder(); k++) { // TODO: what is GetOrder?
        h = HoleProduct(A, B, k);
        cout << "\nk = " << k << ":\n h =" << endl;
        h.print();
        cout << "Trace: " << h.Trace() << endl;
        // TODO: do we have machinery for B = hk * A ? Use FactorMatrixcd
    }
    return h;
}

void reshape (Tensorcd A) {
    cout << "\nreshape:\n" << endl;
    cout << "A.Dim() = " << endl;
    A.Dim().print();
    vector<size_t> dims = {2, 3, 2, 2}; // fourth order order tensor
    size_t ntensor = 2; // group of two tensors
    TensorDim tdim(dims, ntensor);
    assert(tdim.GetDimPart() == A.Dim().GetDimPart());
    cout << "d = " << tdim.GetDimPart() << endl;
    A.Reshape(tdim);
    cout << "A.Dim() = " << endl;
    A.Dim().print();
}

int main() {
    Tensorcd A = create_tensor();
    A = fill_tensor(A);
    Tensorcd B = A;
    Matrixcd w = dot_product(A, B);
    Matrixcd h = hole_product(A, B);
    reshape(A);
}