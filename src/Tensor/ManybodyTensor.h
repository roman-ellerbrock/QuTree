#include "Tensor/Tensor.h"
#include "Tensor/TensorBLAS2.h"
#include "Tensor/TensorLapack.h"

TensorShape shapeStart(const TensorShape& target_shape, const TensorShape& large_shape) {
    vector<size_t> ns;
    for (size_t k = 0; k < target_shape.order(); ++k) {
        ns.push_back(large_shape[k]-target_shape[k]);
    }
    return TensorShape(ns);
}

void shiftIdxs(vector<size_t>& idx, const TensorShape& shift) {
    for (size_t k = 0; k < shift.order(); ++k) {
        idx[k] += shift[k];
    }
}

TensorShape delta(const TensorShape& start, const TensorShape& end) {
    vector<size_t> Delta;
    for (size_t k = 0; k < end.order(); ++k) {
        Delta.push_back(end[k] - start[k]);
    }
    return TensorShape(Delta);
}

Tensord slice(const Tensord& A, const TensorShape& start, const TensorShape& end) {
    TensorShape Delta = delta(start, end);
    Tensord S(Delta);
    vector<size_t> idx(Delta.order());
    for (size_t I = 0; I < Delta.totalDimension(); ++I) {
        indexMapping(idx, I, Delta);
        shiftIdxs(idx, start);
        size_t J = indexMapping(idx, A.shape_);
        S(I) = A(J);
    }
    return S;
}

/*Tensord slice(const Tensord& A, const TensorShape& shape_slice, const TensorShape& shift = TensorShape()) {
    Tensord S(shape_slice);
    TensorShape Delta = delta(A.shape_, shape_slice);

    vector<size_t> idx(shape_slice.order());
    for (size_t I = 0; I < shape_slice.totalDimension(); ++I) {
        indexMapping(idx, I, shape_slice);
        for (auto x : idx) {
            cout << x << " ";
        }cout << " - ";
        shiftIdxs(idx, Delta);
        shiftIdxs(idx, shift);
        size_t J = indexMapping(idx, A.shape_);
        for (auto x : idx) {
            cout << x << " ";
        }cout << "\n";
        S(I) = A(J);
    }
    return S;
}*/

class ManybodyTensor {
public:

    ManybodyTensor() = default;
    ~ManybodyTensor() = default;

    void hosvd(const Tensord& A) {
        C_ = A;
        U_.clear();
        for (size_t k = 0; k < A.shape_.order(); ++k) {
            auto rho = contraction(A, A, k);
            auto x = heev(rho);
            U_.push_back(x.U());
            C_ = matrixTensor(adjoint(U_[k]), C_, k);
        }
        shape_ = A.shape_;
    }

    Tensord recompose() const {
        Tensord A = C_;
        for (size_t k = 0; k < U_.size(); ++k) {
            A = matrixTensor(U_[k], A, k);
        }

        for (size_t k = 0; k < U1_.size(); ++k) {
            Tensord Ak = C1_[k];
            Ak = matrixTensor(U1_[k], Ak, k);
            for (size_t l = 0; l < U_.size(); ++l) {
                if (k == l) { continue; }
                Ak = matrixTensor(U_[l], Ak, l);
            }
            A += Ak;
        }

        return A;
    }

    void slicing(const TensorShape& n, size_t order = 1) {
        const TensorShape N = C_.shape_;
        const TensorShape Delta = delta(n, N);
        
        /// C.shape_ = (N, N, N)
        /// e.g. shape = (n, n, n)
        /// Delta = (N-n, N-n, N-n)
        Tensord A = C_;
        vector<Tensord> Us = U_;

        C_ = slice(A, Delta, N);
        U_.clear();
        for (size_t k = 0; k < C_.shape_.order(); ++k) {
            TensorShape start({0, Delta[k]});
            TensorShape end({N[k], N[k]});
            U_.push_back(slice(Us[k], start, end));
        }

        /// first order correction
        C1_.clear();
        U1_.clear();
        if (order == 0) { return; }
        for (size_t k = 0; k < C_.shape_.order(); ++k) {
            {
                TensorShape start = Delta;
                start.setDimension(0, k);
                TensorShape end = N;
                end.setDimension(Delta[k], k);
                C1_.push_back(slice(A, start, end));
            }

            {
                TensorShape start({0, 0});
                TensorShape end({N[k], Delta[k]});
                U1_.push_back(slice(Us[k], start, end));
            }
        }
    }

    size_t size() const {
        size_t s = C_.shape_.totalDimension();
        for (const auto& u : U_) {
            s += u.shape_.totalDimension();
        }
        for (const auto& c1 : C1_) {
            s += c1.shape_.totalDimension();
        }
        for (const auto& u1 : U1_) {
            s += u1.shape_.totalDimension();
        }
        return s;
    }

    Tensord C_; // core tensor
    vector<Matrixd> U_; // transformatinon to core tensor

    vector<Tensord> C1_; // first order correction slices
    vector<Matrixd> U1_; // correction matrix transformation

    TensorShape shape_;
};


void contraction(Matrixd& S, const ManybodyTensor& A, const ManybodyTensor& B, size_t k);