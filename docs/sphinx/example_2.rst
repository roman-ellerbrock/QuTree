===================================
Example 2: Matrix Class
===================================

The :code:`Matrix` class represents Matrices on and
is, on purpose, separated from the Tensors class, since it
requires less overhead and typically has usages.

Fundamental Usage
=================

The following examples demonstrate how matrices are constructed
and show some basic usage examples.

    Matrixcd A(3, 3);
    A(0, 0) = 2;
    A(4) = 5;
    Matrixcd B = A;
    Matrixcd C = 0.5 * (A + B) * A.Transpose();
    Matrixcd D = 0.5 * (A + A.Adjoint());

Spectral Decompositions
=======================

Hermitian (or symmetric) matrices can be diagonalized using::

    SpectralDecompositioncd spec = Diagonalize(D);
    Matrixcd U = spec.first;
    Vectord ev = spec.second;
    // Rebuild matrix from its spectral decomposition
    Matrixcd Dnew = toMatrix(spec);
    auto Dreg = Regularize(D, 1e-6);

The decomposed matrix can be used to calculate squares, the inverse
and so on::

    Matrixcd Dsq = sqrt(spec);
    Matrixcd Dinv = BuildInverse(spec);

A singular value decomposition can be performed by::

    SVDcd svd_a = svd(A);
    Matrixcd Anew = toMatrix(svd_a);

Other high-level operations can be performed in a similar fashion, e.g. solving SLE or
QR-decomposition.

Eigen Interface
===============

QuTree provides an interface for :code:`Matrix` to Eigen's matrix class::

    Eigen::MatrixXcd Aeigen = toEigen(A);
    Matrixcd Aqutree = toQutree(Aeigen);


Further Functionality
====================

There are various other functions provided by the :code:`Matrix` class.
Please refer to the source code and unit tests for deeper insight.
If a high-level function is not provided, you can most likely use the
Eigen routines via the interface.
