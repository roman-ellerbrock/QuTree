==================
Example 1: Tensors
==================

This example shows manipulation of the :code:`Tensor` class as shown in :code:`examples/examples_1.cpp`.

The :code:`Tensor` is templated over data type, with two common typedefs provided::

    typedef Tensor<complex<double>> Tensorcd;
    typedef Tensor<double> Tensord;

Tensor Creation
===============
The most basic construction of a :code:`Tensor` is from a :code:`TensorShape`::

    TensorShape shape({2, 3, 4, 5});
    Tensorcd A(shape); // create tensor, all entries set to Zero

This will create a 4th order tensor with dimensions 2, 3, 4 and 5, respectively. By default,
memory is allocated and the coefficients are set to zero.

Assigning Elements and Inspection
=================

The elements of a :code:`Tensor` can be set and used via the bracket operators.The most common
bracket operator addresses the tensor elements as linear memory::

    for (size_t i = 0; i < A.shape().totalDimension(); ++i) {
        A(i) = 3. * i;
    }
    A.print();
    Tensorcd B = A;

The :code:`print()` function can be used to print the tensor entries in ASCII-format. You can
either pass a stream to print or it will print the output to :code:`std::cout` by default.
At the end, we create a second Tensor B with the same size and entries as A.

Tensor Contractions
===================

A central Tensor operation is the Tensor :code:`Contraction`. A tensor contraction of
two tensors A and B performs a summation over every index, except for a chosen index - the
so-called hole-index. We can perform a contraction on tensors A and B via::

    Matrixcd B = Contraction(A, B, 0);

This will return a Matrix with the dimension 2x2, since the shape of the tensor was
2,3,4 and 5 and we picked the first index to be excluded from the summation.

Dot Products
============

In many applications, we will store the expansion coefficients of a set of multidimensional
vectors via the Tensor class. Let's say we have 5 vectors on a direct product
vector space, each with the dimensions 2*3*4=24. We can then interpret the tensors A
and B as the numerical representation of such vector elements.
We can calculate the dot-product between the five vectors stored in A with the five vectors
stored in B using the Dot-product routine::

    Matrixcd S = A.DotProduct(B);
    S.print();

Here :code:`S` will be a 5x5 "overlap" matrix containing the dot-product values of a
standard (Euclidean) dot-product.
Note that the dot-product is equivalent to a tensor contraction where the last index
is excluded from the summation.

Matrix Tensor Products
======================

Matix Tensor Products can be performed via::

    Matrixcd M(3, 3);
    Tensorcd C = MatrixTensor(M, A, 1);

Here "1" indicates that the matrix is applied to the second index in A. The routine will check
whether the dimensions of M matches that active index in A upon calling the routine.
There also exist versions that allow to apply adjoint matrices (see :code:`multATB` function)
or that change the order of M and A which is equivalent to applying a transpose matrix M
to the Tensor A (see :code:`TensorMatrix(A, M, 1)`). Please refer to the source code for more
information.

Operators
=========

There exist various operator overloadings of tensors::

    Tensorcd D = C - B;
    D *= 0.5;
    D += C;

and so on.

Reshaping and Re-interpretation
===============================

Tensors and Matrices can be reinterpreted without copying memory, simply by changing
the dimensions. A commonly requested option is the Tensor reshaping::

    TensorShape newshape({5,2,12});
    A.Reshape(newshape);

will reinterpret the linear memory, previously interpreted according to the dimensions
2,3,4 and 5 now by 5,2,12. The reshaping only works if the total number of coefficients
it the same in both shapes. If the number of coefficients is not the same, the :code:`adjust`
function can be used to resize the tensor.

Creating Tensors from Externally Allocated Memory
=================================================

When using QuTree in existing frameworks the memory of tensors
can be previously allocated. At some point one may wish to interpret
a chunk of memory as a tensor from now on. This can be done by
passing the pointer to the tensor contructor::

    complex<double> *ptr = new complex<double>[120];
    bool zero = true;
    bool ownership = true;
    Tensorcd E(shape, ptr, ownership, zero);

By default, the ownership of the memory is transferred to the
tensor and the elemtns are set to zero.

More Examples
=============

There are various other functions of the :code:`Tensor` class. Giving a comprehensive overview
overview would inflate this example section a lot. Instead, we refer the reader to the source code
and the unit test section for more information. The :code:`Tensor.h` contains the
most fundamental routines and :code:`Tensor_Extension` contains more routines that
are sometimes used but not very frequently.