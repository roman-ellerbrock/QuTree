================================
Tutorial 1: Tensors and Matrices
================================

The aim of this exercise is to learn fundamental usage of QDlib.

Information
===========
This programming tasks in this exercise exclusively use the classes
TensorShape, Tensor and Matrix.
The class declarations can be found in (relative to project directory):

- QuTree/include/Core/TensorShape.h
- QuTree/include/Core/Tensor.h
- QuTree/include/Core/Matrix.h

Additional information is required solve the exercises. Use example 1 and 2
to get more information about the classes.

Exercise 2
==========

We want to create a Matrix and apply it to a Tensor.

The class Matrix, like the class Tensor, is a template class that
can be initialized using different base types.
For simplicitly you can use the base type complex<double>, i.e. the classes
Tensorcd and Matrixcd.

1.) Create a TensorsShape named *shape* with 5 indices. The upper (last) index' dimension
    shall be 3. The lower indices' dimensions shall be 2, 3, 4 and 5.

2.) Create 2 Tensors, A and B, with the dimensions given by *shape*. Initialize the
    coefficients of the Tensor using varying values and use a Gram-Schmidt to orthonormalize
    the entries.

3.) Create a Matrix named *mat* acting on an arbitrary index of A and B.
    Fill it with entries (create, e.g., a unit-operator, null-operator, ..)
    and apply it to A or B.
    Investigate the resulting Tensors C = mat * A.

4.) Calculate and investigate the dot-product (A, C) = (A, mat * A).

If you are interested to go even deeper into QuTree and understand the implementation of
the Tensor class, you can try to solve the following exercises:

5.) Write a function that performs a product of a SPOcd and a Tensorcd.

6.) Write a function that performs a Tensor-Hole product for a Tensorcd with a hole in
    index k.

The functions "MatrixTensor" and "Contraction" in Tensor_Implementation.h perform
the tasks 8.) and 9.), respectively. However, these functions are hard to read, because they
avoid the usage of bracket operators due to their overhead. In the present exercise you can use
all bracket operators. This significantly facilitates writing these functions.
The function "OldStateAveragedHoleProduct" in Tensor_Extension.h has an unoptimized version
of a hole-product that can be used as a reference. It is the intended solution for 9.).


