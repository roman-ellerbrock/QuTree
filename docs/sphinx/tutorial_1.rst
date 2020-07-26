================================
Tutorial 1: Tensors
================================

The aim of this exercise is to gain insight into the TensorShape and Tensor classes
that come with QuTree.

Information
-----------
This programming tasks in this exercise exclusively use the
classes TensorShape and Tensor. The class declarations can be found in (relative
to project directory):

- QuTree/include/Core/TensorShape.h
- QuTree/include/Core/Tensor.h

Additional information about QuTree will be
central to solve this exercise. Please refer to Exercise 1.

Exercise 1
-----------

We want to create a Tensor and perform some operations with this Tensor. The exercises
below should offer an introduction and some ideas on how to work with Tensors. Feel free
to design your own exercises and use the QuTree classes to write your own tests.

The class Tensor is a template class that can be initialized using different base types.
For simplicity, you can use the base type complex<double>, i.e. the class Tensorcd.

1.) Create a TensorShape named shape with 4 lower indices. The upper index' dimension (ntensor) shall be 3. The lower indices' dimensions shall be 2, 3, 4 and 5.

2.) Create 2 Tensors, A and B, (type: double or complex<double>) with the dimensions given by shape. Initialize the coefficients of the Tensor using varying values.

3.) Calculate the dot-products (A, A), (A, B), (B, B) and print the result.

4.) Use the Gram-Schmidt method to orthonormalize Tensors A and B. Print the resulting Tensors. Then repeat 4.) using the orthonormal Tensors.

5.) Calculate all possible hole-products of Tensor A and B with itself, respectively. Print the results. Calculate and analyse the trace of the resulting matrices.

