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

    TensorShape tdim({2, 3, 4, 2});
    Tensorcd A(tdim); // create tensor, all entries set to Zero


Tensor Inspection
=================

Various elements of the 

Assigning Elements
==================


Dot Products
============


Hole Products
=============