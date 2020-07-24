====================
Example 5: Operators
====================

QuTree uses abstract operators using class and function representations.
A common representation is the sum-of-product form.
The most fundamental operator is the :code:`LeafOperator` which acts on
a single leaf only. :code:`MultiLeafOperator` represent products of
:code:`LeafOperator` objects. Last, the :code:`SumOfProductOperator` is
a sum of :code:`MultiLeafOperator`.

Leaf Operators
==============

Leaf operators are the building block for higher level operators.
They defined by hand, for example by creating classes that
specify how the operator acts, which memory it uses and so on; alternatively,
the LeafInterface provides some operators whiches implementation depends on the
chosen basis. Matrices can also always be used as LeafOperators.
Here are some examples of how to construct LeafOperators

.. code-block:: C++

    // Fundamental operators from LeafInterface
    LeafOperatorcd& x = &LeafInterface::applyX;
    LeafOperatorcd& x2 = &LeafInterface::applyX2;
    LeafOperatorcd& p = &LeafInterface::applyP;
    LeafOperatorcd& kin = &LeafInterface::applyKin;

    // Create a LeafMatrix from a matrix
    Matrixcd M(2, 2);
    LeafMatrixcd lm(M);

Multi-Leaf Operators (MLO or Product Operator)
==============================================

The :code:`MultiLeafOperator` is the next level of operators. They are products
of :code:`LeafOperator` objects that act on specific leafs. They are stored as
a vector of Leaf operators and corresponding indices that tell which leaf the
operator is applied to. Therefore a :code:`push_back` can be used to multiply
operators from the left.

.. code-block:: C++

    MLOcd M;
    M.push_back(x, 0);
    M.push_back(x, 1);
    M.push_back(x2, 2);

This creates a product operator :code:`x(0)*x(1)*x^2(2)`.
Product operators can be multiplied

.. code-block:: C++

    MLOcd N(x, 0);
    MLOcd O = M * N;

Sum-of-Product operators (SOP)
==============================

Sum-of-Product Operators (SOP) can be constructed by summing product operators
and use operator overloading

.. code-block:: C++

    SOPcd S;
    S.push_back(M, 1.);
    S.push_back(N, 1.);
    SOPcd twoS = S + S;
