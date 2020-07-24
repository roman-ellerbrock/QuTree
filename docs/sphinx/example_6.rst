===================================
Example 6: Operator representations
===================================

Objects of the classes :code:`LeafOperator`, :code:`MultiLeafOperator` and
:code:`SumOfProductOperator` can be represented in the basis provided by a
:code:`TensorTree`. Corresponding to every active node in the :code:`Tree`,
the operator is represented by a Matrix. The representation is only required for
nodes where the representation is not trivial, i.e. not a identity matrix.
Active nodes are detected automatically by QuTree and the resulting object
belongs to the class :code:`SparseMatrixTree`.
Similar to the :code:`DotProduct` and :code:`Contraction` routines for TensorTrees,
there is a bottom-up and top-down building routines for Operators.

The bottom-up representation of an operator can be build via

.. code-block:: C++

    LeafOperatorcd& x2 = &LeafInterface::applyX2;
    LeafOperatorcd& p = &LeafInterface::applyP;
    MLOcd h(p, 0);
    h.push_back(x2, 2);

    Tree tree = TreeFactory::BalancedTree(8, 2, 2);
    mt19937 gen; // RNG
    TensorTreecd Psi(gen, tree);

    SparseMatrixTreecd Hmat = TreeFunctions::Represent(h, Psi, tree);

The top-down representation of an operator (sometimes called mean-field matrices)
requires previously build bottom-up matrices and is build by

.. code-block:: C++

    SparseMatrixTreecd Hmean = TreeFunctions::Contraction(h, Psi, Hmat, tree);

The operator representation of a single LeafOperator or MultiLeafOperator results
in a single :code:`SparseMatrixTree`; a :code:`SumOfProductOperator` representation
requires a :code:`SparseMatrixTree` for every summand, the list of those
objects is saved as a vector of SparseMatrixTrees, also called :code:`SparseMatrixTrees`.

.. code-block:: C++

    SOPcd H;
    H.push_back(h, 1.);
    H.push_back(h, 0.5);

    SparseMatrixTreescd Hmats = TreeFunctions::Represent(H, Psi, tree);
    SparseMatrixTreescd Hmeans = TreeFunctions::Contraction(H, Psi, Hmats, tree);

SparseTrees
===========

A :code:`SparseTree` is a directed subset of nodes in a tree. SparseTrees manage
