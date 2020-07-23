=======================
Example 4: TensorTree
=======================

These examples will demonstrate usage of the the :code:`TensorTree` class.
It is not only one of the most central classes of the package but also
a good example for an inherited :code:`NodeAttribute`, which is a base class
for objects that have an object attributed to every node in a tree.
The examples will also cover MatrixTree which have a Matrix associated to every
node.

TensorTree Contruction
======================

A :code:`TensorTree` can be constructed using a :code:`Tree`::

    // Create a TensorTree filled with zero coefficients
    Tree tree = TreeFactory::BalancedTree(8, 2, 2);
    TensorTreecd Psi(tree);

    // Create a randomly occupied tensor tree
    mt19937 gen; // RNG
    TensorTreecd Chi(gen, tree);
    Chi.print(tree);

Tensor Tree Handling
====================

The tensors in the :code:`TensorTree` can be directly addressed using
a swipe over the tree::

    for (const Node& node : tree) {
        const Tensorcd& Phi = Psi[node];
        auto x = Phi.DotProduct(Phi);
    }

This way of addressing objects associated to a node can be used for any class
that inherits from :code:`NodeAttribute`, e.g. :code:`MatrixTree` as well.

Dot-Products and Contractions
=============================

Similar to tensors, TensorTrees often represent vectors in high-dimensional
vector spaces. Dot-products between :code:`TensorTree` objects can be calculated using::

    MatrixTree S = DotProduct(Chi, Psi, tree);

Contractions of TensorTrees can be caluclated using::

    MatrixTree q = Contraction(Chi, Psi, S, tree);

If the TensorTree uses an orthonormal basis, the self-contraction can be performed more
efficiently::

    MatrixTree rho = Contraction(Psi, tree, true);

where the flag "true" indicated whether the basis is orthogonal or not.

Note that in the context of quantum dynamics applications, a TensorTree often represents
a wavefunction, the dot-product represent the overlap of two wavefunctions and
a Contraction is used to calculate reduced density matrices.
