================================
Tutorial 3: Trees
================================

The aim of this exercise is to learn the fundamental usage of the Tree
representations in QuTree, learn about the node-handling and learn how
to access dimensions.

Information
===========
This programming tasks in this exercise exclusively use the classes Tree,
GetNode and TensorDim. The class declarations can be found in (relative to project directory):

- QuTree/TreeShape/Tree.h
- QuTree/TreeShape/Node.h
- QuTree/TreeShape/Leaf.h

Example 3 will be helpful for completing this tutorial.
The relationship between diagrammatric multilayer MCTDH basis representations [1] and
its corresponding C++ classes is shown in Fig.1.
Both, the diagram and the header-files can be useful for completing this exercise.
Thus, it is useful to take a look at them before and during completing the exercises.

Exercise 3
===========

We calcualte the number of coefficients are required to
represent the multilayer-wavefunction of different chemical and model systems.
The following basis sets will be investigated:

- nocl.basis
  + This basis can be used to calculate the photo-dissociation spectrum of NOCl[2].
- ch3.basis
  + This basis can be used to calculate the lowest 36 vibrational eigenstates of methyl[3].
- binary.constant.basis
  + This is a model of a binary multilayer tree with a constant number of SPFs when accending
    the tree. A standard model when investigating scaling in multilayer MCTDH.
- binary.constant.basis
  + This is another model of a binary multilayer tree. Here the number of SPFs grows linearly
    with each layer. Another standard model to investigate scaling.

1.) Write a tool that calculates the total number of coefficients for
    a given basis representations. Test the program using the basis sets
    supplied by "nocl.basis" and "ch3.basis".

2.) The tool should now provide more detailed information and should tell
    the number of toplayer, bottomLayer- and intermediate layer coefficients, respectively.
    Adjust your tool accordingly.

3.) Investigate how many coefficients are needed to represent wavefunctions corresponding
    to the different basis sets. Study the relative number of toplayer, intermediate-layer and
    bottomLayer coefficients for the different models.

4.) Calculate the required memory space to save the wavefunction, if each coefficients is
    saved as a complex<double> in binary format (see sizeof() function in cppreference.com).

5.*) Investigate scaling of the number of coefficients with increases system size for the
    binary tree model systems. Create basis files for systems with 8, 16, 32, 64 degrees of
    freedom, evaluate the number of coefficients and make a diagram. Compare it to the
    analytic equations to validate your results. In Ref. [1] a detailed study can be found.

Inputs
======
nocl.in

.. code-block:: YAML


ch3.in
.. code-block:: YAML

References
==========
- [1] U. Manthe, J. Chem. Phys., 128, 164116 (2008)
- [2] U. Manthe, H. D. Meyer, J. Chem. GetLeaf., 97, 3199 (1992)
- [3] U. Manthe, J. Chem. GetLeaf, 130, 054109 (2009)

