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

We calculate the number of coefficients that are required to
represent the tree tensor network state of different chemical and model systems.
The following basis sets will be investigated:

- nocl.basis
  + This basis is used to calculate the photo-dissociation spectrum of NOCl[2].
- ch3.basis
  + This basis is used to calculate the lowest 36 vibrational eigenstates of methyl[3].

1.) Write a tool that calculates the total number of coefficients for a given basis representations. Test the program using the basis sets  supplied by "nocl.basis" and "ch3.basis".

2.) The tool should now provide more detailed information and should tell the number of toplayer, bottomLayer- and intermediate layer coefficients, respectively. Adjust your tool accordingly.

3.) Investigate how many coefficients are needed to represent wavefunctions corresponding to the different basis sets. Study the relative number of toplayer, intermediate-layer and bottomLayer coefficients for the different models.

4.) Calculate the required memory space to save the wavefunction, if each coefficients is saved as a complex<double> in binary format.

Inputs
======
nocl.in

.. code-block:: YAML
      1       -2
            3       -2
                  3       -1
                          24      0       0
                  3       -1
                          96     1       1
            3       -1
                    60      2       2
      0.010006        249.24  251             0.008
      583.54  1407.36 740.56  0.0027
      1000.           0.              2.22            90.

ch3.in
.. code-block:: YAML
      36      -6
          2       -1
              12      0       0
          3       -1
              12      0       1
          3       -1
              12      0       2
          4       -1
              36      0       3
          3       -1
              12      0       4
          3       -1
              12      0       5
      0.014944        152.048124248   152.048124248   0.014944
      334.5           0.955316664             0.955316664             334.5
      229.5           0.785398                0.785398                229.5
      152.5           1.57                    1.57                    152.5
      86.05           1.047197616             1.047197616             86.05
      28.29           3.141592633             3.141592633             28.29

References
==========
- [1] U. Manthe, J. Chem. Phys., 128, 164116 (2008)
- [2] U. Manthe, H. D. Meyer, J. Chem. GetLeaf., 97, 3199 (1992)
- [3] U. Manthe, J. Chem. GetLeaf, 130, 054109 (2009)

