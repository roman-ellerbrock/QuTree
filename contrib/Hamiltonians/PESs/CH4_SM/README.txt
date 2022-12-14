Supplementary material for:

Title: A highly accurate ab initio potential energy surface for methane

Authors: A. Owens, S. N. Yurchenko, A. Yachmenev, J. Tennyson, and W. Thiel


1. xy4_PES.f90 : Source fortran code to construct potential energy surface for methane.
                 Program will print potential energy at chosen geometries for respective potential energy surface.


2. ch4_cbs-f12-hl_PES.inp : Input file for the CBS-F12^{HL} PES to use with fortran program xy4_PES.f90.

 Note: The order of the geometry specification in the input file after the parameter set is:

 r(C-H1)/Ang  r(C-H2)/Ang  r(C-H3)/Ang  r(C-H4)/Ang  alpha12(H1-C-H2)/deg  alpha13(H1-C-H3)/deg

 alpha14(H1-C-H4)/deg  alpha23(H2-C-H3)/deg  alpha24(H2-C-H4)/deg  alpha34(H3-C-H4)/deg


3. ch4_cbs-f12-hl_PES.out : Example output file 


4. ch4_cbs-f12-hl_TROVE_Pmax-14_energies.txt : List of computed energy levels from TROVE.

 Order of the file is:  

 Symmetry | Running number | Energy (cm-1) | TROVE quantum numbers (x9) | Basis function contribution |

 e.g.  A1       1      0.000000    0   0   0   0   0   0   0   0   0   0.99

 Note that calculations were carried out with P_max=14 and therefore may not be converged at higher energies.
