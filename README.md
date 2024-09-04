To implement OCCAM, first compile "K_medians_K_neq_p_MEX.cpp" with Armadillo (http://arma.sourceforge.net/) in MATLAB. The commandline I used under Linux was:
    mex -L ***SOMEPATH***/armadillo/armadillo-4.100.0/include/armadillo -Iarmadillo ./cpp_subroutine_MEX/K_medians_K_neq_p_MEX.cpp
Linux users: you can attempt to directly use the mex file I compiled "K_medians_K_neq_p_MEX.mexa64".

Then follow the instruction lines in "OCCAM.m" to run the method.

Reference: Detecting overlapping communities in networks using spectral methods, by Yuan Zhang, Elizaveta Levina and Ji Zhu, Arxiv: 1412.3432