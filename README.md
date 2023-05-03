# Vibration-of-a-loaded-side-cracked-rectangular-plate
This MATLAB code verifies the work of Zeng, et al, in 2016: https://www.sciencedirect.com/science/article/pii/S026382311630088X

The Rayleigh-Ritz Method.
Kirchhoff plate theory.

The scripts of numerical integration were referred to the following links:
Quadrature weights and points for numerically integrating over a triangle: https://au.mathworks.com/matlabcentral/fileexchange/72131-quadtriangle
Quadrature weights and points for numerically integrating over a square: https://au.mathworks.com/matlabcentral/fileexchange/72151-quadsquare-quadrature-rules-for-the-square?s_tid=prof_contriblnk
Uniform triangulation refinement: https://au.mathworks.com/matlabcentral/fileexchange/97737-trirefine?s_tid=prof_contriblnk

The stress distribution is calculated using XFEM.
The two main scripts employ code-writing techniques such as parallel computation and computation acceleration. 
