# DG_lyc
A standard P2 2-D Discontinuous-Galerkin-Method program written by Fortran. This program includes:

1. Rectangular grid

2. P^2 Basis：{1, phi_1(X), phi_1(Y), phi_2(X), phi_1(X)*phi_1(Y), phi_2(Y)}.
Where phi_0 = 1, phi_1 = X, phi_2 = X^2 - 1/3， X = (x - x_i)/(Delta x/2)

3. Initial value determined by L^2 projection

4. Three Fluxes: L-F, HLL and HLLC

5. One step Euler forward and RK3 time iteration modes

6. The TVB-Limiter (component wise)

7. Examples: Smooth Vortex, Orszag Tang Vortex, Rotor

8. L^2 error calculation

9. Save and draw the numerical solution(MATLAB)

You can modify the text 'com.txt' and the subroutine 'init_data.f90' to change above parameters.
