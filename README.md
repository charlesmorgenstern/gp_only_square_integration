# gp_only_square_integration
Integrate the area of a square with gp-only integration

Integrate the area of a square [-.75pi,.75pi]x[-.75pi,.75pi] using gradient path only integration.
Compare the results using a radial and nonradial function to generate the gradient paths.
Check the convergence with number of paths and step size for computing paths and integrating.

nonradial_square_area.jl is the nonradial function library with f(x,y)=2.0+sin(x+pi/2)+sin(y+pi/2).

radial_square_area.jl is the radial function library with f(x,y)=2-x^2-y^2.

The two libraries are identical except the function used.

plotpaths(n,r,h2) plots the gradient paths, function, and domain of integration. Uses n paths starting on a circle of radius r. Path are calculated using Euler's method with step size h2.

getpaths(n,r,h2) generates an array of gradient paths. Uses n paths starting on a circle of radius r. Path are calculated using Euler's method with step size h2.

curv(r) generates the curvature of an isosurface at point r.

int_area(n,h,intradius) calculates the area of the square using gp-only integration. Uses n paths starting on a circle of radius intradius. Path are calculated using Euler's method with step size h. Integration of curvature is done with right hand rule using step size h. Integration of l1(s) is done with trapezoidal rule with step size h. Returns relative error.

errortablen() generates the error and convergence table for varying n, the number of gradient paths.

errortableh() generates the error and convergence table for varying h, the step size for calculating paths and integrating.
