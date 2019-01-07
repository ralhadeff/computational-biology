Monte Carlo and Brownian Dyanamics (Langevin Dynamics) programs for ad-hoc modeling

The 2 dimensions LD simulations is for the Renormalization approach (see https://pubs.acs.org/doi/abs/10.1021/jp1056122)
The complex MC_LD code is a specific implementation (with a hard-coded force-field) for mysin 5 (see https://www.pnas.org/content/114/39/10426)
A 'clean' version of the code is provided, where the FF is deleted, and can be implemented by a user for any purpose.

Please note that the 2D LD code uses the same friction coefficient in x andy axes (this code was discontinued and this issue was never addressed)

example input files are provided
Please contact me (ralhadef@usc.edu) if you have any questions on how to run/modify/interpret results

Markus parabols:

f(x) = a1(x-b1)+c1 
g(x) = a2(x-b2)+c2
e(x) = [f+g-sqrt( (f-g)^2 - 4H12^2 )]/2
where H12 needs to be a small number that makes the conjunction point derivable.
delta-G = delta-c (energy difference between start and end) = c2-c1
delta-r = delta-b (distance in the x coordinate between the start and end points b2-b1)
barrier = c1+(dG+lambda)^2/4lambda-H12
lambda  = (b2-b1)^2 = dr^2
intercept between f and g = [(c2-c1)/(b2-b1)+b1+b2]/2
(NOTE: code uses a simplified function that is faster to derive, and produces the same behavior in MC/LD simulations)
