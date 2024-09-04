Written by Murdock Aubry:

1.  The file levin2d.f90 contains code for evaluating integrals of the form

         b     d
     \int  \int f(x, y) exp( i g(x, y) ) dx dy,
         a     C
     using the two dimensional adaptive Levin method.

2.  The file levin2d.f90 contains code for evaluating integrals of the form

         b     d
     \int  \int f(x, y) exp( i g(x, y) ) dx day
         a     C
     using the two dimensional adaptive Gauss-Legendre scheme.

Written by James Bremer at the University of Toronto:

3.  The file chebyshev.f90 contains code for constructing and manipulating
univariate Chebyshev expansions.

4.  The file legendre.f90 contains code for constructing and manipulating
univariate Legendre expansions and Gauss-Legendre quadrature rules.

5.  The file bicheb.f90 contains code for constructing and manipulating
bivariate Chebyshev expansions.

6.  The file chebpw.f90 contains code for constructing and manipulating
univariate piecewise Chebyshev expansions on intervals.

7.  The file levin.f90 contains code for evaluating integrals of the form

         b  
     \int   f(x) exp( g(x) ) dx
         a

with f(x) and g(x) smooth via an "adaptive Levin method."  Its principal claim 
to fame is that it is often much more efficient than standard methods when
g(x) is of large magnitude.
