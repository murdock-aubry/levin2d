!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing and manipulating bivariate Chebyshev
!  expansions of the forms
!
!                      sum              a       T (x) T (y) ,                               (1)
!                    i+j <= norder       ij      i     j    
!
!
!                      sum              a       T (x) T (y)                                 (2)
!               i^2+j^2 <= norder^2      ij      i     j    
!
!  or 
!                       sum             a       T (x) T (y).                                (3)
!               0 <= i,j <= norder       ij      i     j    
!
!  The routines here take a parameter called itype which controls which type of
!  expansion is being used.  
!  
!  The following subroutines should be regarded as publicly callable:
!
!    bicheb_nterms - return the number of terms in an expansion 
!         
!    bicheb_grid - return the nodes of the (norder+1)-point Chebyshev extrema grid
!      on the interval [-1,1]
!
!    bicheb_coefs - given the values of an expansion at the nodes of the tensor product 
!      of the extrema grid returned by bicheb_grid, calculate the coefficients in the 
!      expansion 
!
!    bicheb_eval - evaluate an expansion at a specified point given its coefficients
!
!    bicheb_interp - evaluate an expansion of either form given its values at the
!      nodes of the tensor product of the extrema grid using barycentric Lagrange
!      interpolation
!
!    bicheb_acoefs - return the values-to-coefficients matrix which takes the 
!      vector of values
!
!         (    f(x_0,x_1)   ) 
!         (    f(x_1,x_1)   ) 
!         (        ...      ) 
!         ( f(x_norder,x_1) ) 
!
!     to the vector of coefficients

!
!    bicheb_adiff - return two spectral differentiation matrices, one which calculates
!      derivative w.r.t. x and one that computes the derivative w.r.t. y
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bicheb

interface         bicheb_coefs
module procedure  bicheb_coefs1
module procedure  bicheb_coefs2
end interface     bicheb_coefs

interface         bicheb_eval
module procedure  bicheb_eval1
module procedure  bicheb_eval2
end interface     bicheb_eval

interface         bicheb_interp
module procedure  bicheb_interp1
module procedure  bicheb_interp2
end interface     bicheb_interp

contains

subroutine bicheb_nterms(itype,norder,nterms)
implicit double precision (a-h,o-z)
!
!  Return the number of terms in an expansion of the form (1) or (2).
!
!  Input parameters:
!    itype - an integer parameter which controls which form the expansion takes;
!        itype = 1  indicates an expansion of the first kind, i.e., one of the form (1)
!        itype = 2  indicates an expansion of the second kind, i.e., one of the form (2)   
!        itype = 3  indicates an expansion of the third kind, i.e., one of the form (3)   
!    norder - the order of the expansion
!
!  Output parameters:
!    nterms - the number of terms in the expansion
!

if (itype .eq. 1) then

nterms = 0
do i=0,norder
do j=0,norder-i
nterms = nterms+1
end do
end do


elseif (itype .eq. 2) then

nterms = 0
do i=0,norder
do j=0,floor(sqrt(norder**2-i**2+0.0d0))
nterms = nterms+1
end do
end do

else

nterms = 0
do i=0,norder
do j=0,norder
nterms = nterms+1
end do
end do

endif

end subroutine


subroutine bicheb_grid(norder,xs)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out) :: xs(:)
!
!  Return the (norder+1)-point Chebyshev extrema grid on the interval
!  (-1,1).
!
!  Input parameters:
!    norder - indicates the number of points in the grid
!
!  Output parameters:
!    xs - an array of length norder+1 specifying the points in the grid
!
!
data pi / 3.14159265358979323846264338327950288d0 /

allocate(xs(0:norder))
do i=0,norder
xs(i) = cos(pi*(norder-i)/(norder+0.0d0))
end do
end subroutine


subroutine bicheb_coefs1(itype,norder,vals,coefs)
implicit double precision (a-h,o-z)
double precision :: vals(0:norder,0:norder)
double precision :: coefs(:)
!
!  Calculate the coefficients in an expansion given its values at the nodes of the tensor 
!  product of the Chebyshev extrema grid returned by bicheb_grid.
!
!  Input parameters:
!    itype - an integer parameter which controls which form the expansion takes;
!        itype = 1  indicates an expansion of the first kind, i.e., one of the form (1)
!        itype = 2  indicates an expansion of the second kind, i.e., one of the form (2)   
!        itype = 3  indicates an expansion of the third kind, i.e., one of the form (3)   
!    norder - the order of the expansion (1)
!    vals - an array with shape (0:norder,0:norder) whose
!      i,j entry gives the value f(x_i, x_j) of the expansion (1),
!      where x_0,...,x_norder are the nodes of the discretization grid
!      returned by bicheb_grid
!
!  Output parameters:
!    coefs - an array of length nterms containing the coefficients in the
!      expansion
!
data pi / 3.14159265358979323846264338327950288d0 /


idx = 0
dd  = (2.0d0/(norder+0.0d0))**2

if (itype .eq. 1) then

do i=0,norder
do j=0,norder-i

  dsum = 0
  do l1=0,norder

  dd1 = 1
  if (l1 .eq. 0 .OR. l1 .eq. norder) dd1 = 0.5d0
  dd1 = dd1 * cos(pi*(norder-l1)*i/(norder+0.0d0))
  do l2=0,norder

  dd2 = 1
  if (l2 .eq. 0 .OR. l2 .eq. norder) dd2 = 0.5d0
  dd2 = dd2 * cos(pi*(norder-l2)*j/(norder+0.0d0))


  dsum = dsum + vals(l1,l2)*dd1*dd2*dd

  end do
  end do

  if ( i == 0 .OR. i == norder) dsum = dsum/2 
  if ( j == 0 .OR. j == norder) dsum = dsum/2 

  idx        = idx + 1
  coefs(idx) = dsum

end do
end do

elseif (itype .eq. 2) then

do i=0,norder
do j=0,floor(sqrt(norder**2-i**2+0.0d0))

  dsum = 0
  do l1=0,norder

  dd1 = 1
  if (l1 .eq. 0 .OR. l1 .eq. norder) dd1 = 0.5d0
  dd1 = dd1 * cos(pi*(norder-l1)*i/(norder+0.0d0))
  do l2=0,norder

  dd2 = 1
  if (l2 .eq. 0 .OR. l2 .eq. norder) dd2 = 0.5d0
  dd2 = dd2 * cos(pi*(norder-l2)*j/(norder+0.0d0))


  dsum = dsum + vals(l1,l2)*dd1*dd2*dd

  end do
  end do

  if ( i == 0 .OR. i == norder) dsum = dsum /2 
  if ( j == 0 .OR. j == norder) dsum = dsum /2 

  idx        = idx + 1
  coefs(idx) = dsum

end do
end do

else

do i=0,norder
do j=0,norder

  dsum = 0
  do l1=0,norder

  dd1 = 1
  if (l1 .eq. 0 .OR. l1 .eq. norder) dd1 = 0.5d0
  dd1 = dd1 * cos(pi*(norder-l1)*i/(norder+0.0d0))
  do l2=0,norder

  dd2 = 1
  if (l2 .eq. 0 .OR. l2 .eq. norder) dd2 = 0.5d0
  dd2 = dd2 * cos(pi*(norder-l2)*j/(norder+0.0d0))


  dsum = dsum + vals(l1,l2)*dd1*dd2*dd

  end do
  end do

  if ( i == 0 .OR. i == norder) dsum = dsum /2 
  if ( j == 0 .OR. j == norder) dsum = dsum /2 

  idx        = idx + 1
  coefs(idx) = dsum

end do
end do

endif

end subroutine


subroutine bicheb_coefs2(itype,norder,vals,coefs)
implicit double precision (a-h,o-z)
double complex :: vals(0:norder,0:norder), coefs(:)
!
!  Calculate the coefficients in an expansion given its values at the nodes of the tensor 
!  product of the Chebyshev extrema grid returned by bicheb_grid.
!
!  Input parameters:
!    itype - an integer parameter which controls which form the expansion takes;
!        itype = 1  indicates an expansion of the first kind, i.e., one of the form (1)
!        itype = 2  indicates an expansion of the second kind, i.e., one of the form (2)   
!        itype = 3  indicates an expansion of the third kind, i.e., one of the form (3)   
!    norder - the order of the expansion (1)
!    vals - an array with shape (0:norder,0:norder) whose
!      i,j entry gives the value f(x_i, x_j) of the expansion (1),
!      where x_0,...,x_norder are the nodes of the discretization grid
!      returned by bicheb_grid
!
!  Output parameters:
!    coefs - an array of length nterms containing the coefficients in the
!      expansion
!
data pi / 3.14159265358979323846264338327950288d0 /

double complex :: dsum

idx = 0
dd  = (2.0d0/(norder+0.0d0))**2

if (itype .eq. 1) then

do i=0,norder
do j=0,norder-i

  dsum = 0
  do l1=0,norder

  dd1 = 1
  if (l1 .eq. 0 .OR. l1 .eq. norder) dd1 = 0.5d0
  dd1 = dd1 * cos(pi*(norder-l1)*i/(norder+0.0d0))
  do l2=0,norder

  dd2 = 1
  if (l2 .eq. 0 .OR. l2 .eq. norder) dd2 = 0.5d0
  dd2 = dd2 * cos(pi*(norder-l2)*j/(norder+0.0d0))


  dsum = dsum + vals(l1,l2)*dd1*dd2*dd

  end do
  end do

  if ( i == 0 .OR. i == norder) dsum = dsum /2 
  if ( j == 0 .OR. j == norder) dsum = dsum /2 

  idx        = idx + 1
  coefs(idx) = dsum

end do
end do

elseif (itype .eq. 2) then

do i=0,norder
do j=0,floor(sqrt(norder**2-i**2+0.0d0))

  dsum = 0
  do l1=0,norder

  dd1 = 1
  if (l1 .eq. 0 .OR. l1 .eq. norder) dd1 = 0.5d0
  dd1 = dd1 * cos(pi*(norder-l1)*i/(norder+0.0d0))
  do l2=0,norder

  dd2 = 1
  if (l2 .eq. 0 .OR. l2 .eq. norder) dd2 = 0.5d0
  dd2 = dd2 * cos(pi*(norder-l2)*j/(norder+0.0d0))


  dsum = dsum + vals(l1,l2)*dd1*dd2*dd

  end do
  end do

  if ( i == 0 .OR. i == norder) dsum = dsum /2 
  if ( j == 0 .OR. j == norder) dsum = dsum /2 

  idx        = idx + 1
  coefs(idx) = dsum

end do
end do

else

do i=0,norder
do j=0,norder

  dsum = 0
  do l1=0,norder

  dd1 = 1
  if (l1 .eq. 0 .OR. l1 .eq. norder) dd1 = 0.5d0
  dd1 = dd1 * cos(pi*(norder-l1)*i/(norder+0.0d0))
  do l2=0,norder

  dd2 = 1
  if (l2 .eq. 0 .OR. l2 .eq. norder) dd2 = 0.5d0
  dd2 = dd2 * cos(pi*(norder-l2)*j/(norder+0.0d0))


  dsum = dsum + vals(l1,l2)*dd1*dd2*dd

  end do
  end do

  if ( i == 0 .OR. i == norder) dsum = dsum /2 
  if ( j == 0 .OR. j == norder) dsum = dsum /2 

  idx        = idx + 1
  coefs(idx) = dsum

end do
end do

endif

end subroutine


subroutine bicheb_eval1(itype,norder,coefs,x,y,val)
implicit double precision (a-h,o-z)
double precision :: coefs(:)
!
!  Evaluate an expansion at specified point in the rectangle [-1,1] x [-1,1].
!
!  Input parameters:
!    itype - an integer parameter which controls which form the expansion takes;
!        itype = 1  indicates an expansion of the first kind, i.e., one of the form (1)
!        itype = 2  indicates an expansion of the second kind, i.e., one of the form (2)   
!        itype = 3  indicates an expansion of the third kind, i.e., one of the form (3)   
!    norder - the order of the expansion 
!    coefs - an array containing the coefficients of the expansion, as
!      ordered as in the bicheb_coefs routine
!    (x,y) - the point at which to evaluate the expansion
!
!  Output parameters:
!    val - the desired value of the expansion
!
!
double precision :: polsx(0:norder), polsy(0:norder)
data pi / 3.14159265358979323846264338327950288d0 /

polsx(0) = 1
polsx(1) = x
polsy(0) = 1
polsy(1) = y
do i=1,norder-1
polsx(i+1) = 2*x*polsx(i)-polsx(i-1)
polsy(i+1) = 2*y*polsy(i)-polsy(i-1)
end do

if (itype .eq. 1) then

idx = 0
!dd  = (2.0d0/(norder+0.0d0))**2
val = 0

do i=0,norder
do j=0,norder-i
idx = idx + 1 
val = val + polsx(i)*polsy(j)*coefs(idx)
end do
end do

elseif (itype .eq. 2) then

idx = 0
dd  = (2.0d0/(norder+0.0d0))**2
val = 0

do i=0,norder
do j=0,floor(sqrt(norder**2-i**2+0.0d0))
idx = idx + 1 
val = val + polsx(i)*polsy(j)*coefs(idx)
end do
end do

else

idx = 0
dd  = (2.0d0/(norder+0.0d0))**2
val = 0

do i=0,norder
do j=0,norder
idx = idx + 1 
val = val + polsx(i)*polsy(j)*coefs(idx)
end do
end do

endif

end subroutine


subroutine bicheb_eval2(itype,norder,coefs,x,y,val)
implicit double precision (a-h,o-z)
double complex :: coefs(:)
!
!  Evaluate an expansion at specified point in the rectangle [-1,1] x [-1,1].
!
!  Input parameters:
!    itype - an integer parameter which controls which form the expansion takes;
!        itype = 1  indicates an expansion of the first kind, i.e., one of the form (1)
!        itype = 2  indicates an expansion of the second kind, i.e., one of the form (2)   
!        itype = 3  indicates an expansion of the third kind, i.e., one of the form (3)   
!    norder - the order of the expansion 
!    coefs - an array containing the coefficients of the expansion, as
!      ordered as in the bicheb_coefs routine
!    (x,y) - the point at which to evaluate the expansion
!
!  Output parameters:
!    val - the desired value of the expansion
!
!
double precision :: polsx(0:norder), polsy(0:norder)
double complex   :: val

data pi / 3.14159265358979323846264338327950288d0 /

polsx(0) = 1
polsx(1) = x
polsy(0) = 1
polsy(1) = y
do i=1,norder-1
polsx(i+1) = 2*x*polsx(i)-polsx(i-1)
polsy(i+1) = 2*y*polsy(i)-polsy(i-1)
end do

if (itype .eq. 1) then

idx = 0
!dd  = (2.0d0/(norder+0.0d0))**2
val = 0

do i=0,norder
do j=0,norder-i
idx = idx + 1 
val = val + polsx(i)*polsy(j)*coefs(idx)
end do
end do

elseif (itype .eq. 2) then

idx = 0
dd  = (2.0d0/(norder+0.0d0))**2
val = 0

do i=0,norder
do j=0,floor(sqrt(norder**2-i**2+0.0d0))
idx = idx + 1 
val = val + polsx(i)*polsy(j)*coefs(idx)
end do
end do

else

idx = 0
dd  = (2.0d0/(norder+0.0d0))**2
val = 0

do i=0,norder
do j=0,norder
idx = idx + 1 
val = val + polsx(i)*polsy(j)*coefs(idx)
end do
end do

endif

end subroutine



subroutine bicheb_interp1(norder,vals,x,y,val)
implicit double precision (a-h,o-z)
double precision :: vals(0:norder,0:norder)
!
!  Use barcyentric Lagrange interpolation to evaluate an expansion of either form
!  given its values at the tensor product of the Chebyshev extrema grid.
!
!  Input parameters:
!    norder - the order of the expansion
!    vals - an array whose i,j entry gives the value of the expansion at the
!      point (x_i, x_j), where x_0, ..., x_norder are the nodes of the
!      Chebyshev extrema grid
!
!    (x,y) - the point at which to evaluate 
!
!
!  Output parameters:
!    val - the value of the expansion
!

double precision xs(0:norder), whtsx(0:norder), whtsy(0:norder)

data pi / 3.14159265358979323846264338327950288d0 /

eps0  = epsilon(0.0d0)
dsum1 = 0
dsum2 = 0
dsign = 1.0d0
idx_x = -1
idx_y = -1

do i=0,norder

xs(i)  = cos(pi*(norder-i)/(norder+0.0d0))
diffx  = x-xs(i)
diffy  = y-xs(i)

if (abs(diffx) .lt. eps0) idx_x = i
if (abs(diffy) .lt. eps0) idx_y = i

dd    = dsign
if (i .eq. 0 .OR. i .eq. norder) dd = dd/2

whtsx(i) = dd/diffx
whtsy(i) = dd/diffy
dsum1    = dsum1 + whtsx(i)
dsum2    = dsum2 + whtsy(i)
dsign    = -dsign
end do

if (idx_x .ge. 0) then

if(idx_y .ge. 0) then
val = vals(idx_x,idx_y)
return
else
goto 1000
endif

endif

if(idx_y .ge. 0) goto 2000

dsum3 = 0
do i=0,norder
do j=0,norder
dsum3 = dsum3 + vals(i,j)*whtsx(i)*whtsy(j)
end do
end do

val = dsum3/(dsum1*dsum2)

return

! x is equal to one of the interpolation nodes but y is not
1000 continue

dsum3 = 0
do i=0,norder
dsum3 = dsum3 + vals(idx_x,i)*whtsy(i)
end do
val = dsum3/dsum2

return

! y is equal to one of the interpolation nodes but x is not
2000 continue

dsum3 = 0
do i=0,norder
dsum3 = dsum3 + vals(i,idx_y)*whtsx(i)
end do
val = dsum3/dsum1
return

end subroutine


subroutine bicheb_interp2(norder,vals,x,y,val)
implicit double precision (a-h,o-z)
double complex :: vals(0:norder,0:norder)
!
!  Use barcyentric Lagrange interpolation to evaluate an expansion of either form
!  given its values at the tensor product of the Chebyshev extrema grid.
!
!  Input parameters:
!    norder - the order of the expansion
!    vals - an array whose i,j entry gives the value of the expansion at the
!      point (x_i, x_j), where x_0, ..., x_norder are the nodes of the
!      Chebyshev extrema grid
!
!    (x,y) - the point at which to evaluate 
!
!
!  Output parameters:
!    val - the value of the expansion
!

double precision :: xs(0:norder), whtsx(0:norder), whtsy(0:norder)
double complex   :: dsum3

data pi / 3.14159265358979323846264338327950288d0 /

eps0  = epsilon(0.0d0)
dsum1 = 0
dsum2 = 0
dsign = 1.0d0
idx_x = -1
idx_y = -1

do i=0,norder

xs(i)  = cos(pi*(norder-i)/(norder+0.0d0))
diffx  = x-xs(i)
diffy  = y-xs(i)

if (abs(diffx) .lt. eps0) idx_x = i
if (abs(diffy) .lt. eps0) idx_y = i

dd    = dsign
if (i .eq. 0 .OR. i .eq. norder) dd = dd/2

whtsx(i) = dd/diffx
whtsy(i) = dd/diffy
dsum1    = dsum1 + whtsx(i)
dsum2    = dsum2 + whtsy(i)
dsign    = -dsign
end do

if (idx_x .ge. 0) then

if(idx_y .ge. 0) then
val = vals(idx_x,idx_y)
return
else
goto 1000
endif

endif

if(idx_y .ge. 0) goto 2000

dsum3 = 0
do i=0,norder
do j=0,norder
dsum3 = dsum3 + vals(i,j)*whtsx(i)*whtsy(j)
end do
end do

val = dsum3/(dsum1*dsum2)

return

! x is equal to one of the interpolation nodes but y is not
1000 continue

dsum3 = 0
do i=0,norder
dsum3 = dsum3 + vals(idx_x,i)*whtsy(i)
end do
val = dsum3/dsum2

return

! y is equal to one of the interpolation nodes but x is not
2000 continue

dsum3 = 0
do i=0,norder
dsum3 = dsum3 + vals(i,idx_y)*whtsx(i)
end do
val = dsum3/dsum1
return

end subroutine


subroutine bicheb_acoefs(itype,norder,acoefs)
implicit double precision (a-h,o-z)
double precision, allocatable     :: acoefs(:,:)
!
!  Construct the ( nterms, (norder+1)^2 ) matrix which takes the values of an expansion
!  at the Chebyshev tensor product nodes to the coefficients in the expansion.
! 
!  Input parameters:
!    itype - an integer parameter which controls which form the expansion takes;
!        itype = 1  indicates an expansion of the first kind, i.e., one of the form (1)
!        itype = 2  indicates an expansion of the second kind, i.e., one of the form (2)   
!        itype = 3  indicates an expansion of the third kind, i.e., one of the form (3)   
!    norder - the order of the expansion (1)
!
!  Output parameters:
!    acoefs - the ( nterms, (norder+1)^2 ) matrix taking the 
!
!



end subroutine

end module
