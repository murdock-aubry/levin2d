!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for forming and manipulating bivariate Chebyshev expansions
!  of the form
!
!             n-1   n-1
!    f(x,y) = sum   sum  a_ij T  (x)  T  (y)                                             (1)
!             i=0   j=0        i       j
!
!  This code represents such expansions using two different mechanisms: via the
!  vector
!
!     a_0,0
!     a_1,0
!     a_2,0
!     ...                                                                                (2)
!     a_n-1,0
!     a_0,1
!     ...
!     a_n-1,1
!     ...
!     a_0,n-1
!     ...
!     a_n-1,n-1
!
!  of expansion coefficients and via the vector 
!
!    f(x_1,y_1)
!    f(x_2,y_1)
!       ...
!    f(x_n,y_1)
!    f(x_1,y_2)
!       ...
!    f(x_n,y_2)                                                                          (3)
!       ...
!    f(x_1,y_n)
!       ...                                                                              
!    f(x_n,y_n)
!
!  of values of the expansion at the nodes of the Chebyshev tensor product
!  discretization grid
!
!      (x_i,y_i) :  i = 1, ..., n                                                        (4)
!
!  The following subroutines should be regarded as publicly callable:
!
!    bicheb2_nodes - return the nodes of the tensor product grid
!
!    bicheb2_interp - evaluate an expansion of the form (1) given the vector of 
!     values (2)
!
!    bicheb2_eval - evaluate an expansion of the form (1) given the vector of
!      coefficients (2)
!
!    bicheb2_xslice - 
!
!    bicheb2_yslice - 
!
!    bicheb2_coefsmatrix - return the values to coefficients matrix
!
!    bicheb2_diffxmatrix - return the spectral differentiation matrix which
!       calculates the values of the partial derivative of (1) w.r.t. x on
!       the tensor product grid
!
!    bicheb2_diffxxmatrix - return the spectral differentiation matrix which
!       calculates the values of derivative of d^2 f / dx^2 on the tensor
!       product grid
!
!    bicheb2_diffymatrix - return the spectral differentiation matrix which
!       calculates the values of the partial derivative of (1) w.r.t. y on
!       the tensor product grid
!
!    bicheb2_diffyymatrix - return the spectral differentiation matrix which
!       calculates the values of derivative of d^2 f / dy^2 on the tensor
!       product grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bicheb2
use utils
use chebyshev


interface        bicheb2_eval
module procedure bicheb2_eval1
module procedure bicheb2_eval2
end interface    bicheb2_eval

interface        bicheb2_interp
module procedure bicheb2_interp1
module procedure bicheb2_interp2
end interface    bicheb2_interp

interface        bicheb2_xslice
module procedure bicheb2_xslice1
module procedure bicheb2_xslice2
end interface    bicheb2_xslice


interface        bicheb2_yslice
module procedure bicheb2_yslice1
module procedure bicheb2_yslice2
end interface    bicheb2_yslice
 
contains


subroutine bicheb2_xslice1(n,vals,y,xvals)
implicit double precision (a-h,o-z)
double precision                           :: vals(:)
double precision, allocatable, intent(out) :: xvals(:)
!
!
!
!  Input parameters:
!
!  Output parameters:
!
!
!

data pi / 3.14159265358979323846264338327950288d0 /

allocate(xvals(n))
h      = pi/(n-1)
do i=1,n
x     = cos(h*(n-i)) 
call bicheb2_interp(n,vals,x,y,xvals(i))
end do

end subroutine



subroutine bicheb2_xslice2(n,vals,y,xvals)
implicit double precision (a-h,o-z)
double complex                           :: vals(:)
double complex, allocatable, intent(out) :: xvals(:)
!
!
!
!  Input parameters:
!
!  Output parameters:
!
!
!

data pi / 3.14159265358979323846264338327950288d0 /

allocate(xvals(n))
h      = pi/(n-1)
do i=1,n
x     = cos(h*(n-i)) 
call bicheb2_interp(n,vals,x,y,xvals(i))
end do

end subroutine



subroutine bicheb2_yslice1(n,vals,x,yvals)
implicit double precision (a-h,o-z)
double precision                           :: vals(:)
double precision, allocatable, intent(out) :: yvals(:)
!
!
!
!  Input parameters:
!
!  Output parameters:
!
!
!

data pi / 3.14159265358979323846264338327950288d0 /

allocate(yvals(n))
h      = pi/(n-1)
do i=1,n
y     = cos(h*(n-i)) 
call bicheb2_interp(n,vals,x,y,yvals(i))
end do

end subroutine


subroutine bicheb2_yslice2(n,vals,x,yvals)
implicit double precision (a-h,o-z)
double complex                           :: vals(:)
double complex, allocatable, intent(out) :: yvals(:)
!
!
!
!  Input parameters:
!
!  Output parameters:
!
!
!

data pi / 3.14159265358979323846264338327950288d0 /

allocate(yvals(n))
h      = pi/(n-1)
do i=1,n
y     = cos(h*(n-i)) 
call bicheb2_interp(n,vals,x,y,yvals(i))
end do

end subroutine


subroutine bicheb2_eval1(n,coefs,x,y,val)
implicit double precision (a-h,o-z)
double precision      :: coefs(:), val
!
!  Evaluate a real-valued expansion of the form (1) at a specified point (x,y) given
!  the vector (2) of coefficients.
!
!  Input parameters:
!    n - parameter controlling the number of terms in the expansion (1)
!    coefs - the vector of expansion coefficients (2)
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the caclulate value of (1)
!
double precision      :: polsx(n), polsy(n)

call chebs(n,x,polsx)
call chebs(n,y,polsy)

val = 0
idx = 0

do j=1,n
do i=1,n
idx = idx+1
val = val + coefs(idx)*polsx(i)*polsy(j)
end do
end do

end subroutine


subroutine bicheb2_eval2(n,coefs,x,y,val)
implicit double precision (a-h,o-z)
double complex      :: coefs(:), val
!
!  Evaluate a complex-valued expansion of the form (1)  at a specified point (x,y) given the
!  vector (2) of coefficients.
!
!  Input parameters:
!    n - parameter controlling the number of terms in the expansion (1)
!    coefs - the vector of expansion coefficients (2)
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the caclulate value of (1)
!

double precision      :: polsx(n), polsy(n)
call chebs(n,x,polsx)
call chebs(n,y,polsy)

val = 0
idx = 0

do j=1,n
do i=1,n
idx = idx+1
val = val + coefs(idx)*polsx(i)*polsy(j)
end do
end do

end subroutine


subroutine bicheb2_interp1(n,vals,x,y,val)
implicit double precision (a-h,o-z)
double precision          :: vals(:), val
!
!  Evaluate a real-valued expansion of the form (1)  at a specified point (x,y) given the
!  vector (3) of values
!
!  Input parameters:
!    n - parameter controlling the number of terms in the expansion (1)
!    coefs - the vector (3) of values of the expansion at the tensor product grid
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the caclulate value of (1)
!

double precision :: whtsx(1:n), whtsy(1:n)

data pi / 3.14159265358979323846264338327950288d0 /

eps0   = epsilon(0.0d0)
h      = pi/(n-1)
idx_x  = -1
idx_y  = -1
dsign  = 1.0d0

do i=1,n

xx     = cos(h*(n-i)) 
diffx  = x-xx
diffy  = y-xx

if (abs(diffx) .eq. 0) idx_x = i
if (abs(diffy) .eq. 0) idx_y = i

dd    = dsign
if (i .eq. 1 .OR. i .eq. n) dd = dd/2

whtsx(i) = dd/diffx
whtsy(i) = dd/diffy
dsign    = -dsign
end do


if (idx_x .gt. 0 .AND. idx_y .gt. 0) then
val = vals(idx_x + (idx_y-1)*n)
return
endif

if (idx_x .gt. 0) then
dsum3 = 0
do i=1,n
dsum3 = dsum3 + vals(idx_x + (i-1)*n)*whtsy(i)
end do
val = dsum3/sum(whtsy)
return
endif

if (idx_y .gt. 0) then
dsum3 = 0
do i=1,n
dsum3 = dsum3 + vals(i + (idx_y-1)*n)*whtsx(i)
end do
val = dsum3/sum(whtsx)
return
endif


dsum3 = 0
do j=1,n
i1 = 1+(j-1)*n
i2 = n+(j-1)*n
dsum3 = dsum3 + sum(vals(i1:i2)*whtsx)*whtsy(j)
end do
val = dsum3/(sum(whtsx)*sum(whtsy))
return

end subroutine

subroutine bicheb2_interp2(n,vals,x,y,val)
implicit double precision (a-h,o-z)
double complex            :: vals(:), val
!
!  Evaluate a complex-valued expansion of the form (1)  at a specified point (x,y) given the
!  vector (3) of values
!
!  Input parameters:
!    n - parameter controlling the number of terms in the expansion (1)
!    coefs - the vector (3) of values of the expansion at the tensor product grid
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the caclulate value of (1)
!

double precision :: whtsx(1:n), whtsy(1:n)

data pi / 3.14159265358979323846264338327950288d0 /

eps0   = epsilon(0.0d0)
h      = pi/(n-1)
idx_x  = -1
idx_y  = -1
dsign  = 1.0d0

do i=1,n

xx     = cos(h*(n-i)) 
diffx  = x-xx
diffy  = y-xx

if (abs(diffx) .eq. 0) idx_x = i
if (abs(diffy) .eq. 0) idx_y = i

dd    = dsign
if (i .eq. 1 .OR. i .eq. n) dd = dd/2

whtsx(i) = dd/diffx
whtsy(i) = dd/diffy
dsign    = -dsign
end do

if (idx_x .gt. 0 .AND. idx_y .gt. 0) then
val = vals(idx_x + (idx_y-1)*n)
return
endif

if (idx_x .gt. 0) then
dsum3 = 0
do i=1,n
dsum3 = dsum3 + vals(idx_x + (i-1)*n)*whtsy(i)
end do
val = dsum3/sum(whtsy)
return
endif

if (idx_y .gt. 0) then
dsum3 = 0
do i=1,n
dsum3 = dsum3 + vals(i + (idx_y-1)*n)*whtsx(i)
end do
val = dsum3/sum(whtsx)
return
endif

dsum3 = 0
do j=1,n
i1 = 1+(j-1)*n
i2 = n+(j-1)*n
dsum3 = dsum3 + sum(vals(i1:i2)*whtsx)*whtsy(j)
end do
val = dsum3/(sum(whtsx)*sum(whtsy))
return

end subroutine


subroutine bicheb2_coefsmatrix(n,acoefs)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: acoefs(:,:)
!
!  Return the matrix which takes the vector of values (3) to the vector of coefficients (2).
!
!  Input parameters:
!    n - parameter controlling the number of terms in the expansion (1)
!
!  Output parameters:
!    acoefs - the values to coefficients matrix
!
double precision, allocatable  :: acoefs0(:,:), eye0(:,:)
double precision, allocatable  :: a(:,:), b(:,:)

allocate(acoefs(n**2,n**2),eye0(n,n))

eye0 = 0
do i=1,n
eye0(i,i) = 1
end do

call chebyshev_coefsmatrix(n,acoefs0)
call bicheb2_kron(n,acoefs0,eye0,a)
call bicheb2_kron(n,eye0,acoefs0,b)

acoefs = matmul(a,b)

end subroutine


subroutine bicheb2_diffxmatrix(n,adiffx)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)      :: adiffx(:,:)
!
!  Return the matrix which takes the vector of values (3) to the vector of 
!  values of the derivative w.r.t. x of (1).
!
!  Input parameters:
!    n - parameter controlling the number of terms in the expansion (1)
!
!  Output parameters:
!    adiffx - the spectral differentiation matrix

double precision, allocatable :: adiff0(:,:), eye0(:,:)

allocate(eye0(n,n))

eye0 = 0
do i=1,n
eye0(i,i) = 1
end do

call chebyshev_diffmatrix(n,adiff0)
call bicheb2_kron(n,eye0,adiff0,adiffx)

end subroutine


subroutine bicheb2_diffxxmatrix(n,adiffxx)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)      :: adiffxx(:,:)
!
!  Return the matrix which takes the vector of values (3) to the vector of 
!  values of the derivative w.r.t. x of (1).
!
!  Input parameters:
!    n - parameter controlling the number of terms in the expansion (1)
!
!  Output parameters:
!    adiffx - the spectral differentiation matrix

double precision, allocatable :: adiff0(:,:), eye0(:,:), adiff1(:,:)

allocate(eye0(n,n),adiff1(n,n))

eye0 = 0
do i=1,n
eye0(i,i) = 1
end do

call chebyshev_diffmatrix(n,adiff0)
adiff1 = matmul(adiff0,adiff0)
call bicheb2_kron(n,eye0,adiff1,adiffxx)

end subroutine

subroutine bicheb2_diffyymatrix(n,adiffyy)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)      :: adiffyy(:,:)
!
!  Return the matrix which takes the vector of values (3) to the vector of 
!  values of the derivative w.r.t. x of (1).
!
!  Input parameters:
!    n - parameter controlling the number of terms in the expansion (1)
!
!  Output parameters:
!    adiffx - the spectral differentiation matrix

double precision, allocatable :: adiff0(:,:), eye0(:,:), adiff1(:,:)

allocate(eye0(n,n),adiff1(n,n))

eye0 = 0
do i=1,n
eye0(i,i) = 1
end do

call chebyshev_diffmatrix(n,adiff0)
adiff1 = matmul(adiff0,adiff0)
call bicheb2_kron(n,adiff1,eye0,adiffyy)

end subroutine


subroutine bicheb2_diffymatrix(n,adiffy)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)      :: adiffy(:,:)
!
!  Return the matrix which takes the vector of values (3) to the vector of 
!  values of the derivative w.r.t. y of (1).
!
!  Input parameters:
!    n - parameter controlling the number of terms in the expansion (1)
!
!  Output parameters:
!    adiffx - the spectral differentiation matrix

double precision, allocatable :: adiff0(:,:), eye0(:,:)

allocate(eye0(n,n))

eye0 = 0
do i=1,n
eye0(i,i) = 1
end do

call chebyshev_diffmatrix(n,adiff0)
call bicheb2_kron(n,adiff0,eye0,adiffy)

end subroutine


subroutine bicheb2_kron(n,a,b,c)
implicit double precision (a-h,o-z)
double precision                            :: a(:,:), b(:,:)
double precision, allocatable, intent(out)  :: c(:,:)
!
!  Form the kronecker product of the (n,n) matrix a and b; return it in
!  the matrix c.
!
allocate(c(n**2,n**2))

l1 = 1
l2 = n

do i=1,n
k1 = 1
k2 = n
do j=1,n

c(l1:l2,k1:k2) = a(i,j) * b

k1 = k1+n
k2 = k2+n
end do
l1 = l1+n
l2 = l2+n
end do

end subroutine


end module 

