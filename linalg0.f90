!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains certain unoptimized linear algebraic routines.  Many of these
!  routines have two versions, one for real-valued inputs and a second version for
!  complex-valued inputs.
!
!  The raisons d'être for these routines:
!
!    1. They have no external dependencies.
!
!    2. They can easily make use of extended precision arithmetic.
!       
!    3. In the case of small-scale problems, they are often faster than equivalent
!       LAPACK routines owing to a lack of overhead.
!
!  Unless otherwise noted in a routine's documentation, it does NOT destroy input 
!  matrices and vectors.
!
!  The following subroutines are publicly callable:
!    
!    linalg0_solve - solve the system of equations A*x = b with A a well-conditioned
!      square matrix by forming a QR decomposition of the coefficient matrix A
!
!      THIS ROUTINE USES A NONPIVOTED ALGORITHM AND IS ONLY SUITABLE WHEN A IS
!      WELL-CONDITIONED, BUT IT IS GENERALLY FASTER THAN THE OTHER SOLVERS
!
!    linalg0_invert - form the inverse of a square matrix A using a QR decomposition
!
!      THIS ROUTINE USES A NONPIVOTED ALGORITHM AND IS ONLY SUITABLE WHEN A IS 
!      WELL-CONDITIONED, BUT IT IS GENERALLY FASTER THAN THE OTHER INVERTERS
!
!
!
!    linalg0_qrt - construct a rank-revealing QR factorization of the transpose/conjugate
!      of the matrix A; that is, factor it as
!
!        A^*(m,n) = Q(m,k) R(k,n)
!
!      where, Q has orthonormal columns and R is obtained from an upper trapezoidal
!      matrix via a column permutation (in fact, R(:,ipivs) is upper trapezoidal)
!
!    linalg0_qrtsolver - use a pre-exsiting QR decomposition of A^* to solve the system 
!      of equations AX = B with B either a vector or matrix
!
!    linalg0_qrtsolve - this is a convenience routine which forms a QR decomposition
!      of the transpose/cojugate of A and then solves the system AX = B using it without 
!      returning the products of the factorization routine
!
!    linalg0_qrtinverter - uisng a pre-existing QR decomposition of A^*, solve the
!      system of equations AX = I with I an appropriately-sized identity matrix
!
!    linalg0_qrtinvert - this is a convenience routine which forms the QR decomposition
!      of A^* and solves the system AX=I without returning to the user the products of 
!      the factorization routine
!
!
!
!    linalg0_qr - construct a rank-revealing QR factorization of the matrix A;
!      that is, factor it as
!
!        A(n,m) = Q(n,k) R(k,m)
!
!      where, Q has orthonormal columns and R is obtained from an upper trapezoidal
!      matrix via a column permutation (in fact, R(:,ipivs) is upper trapezoidal)
!
!    linalg0_qrsolver - using a pre-existing QR decomposition of A, solve the system of 
!      equations A*X = B with B either vector or matrix
!
!    linalg0_qrsolve - this is a convenience routine which forms a rank-reavling
!      QR decomposition using the qr routine and then solves the system AX = B 
!      for the user without returning the products of the factorization routine
!
!    linalg0_qrinverter - using a pre-existing QR decomposition of A, solve the system of 
!      equations AX = I with I an appropriately-sized identity matrix
!
!    linalg0_qrinvert - this is a convenience routine which forms the qr decomposition
!      of A and uses it to solve the system AX=I without returning to the user the 
!      products of the factorization routine
!
!
!
!    linalg0_utv - form a utv factorization of an input matrix A --- that is, produce
!      a factorization of the type
!
!        A(n,m) = U (n,k) T(k,k) V^*(k,m), 
!
!      where U and V have orthonormal columns and a ROW perumation of T is lower
!      triangular (i.e., T(ipivs,:) is lower triangular)
!
!    linalg0_utvsolver - using a pre-existing UTV decomposition of A, solve the system of 
!      equations A*X = B with B either a vector or matrix
!
!    linalg0_utvsolve - this is a convenience routine which forms a utv decompositionb
!      and then uses it to solve the system AX = B for the user without returning
!      the products of the factorization routine
!
!    linalg0_utvinverter - using a pre-existing UTV decomposition of A, solve the system of 
!      equations A*X = I with I an appropriately-sized identity matrix 
!
!    linalg0_utvinvert - this is a convenience routine which forms the utv decomposition
!      of A and solves the system AX=I without returning to the user the products of the
!      factorization routine
!
!
!      
!    linalg0_svd - form the singular value decomposition of an input matrix
!      A; that is, factor A as 
!
!        A(n,m) = U(n,l) D(l,l) V^* (l,m),
!
!      where l = min(n,m), the matrices U and V have orthonormal columns and D is
!      diagonal with positive real entries, sorted in descending order
!
!    linalg0_svdsolver - use a pre-existing full SVD factorization of a matrix A
!      to solve the system A*X = B with B either a matrix or a vector
!
!    linalg0_svdsolve - this is a convenience routine which forms the SVD and then
!      uses it to solve a system A*X = B without 
!
!
!
!    linalg0_jaceig - use Jacobi iterations to diagonalize a Hermitian or symmetric
!      matrix
!
!    linalg0_jacsvd - use Jacobi iterations to compute the singular value decomposition
!      of a square matrix
!
!    linalg0_jacqr - use Jacobi rotations to compute the QR decomposition of a 
!      rectangular matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!  Observations / To do:
!    
!    - Check to make sure that fortran is not making copies of objects on the stack
!      unnecessarily
!     
!    - Write nonsymmetric eigendecomposition/Schur decompostion code
!
!    - Write the svd_inverter and svd_invert routines
!
!    - Add the "truncated "SVD computed via UTV followed by Jacobi rotations
!     
!    - Eliminate some of the many inefficiencies in these codes, it would be
!      particularly useful to speed up the QRT routines
!
!    - Write a "linalg0_sings" or linalg0_jacsings" routine for computing only
!      the singular values of a rectangular matrix ... do it in place with
!      Givens rotations to zero the nonsquare portion followed by Jacobi iterations.
!
!    - Add interpolative decomposition routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module linalg0

use utils

interface        linalg0_solve
module procedure linalg0_solve_r
module procedure linalg0_solve_c
end interface    linalg0_solve

interface        linalg0_invert
module procedure linalg0_invert_r
module procedure linalg0_invert_c
end interface    linalg0_invert



interface        linalg0_qrt
module procedure linalg0_qrt_c
module procedure linalg0_qrt_r
end interface    linalg0_qrt

interface        linalg0_qrtsolver
module procedure linalg0_qrtsolver_r1
module procedure linalg0_qrtsolver_c1
module procedure linalg0_qrtsolver_r2
module procedure linalg0_qrtsolver_c2
end interface    linalg0_qrtsolver

interface        linalg0_qrtsolve
module procedure linalg0_qrtsolve_r1
module procedure linalg0_qrtsolve_c1
module procedure linalg0_qrtsolve_r2
module procedure linalg0_qrtsolve_c2
end interface    linalg0_qrtsolve

interface        linalg0_qrtinverter
module procedure linalg0_qrtinverter_r
module procedure linalg0_qrtinverter_c
end interface    linalg0_qrtinverter

interface        linalg0_qrtinvert
module procedure linalg0_qrtinvert_r
module procedure linalg0_qrtinvert_c
end interface    linalg0_qrtinvert


interface        linalg0_qr
module procedure linalg0_qr_c
module procedure linalg0_qr_r
end interface    linalg0_qr

interface        linalg0_qrsolver
module procedure linalg0_qrsolver_r1
module procedure linalg0_qrsolver_c1
module procedure linalg0_qrsolver_r2
module procedure linalg0_qrsolver_c2
end interface    linalg0_qrsolver

interface        linalg0_qrsolve
module procedure linalg0_qrsolve_r1
module procedure linalg0_qrsolve_c1
module procedure linalg0_qrsolve_r2
module procedure linalg0_qrsolve_c2
end interface    linalg0_qrsolve

interface        linalg0_qrinverter
module procedure linalg0_qrinverter_r
module procedure linalg0_qrinverter_c
end interface    linalg0_qrinverter

interface        linalg0_qrinvert
module procedure linalg0_qrinvert_r
module procedure linalg0_qrinvert_c
end interface    linalg0_qrinvert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


interface        linalg0_utv
module procedure linalg0_utv_r
module procedure linalg0_utv_c
end interface    linalg0_utv

interface        linalg0_utvsolve
module procedure linalg0_utvsolve_r1
module procedure linalg0_utvsolve_c1
module procedure linalg0_utvsolve_r2
module procedure linalg0_utvsolve_c2
end interface    linalg0_utvsolve

interface        linalg0_utvsolver
module procedure linalg0_utvsolver_r1
module procedure linalg0_utvsolver_c1
module procedure linalg0_utvsolver_r2
module procedure linalg0_utvsolver_c2
end interface    linalg0_utvsolver


interface        linalg0_utvinverter
module procedure linalg0_utvinverter_r
module procedure linalg0_utvinverter_c
end interface    linalg0_utvinverter

interface        linalg0_utvinvert
module procedure linalg0_utvinvert_r
module procedure linalg0_utvinvert_c
end interface    linalg0_utvinvert


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

interface        linalg0_svd
module procedure linalg0_svd_r
module procedure linalg0_svd_c
end interface    linalg0_svd

interface        linalg0_svdsolver
module procedure linalg0_svdsolver_r1
module procedure linalg0_svdsolver_c1
module procedure linalg0_svdsolver_r2
module procedure linalg0_svdsolver_c2
end interface    linalg0_svdsolver

interface        linalg0_svdsolve
module procedure linalg0_svdsolve_r1
module procedure linalg0_svdsolve_c1
module procedure linalg0_svdsolve_r2
module procedure linalg0_svdsolve_c2
end interface    linalg0_svdsolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


interface        linalg0_jacsvd
module procedure linalg0_jacsvd_r
module procedure linalg0_jacsvd_c
end interface    linalg0_jacsvd

interface        linalg0_jacqr
module procedure linalg0_jacqr_r
module procedure linalg0_jacqr_c
end interface    linalg0_jacqr

interface        linalg0_jaceig
module procedure linalg0_jaceig_r
module procedure linalg0_jaceig_c
end interface    linalg0_jaceig

interface        linalg0_gspiv
module procedure linalg0_gspiv_r
module procedure linalg0_gspiv_c
end interface    linalg0_gspiv

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Quick and dirty solve / invert routines 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg0_solve_r(n,a,rhs)
implicit double precision (a-h,o-z)
double precision   :: a(:,:),rhs(:)
!
!  This subroutine uses a version of QR-decomposition to solve the equation
!  A x = b.  
!
!  THE MATRIX A IS DESTROYED BY THIS ROUTINE.
!
!  Input parameters:
!    a - the (n,n) matrix of coefficients
!    n - an integer specifying the size of the system of 
!    rhs - a vector of length n speciying the rhs of the system
!
!  Output parameters:
!
!   rhs - upon return, the solution of the linear system


double precision :: aa(2),u(2,2)
! 
! transpose the input matrix a 
!

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)**2
size22=size22+a(i,j)**2
end do
end do

!
!  Reduce to upper triangular 
!
do i=1,n-1
do j=n,i+1,-1
aa(1)=a(i,j-1)
aa(2)=a(i,j)

u22=-aa(1)
u12=aa(2)
d=u22**2+u12**2
if(d .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
d=sqrt(d)
u(2,2)=u22/d
u(1,2)=u12/d
u(1,1)=-u(2,2)
u(2,1)=u(1,2)
endif

do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j) 
a(ii,j-1)=d1
a(ii,j)=d2
end do
d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
rhs(j-1)=d1
rhs(j)=d2
end do
end do

!
!  Apply the inverse of the triangular matrix
! 

rhs(n)=rhs(n)/a(n,n)
do i=n-1,1,-1
d=0
do j=n,i+1,-1
d=d+a(j,i)*rhs(j)
end do
rhs(i)=(rhs(i)-d)/a(i,i)
end do

return
end subroutine



subroutine linalg0_solve_c(n,a,rhs)
implicit complex *16 (a-h,o-z)
double complex         :: a(n,n),u(2,2),rhs(n)
! 
!  This subroutine uses a version of QR-decomposition to solve
!  the user-supplied system of linear algebraic equations
!  Ax = y.
!
!  THE MATRIX A IS DESTROYED BY THIS ROUTINE AND THE VECTOR RHS IS
!  OVERWRITTEN.
! 
!  Input parameters:
!    n - the dimensionality of the system being solved
!    a - the (n,n) complex-valued coefficient matrix WHICH WILL
!      BE DESTROYED BY THIS ROUUTINE
!    rhs - the right-hand side of the system to be solved, which
!      is overwritten by this subroutine
! 
!  Output parameters:
!    rhs - the solution of the linear system
!

double precision       ::  dmax,dmin,size22,dd,eps0
double complex         ::  aa(2)
integer, allocatable   ::  ipiv(:)

!
! Transpose the input matrix a
! 

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)*conjg(a(j,i))
size22=size22+a(i,j)*conjg(a(i,j))
end do
end do

! 
! Eliminate the upper right triangle
! 
do i=1,n-1
! 
do j=n,i+1,-1
! 
aa(1)=a(i,j-1)
aa(2)=a(i,j)
u21=-aa(2)
u22=aa(1)
dd=u22*conjg(u22)+u21*conjg(u21)
if(dd .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
dd=sqrt(dd)
u(2,2)=u22/dd
u(2,1)=u21/dd
u(1,1)=-conjg(u(2,2))
u(1,2)=conjg(u(2,1))
endif


do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j)
a(ii,j-1) = d1
a(ii,j)   = d2
end do

d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
rhs(j-1)=d1
rhs(j)=d2
end do
end do

!
!  Apply the inverse of the triangular matrix to rhs
!

rhs(n)=rhs(n)/a(n,n)
do i=n-1,1,-1
d=0
do j=n,i+1,-1
d=d+a(j,i)*rhs(j)
end do
rhs(i)=(rhs(i)-d)/a(i,i)
end do

end subroutine


subroutine linalg0_invert_r(n,a)
implicit double precision (a-h,o-z)
double precision a(:,:)
!
!  Invert the (n,n) real-valued matrix a using a QR decomposition.
!
!  THE MATRIX A IS OVERWRITTEN BY THIS PROCEDURE.
!
!  Input/output parameters:
!    a - upon input the (n,n) matrix to invert, which is replaced by this procedure 
!      with its inverse
!
double precision, allocatable :: b(:,:), work(:)
allocate(b(n,n), work(n))
done=1
zero=0

b = 0.0d0
do i=1,n
b(i,i) = 1.0d0
end do
        
! GS the input matrix
do i=1,n
do j=1,i-1
cd = dot_product(a(i,:),a(j,:))
a(i,:) = a(i,:) - a(j,:)*cd
b(i,:) = b(i,:) - b(j,:)*cd
end do

d = 1/norm2(a(i,:))
a(i,:) = a(i,:)*d
b(i,:) = b(i,:)*d
do j=i+1,n
cd = dot_product(a(i,:),a(j,:))
a(j,:) = a(j,:)-cd*a(i,:)
b(j,:) = b(j,:)-cd*b(i,:)
end do
end do
!
!  Multiply the adjoint of the resulting orthogonal matrix by the
!  triangular one to invert a  
!
do i=1,n
do j=1,n
cd = dot_product(a(:,i),b(:,j))
work(j)=cd
end do
a(:,i) = work
end do
!
!  Transpose A
!
b = transpose(a)
a = b

return
end subroutine



subroutine linalg0_invert_c(n,a)
implicit double precision (a-h,o-z)
double complex :: a(:,:)
!
!  Invert the (n,n) complex-valued matrix a using a QR decomposition.
!
!  THE MATRIX A IS OVERWRITTEN BY THIS PROCEDURE.
!
!  Input/output parameters:
!    a - upon input the (n,n) matrix to invert, which is replaced by this 
!      procedure  with its inverse
!
double complex, allocatable :: b(:,:), work(:)
double complex              :: cd

allocate(b(n,n), work(n))
done=1
zero=0

b = 0.0d0
do i=1,n
b(i,i) = 1.0d0
end do
        
! GS the input matrix
do i=1,n
do j=1,i-1
cd = dot_product(a(i,:),a(j,:))
a(i,:) = a(i,:) - a(j,:)*cd
b(i,:) = b(i,:) - b(j,:)*cd
end do

d = 1/norm2(abs(a(i,:)))
a(i,:) = a(i,:)*d
b(i,:) = b(i,:)*d
do j=i+1,n
cd = dot_product(a(i,:),a(j,:))
a(j,:) = a(j,:)-cd*a(i,:)
b(j,:) = b(j,:)-cd*b(i,:)
end do
end do
!
!  Multiply the adjoint of the resulting orthogonal matrix by the
!  triangular one to invert a  
!
do i=1,n
do j=1,n
cd = dot_product(a(:,i),b(:,j))
work(j)=cd
end do
a(:,i) = work
end do
!
!  Transpose A
!
b = transpose(a)
a = b

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  QRT Factorization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg0_qrt_r(eps,n,m,a,krank,rnorms,ipivs,q,r)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
integer, allocatable, intent(out)          :: ipivs(:)
double precision, allocatable, intent(out) :: rnorms(:)
double precision, allocatable, intent(out) :: q(:,:)
double precision, allocatable, intent(out) :: r(:,:)
!
!  Form a rank-revealing QR decomposition of the transpose of the input matrix A. 
!  That is, factor A^* as
! 
!     A^*(m,n) = Q(m,krank) * R(krank,n)
!
!  where Q has orthonormal columns and R(:,ipivs) is upper trapezoidal.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    krank - the apparent numerical rank of the input matrix
!    rnorms - an array of length m specifying the normalizing factors
!      obtained from the GS procedure
!    ipivs - an integer array of length m specifying the list of pivots
!      from the GS procedure
!    q - the matrix q in the factorization
!    r - the matrix r in the factorization
!

double precision, allocatable   :: at(:,:)

allocate(at(m,n))
at = transpose(a)
call linalg0_qr(eps,m,n,at,krank,rnorms,ipivs,q,r)
return

! allocate(at(m,n))
! at = transpose(a)

! allocate(ipivs(n),rnorms(n))
! call linalg0_gspiv_r(eps,m,n,at,krank,rnorms,ipivs)

! allocate(q(m,krank), r(krank,n))
! q  = at(:,1:krank)
! !r  = matmul(transpose(q),transpose(a))

end subroutine



subroutine linalg0_qrt_c(eps,n,m,a,krank,rnorms,ipivs,q,r)
implicit double precision (a-h,o-z)
double complex                             :: a(:,:)
integer, allocatable, intent(out)          :: ipivs(:)
double precision, allocatable, intent(out) :: rnorms(:)
double complex, allocatable, intent(out)   :: q(:,:)
double complex, allocatable, intent(out)   :: r(:,:)
!
!  Form a rank-revealing QR decomposition of the conjugate of the input matrix A. 
!  That is, factor A as
! 
!     A^*(n,m) = Q(m,krank) * R(krank,n)
!
!  where Q has orthonormal columns and R(:,ipivs) is upper trapezoidal.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    krank - the apparent numerical rank of the input matrix
!    rnorms - an array of length m specifying the normalizing factors
!      obtained from the GS procedure
!    ipivs - an integer array of length m specifying the list of pivots
!      from the GS procedure
!    q - the matrix q in the factorization
!    r - the matrix r in the factorization
!
!

double complex, allocatable   :: at(:,:)


allocate(at(m,n))
at = conjg(transpose(a))
call linalg0_qr(eps,m,n,at,krank,rnorms,ipivs,q,r)
return

! ! Make a copy of the input array.
! allocate(at(m,n))
! at = transpose(conjg(a))

! allocate(ipivs(n),rnorms(n))
! call linalg0_gspiv_c(eps,m,n,at,krank,rnorms,ipivs)
! allocate(q(m,krank), r(krank,n))
! q  = at(:,1:krank)

! !!! THIS SEEMS TO BE REQUIRED TO PREVENT A COPY OPERATION, AT LEAST ON
! !!! SOME COMPILERS
! at = conjg(transpose(a))
! r  = matmul(conjg(transpose(q)),at)

end subroutine


subroutine linalg0_qrtsolver_r1(n,m,q,r,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double precision                 :: q(:,:), r(:,:), x(:), y(:)
integer                          :: ipivs(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  using a QR decomposition of A^* formed using the linalg0_qrt routine.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    q, r - the factors appearing in the decomposition of A^* formed by 
!     the qrt routine
!    krank - the approximate numerical rank returned by the qrt routine
!    ipivs - the list of pivots returned by the qrt routine
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!
!
double precision, allocatable :: z(:)

!
!  Solve the system R^* z = y 
!

allocate(z(krank))
do i=1,krank
z(i) = y(ipivs(i))
do j=1,i-1
z(i) = z(i) - z(j) * r(j,ipivs(i))
end do
z(i) = z(i) / r(i,ipivs(i))
end do

x = matmul(q,z)

end subroutine


subroutine linalg0_qrtsolver_c1(n,m,q,r,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double complex                :: q(:,:), r(:,:), x(:), y(:)
integer                       :: ipivs(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  using a QR decomposition of A^* formed using the linalg0_qrt routine.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    q, r - the factors appearing in the decomposition of A^* formed by 
!     the qrt routine
!    krank - the approximate numerical rank returned by the qrt routine
!    ipivs - the list of pivots returned by the qrt routine
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!

double complex, allocatable   :: z(:)


!
!  Solve the system R^* z = y 
!

allocate(z(krank))
do i=1,krank
z(i) = y(ipivs(i))
do j=1,i-1
z(i) = z(i) - z(j) * conjg(r(j,ipivs(i)))
end do
z(i) = z(i) / conjg(r(i,ipivs(i)))
end do

x = matmul(q,z)

end subroutine

subroutine linalg0_qrtsolver_r2(n,m,q,r,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double precision              :: q(:,:), r(:,:), x(:,:), y(:,:)
integer                       :: ipivs(:)
!
!  Solve the linear system of equations
!
!    A X = Y
!
!  using a QR decomposition of A^* formed using the linalg0_qrt routine.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    q, r - the factors appearing in the decomposition of A^* formed by 
!     the qrt routine
!    krank - the approximate numerical rank returned by the qrt routine
!    ipivs - the list of pivots returned by the qrt routine
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!

double precision, allocatable   :: z(:,:)

l = size(y,2)

!
!  Solve the system R^* z = y 
!

allocate(z(krank,l))
do i=1,krank
z(i,:) = y(ipivs(i),:)
do j=1,i-1
z(i,:) = z(i,:) - z(j,:) * r(j,ipivs(i))
end do
z(i,:) = z(i,:) / r(i,ipivs(i))
end do


x = matmul(q,z)

end subroutine

subroutine linalg0_qrtsolver_c2(n,m,q,r,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double complex                :: q(:,:), r(:,:), x(:,:), y(:,:)
integer                       :: ipivs(:)
!
!  Solve the linear system of equations
!
!    A X = Y
!
!  using a QR decomposition of A^* formed using the linalg0_qrt routine.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    q, r - the factors appearing in the decomposition of A^* formed by 
!     the qrt routine
!    krank - the approximate numerical rank returned by the qrt routine
!    ipivs - the list of pivots returned by the qrt routine
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!

double complex, allocatable   :: z(:,:)

l = size(y,2)

!
!  Solve the system R^* z = y 
!

allocate(z(krank,l))
do i=1,krank
z(i,:) = y(ipivs(i),:)
do j=1,i-1
z(i,:) = z(i,:) - z(j,:) * conjg(r(j,ipivs(i)))
end do
z(i,:) = z(i,:) / conjg(r(i,ipivs(i)))
end do

x = matmul(q,z)

end subroutine


subroutine linalg0_qrtsolve_r1(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double precision                 :: a(:,:), x(:), y(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  with y a vector by forming a rank-revealing QR decomposition of the 
!  transpose of A.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!
!
double precision, allocatable :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qrt(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrtsolver_r1(n,m,q,r,krank,ipivs,x,y)

end subroutine


subroutine linalg0_qrtsolve_c1(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double complex                  :: a(:,:), x(:), y(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  with y a vector by forming a rank-revealing QR decomposition of the transpose of A.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!
!
double complex, allocatable   :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qrt(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrtsolver(n,m,q,r,krank,ipivs,x,y)

end subroutine

subroutine linalg0_qrtsolve_r2(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double precision                 :: a(:,:), x(:,:), y(:,:)
!
!  Solve the linear system of equations
!
!    A X = Y 
!
!  with Y a matrix by forming a rank-revealing QR decomposition of the transpose of A.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!
!
double precision, allocatable :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qrt(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrtsolver(n,m,q,r,krank,ipivs,x,y)

end subroutine


subroutine linalg0_qrtsolve_c2(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double complex                  :: a(:,:), x(:,:), y(:,:)
!
!  Solve the linear system of equations
!
!    A X = Y
!
!  with Y a matrix by forming a rank-revealing QR decomposition of the transpose of A.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!
!
double complex, allocatable   :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qrt(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrtsolver(n,m,q,r,krank,ipivs,x,y)

end subroutine


subroutine linalg0_qrtinverter_r(n,m,q,r,krank,ipivs,ainv)
implicit double precision (a-h,o-z)
double precision                             :: q(:,:), r(:,:), ainv(:,:)
integer                                      :: ipivs(:)
!
!  Using the QR decomposition of the transpose of A, solve the linear system
!
!    A X = I.
!
!  Input parameters:
!    (n,m) - the dimensions of the input matrix
!    q,r - the factors in the QR decomposition of A^*
!    krank, ipivs - the numerical rank and list of pivots returned by the
!      qrt routine
!
!  Output parameters:
!    ainv - the solution of the system
!
!
double precision, allocatable :: z(:,:)


!
!  Solve the system R^* z = I
!

allocate(z(krank,n))
z = 0
do i=1,krank
z(i,ipivs(i)) = 1
do j=1,i-1
z(i,:) = z(i,:) - z(j,:) * r(j,ipivs(i))
end do
z(i,:) = z(i,:) / r(i,ipivs(i))
end do

!
!  Form the solution
!

ainv = matmul(q,z)

end subroutine


subroutine linalg0_qrtinverter_c(n,m,q,r,krank,ipivs,ainv)
implicit double precision (a-h,o-z)
double complex                             :: q(:,:), r(:,:), ainv(:,:)
integer                                    :: ipivs(:)
!
!  Using the QR decomposition of the transpose of A, solve the linear system
!
!    A X = I.
!
!  Input parameters:
!    (n,m) - the dimensions of the input matrix
!    q,r - the factors in the QR decomposition of A^*
!    krank, ipivs - the numerical rank and list of pivots returned by the
!      qrt routine
!
!  Output parameters:
!    ainv - the solution of the system
!
!

double complex, allocatable :: z(:,:)

!
!  Solve the system R^* z = I
!

allocate(z(krank,n))
z = 0
do i=1,krank
z(i,ipivs(i)) = 1
do j=1,i-1
z(i,:) = z(i,:) - z(j,:) * conjg(r(j,ipivs(i)))
end do
z(i,:) = z(i,:) / conjg(r(i,ipivs(i)))
end do

!
!  Form the solution
!

ainv = matmul(q,z)


end subroutine



subroutine linalg0_qrtinvert_r(eps,n,m,a,ainv)
implicit double precision (a-h,o-z)
double precision                             :: a(:,:), ainv(:,:)
!
! Form a QR decomposition of the transpose of A and use it to solve the linear system
!
!    A X = I.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimensions of the input matrix
!    q,r - the factors in the QR decomposition of A^*
!    krank, ipivs - the numerical rank and list of pivots returned by the
!      qrt routine
!
!  Output parameters:
!    ainv - the solution of the system
!
!

double precision, allocatable   :: q(:,:), r(:,:)
double precision, allocatable   :: rnorms(:)
integer, allocatable            :: ipivs(:)

call linalg0_qrt(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrtinverter(n,m,q,r,krank,ipivs,ainv)

end subroutine


subroutine linalg0_qrtinvert_c(eps,n,m,a,ainv)
implicit double precision (a-h,o-z)
double complex                             :: a(:,:), ainv(:,:)
!
!  Form the QR decommposition of the transpose of A and use it to solve the system.
!
!    A X = I.
!
!  Input parameters:
!    eps - precision for the factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ainv - the solution of the system
!
!
double complex, allocatable   :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qrt(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrtinverter(n,m,q,r,krank,ipivs,ainv)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  QR Factorization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg0_qr_r(eps,n,m,a,krank,rnorms,ipivs,q,r)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
integer, allocatable, intent(out)          :: ipivs(:)
double precision, allocatable, intent(out) :: rnorms(:)
double precision, allocatable, intent(out) :: q(:,:)
double precision, allocatable, intent(out) :: r(:,:)
!
!  Form a rank-revealing QR decomposition of the input matrix A; that is, factor
!  A
! 
!     A(m,n) = Q(n,krank) * R(krank,m)
!
!  where Q has orthonormal columns and R(:,ipivs) is upper trapezoidal.
!
!  Input parameters:
!    eps - the desired precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    krank - the apparent numerical rank of the input matrix
!    rnorms - an array of length m specifying the normalizing factors
!      obtained from the GS procedure
!    ipivs - an integer array of length m specifying the list of pivots
!      from the GS procedure
!    q - the matrix q in the factorization
!    r - the matrix r in the factorization
!

double precision, allocatable   :: b(:,:)

allocate(b(n,m))
b = a

allocate(ipivs(m),rnorms(m))
call linalg0_gspiv_r(eps,n,m,b,krank,rnorms,ipivs)

allocate(q(n,krank), r(krank,m))
q  = b(:,1:krank)
r  = matmul(transpose(q),a)

end subroutine


subroutine linalg0_qr_c(eps,n,m,a,krank,rnorms,ipivs,q,r)
implicit double precision (a-h,o-z)
double complex                             :: a(:,:)
integer, allocatable, intent(out)          :: ipivs(:)
double precision, allocatable, intent(out) :: rnorms(:)
double complex, allocatable, intent(out)   :: q(:,:)
double complex, allocatable, intent(out)   :: r(:,:)
!
!  Form a rank-revealing QR decomposition of the input matrix A
! 
!     A(n,m) = Q(n,krank) * R(krank,m)
!
!  where Q has orthonormal columns and R(:,ipivs) is upper trapezoidal.
!
!  Input parameters:
!    eps - the desired precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    krank - the apparent numerical rank of the input matrix
!    rnorms - an array of length m specifying the normalizing factors
!      obtained from the GS procedure
!    ipivs - an integer array of length m specifying the list of pivots
!      from the GS procedure
!    q - the matrix q in the factorization
!    r - the matrix r in the factorization
!
!

double complex, allocatable   :: b(:,:)


! Make a copy of the input array.

allocate(b(n,m))
b = a

allocate(ipivs(m),rnorms(m))
call linalg0_gspiv_c(eps,n,m,b,krank,rnorms,ipivs)

allocate(q(n,krank), r(krank,m))
q  = b(:,1:krank)
r  = matmul(conjg(transpose(q)),a)

end subroutine


subroutine linalg0_qrsolver_r1(n,m,q,r,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double precision              :: q(:,:), r(:,:), x(:), y(:)
integer                       :: ipivs(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  using a QR decomposition of A formed using the linalg0_qrt routine.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    q, r - the factors appearing in the decomposition of A^* formed by 
!     the qrt routine
!    krank - the approximate numerical rank returned by the qrt routine
!    ipivs - the list of pivots returned by the qrt routine
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!


double precision, allocatable   :: z(:)


!
!  Form the vector z = Q^*y
!

allocate(z(krank))
z = matmul(transpose(q),y)


!
!  Solve the system R x = z with  R(:,ipivs) upper trapezoidal
!

do i=krank,1,-1
do j=krank,i+1,-1
z(i) = z(i) - z(j) * r(i,ipivs(j))
end do
z(i) = z(i) / r(i,ipivs(i))
end do


x(1:m)            = 0
x(ipivs(1:krank)) = z(1:krank)

end subroutine


subroutine linalg0_qrsolver_c1(n,m,q,r,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double complex                :: q(:,:), r(:,:), x(:), y(:)
integer                       :: ipivs(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  using a QR decomposition of A formed using the linalg0_qrt routine.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    q, r - the factors appearing in the decomposition of A^* formed by 
!     the qrt routine
!    krank - the approximate numerical rank returned by the qrt routine
!    ipivs - the list of pivots returned by the qrt routine
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!


double complex, allocatable   :: z(:)


!
!  Form the vector z = Q^*y
!

allocate(z(krank))
z = matmul(conjg(transpose(q)),y)

!
!  Solve the system R x = z with  R(:,ipivs) upper trapezoidal
!

do i=krank,1,-1
do j=krank,i+1,-1
z(i) = z(i) - z(j) * r(i,ipivs(j))
end do
z(i) = z(i) / r(i,ipivs(i))
end do

x(1:m)            = 0
x(ipivs(1:krank)) = z(1:krank)

end subroutine



subroutine linalg0_qrsolver_r2(n,m,q,r,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double precision              :: q(:,:), r(:,:), x(:,:), y(:,:)
integer                       :: ipivs(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  using a QR decomposition of A formed using the linalg0_qrt routine.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    q, r - the factors appearing in the decomposition of A^* formed by 
!     the qrt routine
!    krank - the approximate numerical rank returned by the qrt routine
!    ipivs - the list of pivots returned by the qrt routine
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!


double precision, allocatable   :: z(:,:)

l = size(y,2)

!
!  Form the vector z = Q^*y
!

allocate(z(krank,l))
z = matmul(transpose(q),y)

! allocate(z2(krank,l))
! z2 = z


!
!  Solve the system R x = z with  R(:,ipivs) upper trapezoidal
!

do i=krank,1,-1
do j=krank,i+1,-1
z(i,:) = z(i,:) - z(j,:) * r(i,ipivs(j))
end do
z(i,:) = z(i,:) / r(i,ipivs(i))
end do


x(1:m,1:l)            = 0
x(ipivs(1:krank),1:l) = z(1:krank,:)

end subroutine


subroutine linalg0_qrsolver_c2(n,m,q,r,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double complex                :: q(:,:), r(:,:), x(:,:), y(:,:)
integer                       :: ipivs(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  using a QR decomposition of A formed using the linalg0_qrt routine.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    q, r - the factors appearing in the decomposition of A^* formed by 
!     the qrt routine
!    krank - the approximate numerical rank returned by the qrt routine
!    ipivs - the list of pivots returned by the qrt routine
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!


double complex, allocatable   :: z(:,:)


l = size(y,2)

!
!  Form the vector z = Q^*y
!

allocate(z(krank,l))
z = matmul(conjg(transpose(q)),y)



!
!  Solve the system R x = z with  R(:,ipivs) upper trapezoidal
!

do i=krank,1,-1
do j=krank,i+1,-1
z(i,:) = z(i,:) - z(j,:) * r(i,ipivs(j))
end do
z(i,:) = z(i,:) / r(i,ipivs(i))
end do


x(1:m,1:l)            = 0
x(ipivs(1:krank),1:l) = z(1:krank,:)

end subroutine



subroutine linalg0_qrsolve_r1(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double precision                 :: a(:,:), x(:), y(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  with y a vector by forming a rank-revealing QR decomposition of the 
!  input matrix A.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!
!
double precision, allocatable :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qr(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrsolver(n,m,q,r,krank,ipivs,x,y)

end subroutine


subroutine linalg0_qrsolve_c1(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double complex                  :: a(:,:), x(:), y(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  with y a vector by forming a rank-revealing QR decomposition of the matrix A.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!
!
double complex, allocatable   :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qr(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrsolver(n,m,q,r,krank,ipivs,x,y)

end subroutine

subroutine linalg0_qrsolve_r2(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double precision                 :: a(:,:), x(:,:), y(:,:)
!
!  Solve the linear system of equations
!
!    A X = Y 
!
!  with Y a matrix by forming a rank-revealing QR decomposition of the matrix A.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!
!
double precision, allocatable :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qr(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrsolver(n,m,q,r,krank,ipivs,x,y)

end subroutine


subroutine linalg0_qrsolve_c2(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double complex                  :: a(:,:), x(:,:), y(:,:)
!
!  Solve the linear system of equations
!
!    A X = Y
!
!  with Y a matrix by forming a rank-revealing QR decomposition of the input
!  matrix A.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the solution of the system
!
!
double complex, allocatable   :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qr(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrsolver(n,m,q,r,krank,ipivs,x,y)

end subroutine


subroutine linalg0_qrinverter_r(n,m,q,r,krank,ipivs,ainv)
implicit double precision (a-h,o-z)
double precision                             :: q(:,:), r(:,:), ainv(:,:)
integer                                      :: ipivs(:)
!
!  Using the QR decomposition of A, solve the linear system
!
!    A X = I.
!
!  Input parameters:
!    (n,m) - the dimensions of the input matrix
!    q,r - the factors in the QR decomposition of A^*
!    krank, ipivs - the numerical rank and list of pivots returned by the
!      qrt routine
!
!  Output parameters:
!    ainv - the solution of the system
!
!
double precision, allocatable :: z(:,:)



allocate(z(krank,n))
z = transpose(q)

!
!  Solve the system R x = z with  R(:,ipivs) upper trapezoidal
!

do i=krank,1,-1
do j=krank,i+1,-1
z(i,:) = z(i,:) - z(j,:) * r(i,ipivs(j))
end do
z(i,:) = z(i,:) / r(i,ipivs(i))
end do


ainv(1:m,1:n)            = 0
ainv(ipivs(1:krank),1:n) = z(1:krank,:)

end subroutine


subroutine linalg0_qrinverter_c(n,m,q,r,krank,ipivs,ainv)
implicit double precision (a-h,o-z)
double complex                             :: q(:,:), r(:,:), ainv(:,:)
integer                                    :: ipivs(:)
!
!  Using a QR decomposition of the matrix A, solve the linear system
!
!    A X = I.
!
!  Input parameters:
!    (n,m) - the dimensions of the input matrix
!    q,r - the factors in the QR decomposition of A^*
!    krank, ipivs - the numerical rank and list of pivots returned by the
!      qrt routine
!
!  Output parameters:
!    ainv - the solution of the system
!
!

double complex, allocatable :: z(:,:)


allocate(z(krank,n))
z = conjg(transpose(q))

!
!  Solve the system R x = z with  R(:,ipivs) upper trapezoidal
!

do i=krank,1,-1
do j=krank,i+1,-1
z(i,:) = z(i,:) - z(j,:) * r(i,ipivs(j))
end do
z(i,:) = z(i,:) / r(i,ipivs(i))
end do


ainv(1:m,1:n)            = 0
ainv(ipivs(1:krank),1:n) = z(1:krank,:)

end subroutine



subroutine linalg0_qrinvert_r(eps,n,m,a,ainv)
implicit double precision (a-h,o-z)
double precision                             :: a(:,:), ainv(:,:)
!
! Form a QR decomposition of A and use it to solve the linear system
!
!    A X = I.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimensions of the input matrix
!    q,r - the factors in the QR decomposition of A^*
!    krank, ipivs - the numerical rank and list of pivots returned by the
!      qrt routine
!
!  Output parameters:
!    ainv - the solution of the system
!
!

double precision, allocatable   :: q(:,:), r(:,:)
double precision, allocatable   :: rnorms(:)
integer, allocatable            :: ipivs(:)

call linalg0_qr(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrinverter(n,m,q,r,krank,ipivs,ainv)

end subroutine


subroutine linalg0_qrinvert_c(eps,n,m,a,ainv)
implicit double precision (a-h,o-z)
double complex                             :: a(:,:), ainv(:,:)
!
!  Form a QR decommposition of A and use it to solve the system.
!
!    A X = I.
!
!  Input parameters:
!    eps - precision for the factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ainv - the solution of the system
!
!
double complex, allocatable   :: q(:,:), r(:,:)
double precision, allocatable :: rnorms(:)
integer, allocatable          :: ipivs(:)

call linalg0_qr(eps,n,m,a,krank,rnorms,ipivs,q,r)
call linalg0_qrinverter(n,m,q,r,krank,ipivs,ainv)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  UTV Factorization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg0_utv_r(eps,n,m,a,krank,rnorms,ipivs,u,t,v)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
double precision, allocatable, intent(out) :: u(:,:), t(:,:), v(:,:)
double precision, allocatable, intent(out) :: rnorms(:)
integer, allocatable, intent(out)          :: ipivs(:)
!
!  Factor an input matrix A as
!
!    A(n,m) = U(n,krank) T(krank,krank) V(krank,m)^t
!
!  where the columns of U and V are orthonormal and T(ipivs,:) is 
!  lower triangular.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    krank - the apparent numerical rank of the input matrix
!    rnorms - an array of length m specifying the normalizing factors
!      obtained from the GS procedure
!    ipivs - an integer array of length m specifying the list of pivots
!      from the GS procedure
!    u, t, v - the matrices i nthe factorization
!
!

double precision, allocatable  :: r(:,:), tt(:,:), rt(:,:)
double precision, allocatable  :: rnorms2(:)

call linalg0_qr_r(eps,n,m,a,krank,rnorms,ipivs,u,r)

allocate(rt(m,krank))
rt = transpose(r)
call linalg0_qr3_r(m,krank,rt,rnorms2,ipivs,v,tt)
t = transpose(tt)

end subroutine


subroutine linalg0_utv_c(eps,n,m,a,krank,rnorms,ipivs,u,t,v)
implicit double precision (a-h,o-z)
double complex                             :: a(:,:)
double complex, allocatable, intent(out)   :: u(:,:), t(:,:), v(:,:)
double precision, allocatable, intent(out) :: rnorms(:)
integer, allocatable, intent(out)          :: ipivs(:)
!
!  Factor an input matrix A as
!
!    A(n,m) = U(n,krank) T(krank,krank) V(krank,m)^*
!
!  where the columns of U and V are orthonormal and T(ipivs,:) is 
!  lower triangular.
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    krank - the apparent numerical rank of the input matrix
!    rnorms - an array of length m specifying the normalizing factors
!      obtained from the GS procedure
!    ipivs - an integer array of length m specifying the list of pivots
!      from the GS procedure
!    u, t, v - the matrices i nthe factorization
!
!

double complex, allocatable    :: r(:,:), tt(:,:), rt(:,:)
double precision, allocatable  :: rnorms2(:)

call linalg0_qr_c(eps,n,m,a,krank,rnorms,ipivs,u,r)

allocate(rt(m,krank))
rt = conjg(transpose(r))
call linalg0_qr3_c(m,krank,rt,rnorms2,ipivs,v,tt)
t = conjg(transpose(tt))


end subroutine



subroutine linalg0_utvsolver_r1(n,m,u,t,v,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double precision        :: u(:,:), t(:,:), v(:,:), x(:), y(:)
integer                 :: ipivs(:)
!
!  Use a preexisting UTV decomposition of a to solve the system
!
!    A(n,m) * X(m) = B(n)
!
!  for X.
!
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double precision, allocatable :: z(:), z2(:)

!
!  Form the vector z = U* v
!
allocate(z(krank),z2(krank))
z  = matmul(transpose(u),y)

do i=1,krank
z2(i) = z(ipivs(i))
do j=1,i-1
z2(i) = z2(i) - z2(j) * t(ipivs(i),j)
end do
z2(i) = z2(i) / t(ipivs(i),i)
end do

!
!  Form the solution Vy
!
x = matmul(v,z2)

end subroutine


subroutine linalg0_utvsolver_c1(n,m,u,t,v,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double complex         :: u(:,:), t(:,:), v(:,:), x(:), y(:)
integer                 :: ipivs(:)
!
!  Use a preexisting UTV decomposition of a to solve the system
!
!    A(n,m) * X(m) = B(n)
!
!  for X.
!
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double complex, allocatable :: z(:), z2(:)

!
!  Form the vector z = U* v
!
allocate(z(krank),z2(krank))
z  = matmul(conjg(transpose(u)),y)

do i=1,krank
z2(i) = z(ipivs(i))
do j=1,i-1
z2(i) = z2(i) - z2(j) * t(ipivs(i),j)
end do
z2(i) = z2(i) / t(ipivs(i),i)
end do

!
!  Form the solution Vy
!
x = matmul(v,z2)

end subroutine


subroutine linalg0_utvsolver_r2(n,m,u,t,v,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double precision        :: u(:,:), t(:,:), v(:,:), x(:,:), y(:,:)
integer                 :: ipivs(:)
!
!  Use a preexisting UTV decomposition of a to solve the system
!
!    A(n,m) * X(m) = B(n)
!
!  for X.
!
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double precision, allocatable :: z(:,:), z2(:,:)

l = size(y,2)
!
!  Form the vector z = U* v
!
allocate(z(krank,l),z2(krank,l))
z  = matmul(transpose(u),y)

do i=1,krank
z2(i,:) = z(ipivs(i),:)
do j=1,i-1
z2(i,:) = z2(i,:) - z2(j,:) * t(ipivs(i),j)
end do
z2(i,:) = z2(i,:) / t(ipivs(i),i)
end do

!
!  Form the solution Vy
!
x = matmul(v,z2)

end subroutine


subroutine linalg0_utvsolver_c2(n,m,u,t,v,krank,ipivs,x,y)
implicit double precision (a-h,o-z)
double complex         :: u(:,:), t(:,:), v(:,:), x(:,:), y(:,:)
integer                 :: ipivs(:)
!
!  Use a preexisting UTV decomposition of a to solve the system
!
!    A(n,m) * X(m) = B(n)
!
!  for X.
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double complex, allocatable :: z(:,:), z2(:,:)

l = size(y,2)

!
!  Form the vector z = U* v
!
allocate(z(krank,l),z2(krank,l))
z  = matmul(conjg(transpose(u)),y)

do i=1,krank
z2(i,:) = z(ipivs(i),:)
do j=1,i-1
z2(i,:) = z2(i,:) - z2(j,:) * t(ipivs(i),j)
end do
z2(i,:) = z2(i,:) / t(ipivs(i),i)
end do

!
!  Form the solution Vy
!
x = matmul(v,z2)

end subroutine



subroutine linalg0_utvsolve_r1(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double precision        :: a(:,:), x(:), y(:)
!
!  Solve the equation
!
!    A(n,m) * X(m) = B(n)
!
!  in a least squares sense for X using a UTV decomposition.
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double precision, allocatable :: u(:,:), v(:,:), t(:,:)
integer, allocatable          :: ipivs(:)
double precision, allocatable :: rnorms(:)

!
!  Form the U T V decomposition
!

call linalg0_utv(eps,n,m,a,krank,rnorms,ipivs,u,t,v)
call linalg0_utvsolver(n,m,u,t,v,krank,ipivs,x,y)

end subroutine


subroutine linalg0_utvsolve_c1(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double complex        :: a(:,:), x(:), y(:)
!
!  Solve the equation
!
!    A(n,m) * X(m) = B(n)
!
!  in a least squares sense for X using a UTV decomposition.
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double complex, allocatable   :: u(:,:), v(:,:), t(:,:)
integer, allocatable          :: ipivs(:)
double precision, allocatable :: rnorms(:)
double complex                :: cd

!
!  Form the U T V decomposition
!

call linalg0_utv_c(eps,n,m,a,krank,rnorms,ipivs,u,t,v)
call linalg0_utvsolver(n,m,u,t,v,krank,ipivs,x,y)

end subroutine


subroutine linalg0_utvsolve_r2(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double precision        :: a(:,:), x(:,:), y(:,:)
!
!  Solve the equation
!
!    A(n,m) * X(m) = B(n)
!
!  in a least squares sense for X using a UTV decomposition.
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double precision, allocatable   :: u(:,:), v(:,:), t(:,:)
integer, allocatable            :: ipivs(:)
double precision, allocatable   :: rnorms(:)


!
!  Form the U T V decomposition
!

call linalg0_utv_r(eps,n,m,a,krank,rnorms,ipivs,u,t,v)
call linalg0_utvsolver(n,m,u,t,v,krank,ipivs,x,y)

end subroutine


subroutine linalg0_utvsolve_c2(eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double complex        :: a(:,:), x(:,:), y(:,:)
!
!  Solve the equation
!
!    A(n,m) * X(m) = B(n)
!
!  in a least squares sense for X using a UTV decomposition.
!
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double complex, allocatable   :: u(:,:), v(:,:), t(:,:)
integer, allocatable          :: ipivs(:)
double precision, allocatable :: rnorms(:)


call linalg0_utv_c(eps,n,m,a,krank,rnorms,ipivs,u,t,v)
call linalg0_utvsolver(n,m,u,t,v,krank,ipivs,x,y)

end subroutine


subroutine linalg0_utvinverter_r(n,m,u,t,v,krank,ipivs,ainv)
implicit double precision (a-h,o-z)
double precision             :: ainv(:,:), u(:,:), t(:,:), v(:,:)
integer                      :: ipivs(:)
!
!  Using an existing UTV decomposition,  solve the system
!
!    A(n,m) * X(m) = I(n,n),
!
!  where I is the (n,n) identity.
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ainv - the pseudoinverse
!

double precision, allocatable   :: z(:,:), z2(:,:)

!
!  Form the vector z = U* v
!
allocate(z(krank,m), z2(krank,m))
z  = transpose(u)
z2 = z(ipivs(1:krank),:)

do i=1,krank
do j=1,i-1
z2(i,:) = z2(i,:) - z2(j,:) * t(ipivs(i),j)
end do
z2(i,:) = z2(i,:) / t(ipivs(i),i)
end do

!
!  Form the solution Vy
!
ainv = matmul(v,z2)

end subroutine

subroutine linalg0_utvinverter_c(n,m,u,t,v,krank,ipivs,ainv)
implicit double precision (a-h,o-z)
double complex               :: ainv(:,:), u(:,:), t(:,:), v(:,:)
integer                      :: ipivs(:)
!
!  Using an existing UTV decomposition,  solve the system
!
!    A(n,m) * X(m) = I(n,n),
!
!  where I is the (n,n) identity.
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ainv - the pseudoinverse
!

double complex, allocatable   :: z(:,:), z2(:,:)


!
!  Form the vector z = U* v
!
allocate(z(krank,m), z2(krank,m))
z  = conjg(transpose(u))
z2 = z(ipivs(1:krank),:)

do i=1,krank
do j=1,i-1
z2(i,:) = z2(i,:) - z2(j,:) * t(ipivs(i),j)
end do
z2(i,:) = z2(i,:) / t(ipivs(i),i)
end do

!
!  Form the solution Vy
!
ainv = matmul(v,z2)

end subroutine

subroutine linalg0_utvinvert_r(eps,n,m,a,ainv)
implicit double precision (a-h,o-z)
double precision              :: a(:,:), ainv(:,:)
!
!  Solve the equation
!
!    A(n,m) * X(m) = I(n,n),
!
!  where I is the (n,n) identity, in a least squares sense for X using a UTV
!  decomposition.
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ainv - the pseudoinverse
!

double precision, allocatable :: u(:,:), v(:,:), t(:,:), z(:,:), z2(:,:)
integer, allocatable          :: ipivs(:)
double precision, allocatable :: rnorms(:)


!
!  Form the U T V decomposition
!

call linalg0_utv(eps,n,m,a,krank,rnorms,ipivs,u,t,v)
call linalg0_utvinverter(n,m,u,t,v,krank,ipivs,ainv)

end subroutine


subroutine linalg0_utvinvert_c(eps,n,m,a,ainv)
implicit double precision (a-h,o-z)
double complex               :: a(:,:), ainv(:,:)
!
!  Solve the equation
!
!    A(n,m) * X(m) = I(n,n),
!
!  where I is the (n,n) identity, in a least squares sense for X using a UTV
!  decomposition.
!
!  Input parameters:
!    eps - the precision for the utv factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ainv - the pseudoinverse
!

double complex, allocatable :: u(:,:), v(:,:), t(:,:), z(:,:), z2(:,:)
integer, allocatable          :: ipivs(:)
double precision, allocatable :: rnorms(:)


!
!  Form the U T V decomposition
!

call linalg0_utv(eps,n,m,a,krank,rnorms,ipivs,u,t,v)
call linalg0_utvinverter(n,m,u,t,v,krank,ipivs,ainv)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Solvers / Inverters based on the SVD
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg0_svdsolver_r1(eps,n,m,u,sigma,vt,x,y)
implicit double precision (a-h,o-z)
double precision       :: u(:,:), vt(:,:), x(:), y(:)
double precision       :: sigma(:)
!
!  Use a pre-existing singular value decomposition of a to solve the
!  system A*X = B.  Singular vectors whose corresponding singular values
!  are below a user-specified thrshold will not be included in the calculations
!
!  Input parameters:
!    eps - the threshold for selecting singular values
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double precision, allocatable :: z(:), z2(:)

!
!  Count the number of singular values above the threshold
!
do krank=1,min(n,m)-1
if (sigma(krank) .lt. eps) exit 
end do

allocate(z(krank))
z  = matmul(transpose(u(:,1:krank)),y)

do i=1,krank
z(i) = z(i) / sigma(i)
end do

x = matmul(transpose(vt(1:krank,:)),z)
end subroutine



subroutine linalg0_svdsolver_c1(eps,n,m,u,sigma,vt,x,y)
implicit double precision (a-h,o-z)
double complex         :: u(:,:), vt(:,:), x(:), y(:)
double precision       :: sigma(:)
!
!  Use a pre-existing singular value decomposition of a to solve the
!  system A*X = B.  Singular vectors whose corresponding singular values
!  are below a user-specified thrshold will not be included in the calculations
!
!  Input parameters:
!    eps - the threshold for selecting singular values
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double complex, allocatable :: z(:), z2(:)

!
!  Count the number of singular values above the threshold
!

do krank=1,min(n,m)-1
if (sigma(krank+1) .lt. eps) exit 
end do

allocate(z(krank))
z  = matmul(conjg(transpose(u(:,1:krank))),y)

do i=1,krank
z(i) = z(i) / sigma(i)
end do

x = matmul(conjg(transpose(vt(1:krank,:))),z)

end subroutine

subroutine linalg0_svdsolver_r2(eps,n,m,u,sigma,vt,x,y)
implicit double precision (a-h,o-z)
double precision       :: u(:,:), vt(:,:), x(:,:), y(:,:)
double precision       :: sigma(:)
!
!  Use a pre-existing singular value decomposition of a to solve the
!  system A*X = B.  Singular vectors whose corresponding singular values
!  are below a user-specified thrshold will not be included in the calculations
!
!  Input parameters:
!    eps - the threshold for singular values
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double precision, allocatable :: z(:,:), z2(:,:)


l = size(y,2)

!
!  Count the number of singular values above the threshold
!
do krank=1,min(n,m)-1
if (sigma(krank+1) .lt. eps) exit 
end do

allocate(z(krank,l))
z  = matmul(transpose(u(:,1:krank)),y)


do i=1,krank
z(i,:) = z(i,:) / sigma(i)
end do

x = matmul(transpose(vt(1:krank,:)),z)

end subroutine


subroutine linalg0_svdsolver_c2(eps,n,m,u,sigma,vt,x,y)
implicit double precision (a-h,o-z)
double complex         :: u(:,:), vt(:,:), x(:,:), y(:,:)
double precision       :: sigma(:)
!
!  Use a pre-existing singular value decomposition of a to solve the
!  system A*X = B.  Singular vectors whose corresponding singular values
!  are below a user-specified thrshold will not be included in the calculations
!
!  Input parameters:
!    eps - the threshold for singular values
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the equation to solve
!
!  Output parameters:
!    x - the obtained solution
!

double complex, allocatable :: z(:,:), z2(:,:)


l = size(y,2)

!
!  Count the number of singular values above the threshold
!
do krank=1,min(n,m)-1
if (sigma(krank+1) .lt. eps) exit 
end do

allocate(z(krank,l))
z  = matmul(conjg(transpose(u(:,1:krank))),y)


do i=1,krank
z(i,:) = z(i,:) / sigma(i)
end do

x = matmul(conjg(transpose(vt(1:krank,:))),z)

end subroutine


subroutine linalg0_svdsolve_r1(ier,eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double precision                  :: a(:,:), x(:), y(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  via an SVD of the input matrix. Singular vectors corresponding
!  to singular values which are below a specified threshold are
!  ignored.
!
!  Input parameters:
!    eps - the theshold for singular values
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    ier - an error return code;
!     ier = 0  indicates successful execution
!     ier = 4  means that the Jacobi iterations failed to
!              converge while trying to compute the SVD
!    x - the solution of the system
!
!
double precision, allocatable   :: u(:,:), vt(:,:)
double precision, allocatable   :: sigma(:)

call linalg0_svd(ier,n,m,a,u,sigma,vt)
if (ier .ne.0 ) return
call linalg0_svdsolver(eps,n,m,u,sigma,vt,x,y)

end subroutine


subroutine linalg0_svdsolve_c1(ier,eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double complex                  :: a(:,:), x(:), y(:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  via an SVD of the input matrix. Singular vectors corresponding
!  to singular values which are below a specified threshold are
!  ignored.
!
!  Input parameters:
!    eps - the theshold for singular values
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    ier - an error return code;
!     ier = 0  indicates successful execution
!     ier = 4  means that the Jacobi iterations failed to
!              converge while trying to compute the SVD
!    x - the solution of the system
!
!
double complex, allocatable   :: u(:,:), vt(:,:)
double precision, allocatable :: sigma(:)

double complex :: ima

call linalg0_svd(ier,n,m,a,u,sigma,vt)
if (ier .ne.0 ) return
call linalg0_svdsolver(eps,n,m,u,sigma,vt,x,y)

end subroutine


subroutine linalg0_svdsolve_r2(ier,eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double precision                  :: a(:,:), x(:,:), y(:,:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  via an SVD of the input matrix. Singular vectors corresponding
!  to singular values which are below a specified threshold are
!  ignored.
!
!  Input parameters:
!    eps - the theshold for singular values
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    ier - an error return code;
!     ier = 0  indicates successful execution
!     ier = 4  means that the Jacobi iterations failed to
!              converge while trying to compute the SVD
!    x - the solution of the system
!
!
double precision, allocatable   :: u(:,:), vt(:,:)
double precision, allocatable   :: sigma(:)

call linalg0_svd(ier,n,m,a,u,sigma,vt)
if (ier .ne.0 ) return
call linalg0_svdsolver(eps,n,m,u,sigma,vt,x,y)

end subroutine


subroutine linalg0_svdsolve_c2(ier,eps,n,m,a,x,y)
implicit double precision (a-h,o-z)
double complex                  :: a(:,:), x(:,:), y(:,:)
!
!  Solve the linear system of equations
!
!    A x = y 
!
!  via an SVD of the input matrix. Singular vectors corresponding
!  to singular values which are below a specified threshold are
!  ignored.
!
!  Input parameters:
!    eps - the theshold for singular values
!    (n,m) - the dimension of the input matrix
!    a -  the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    ier - an error return code;
!     ier = 0  indicates successful execution
!     ier = 4  means that the Jacobi iterations failed to
!              converge while trying to compute the SVD
!    x - the solution of the system
!
!
double complex, allocatable   :: u(:,:), vt(:,:)
double precision, allocatable :: sigma(:)

call linalg0_svd(ier,n,m,a,u,sigma,vt)
if (ier .ne.0 ) return
call linalg0_svdsolver(eps,n,m,u,sigma,vt,x,y)

end subroutine



subroutine linalg0_svdinverter_c(eps,n,m,u,sigma,vt,ainv)
implicit double precision (a-h,o-z)
double complex               :: u(:,:), vt(:,:), ainv(:,:)
double precision             :: sigma(:)
!
!  Form the pseudoinverse of a using a pre-existing SVD of A.  Singular vectors 
!  whose corresponding singular values are below a user-specified thrshold 
!  will not be included in the calculations
!
!
!  Input parameters:
!    eps - the threshold for singular values
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ainv - the pseudoinverse
!

double complex, allocatable   :: z(:,:)


!
!  Count the number of singular values above the threshold
!
krank=0
do i=1,min(n,m)
if (sigma(i) .gt. eps) krank= i
end do

! do krank=1,min(n,m)-1
! if (sigma(krank) .lt. eps) exit 
! end do



allocate(z(krank,n))
z  = conjg(transpose(u(:,1:krank)))
do i=1,krank
z(i,:) = z(i,:)/sigma(i)
end do

ainv = matmul(conjg(transpose(vt(1:krank,:))),z)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Jacobi-based routines for diagonalization and computing SVDs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg0_svd_c(ier,n,m,a,u,sigma,vt)
implicit double precision (a-h,o-z)
double complex                                     :: a(:,:)
double complex, allocatable, intent(out)           :: u(:,:), vt(:,:)
double precision , allocatable, intent(out)        :: sigma(:)
!
!  Form the SVD of an input matrix A; that is, factor A as 
!
!    A(n,m) = U(n,l) D(l,l) V(l,m)^*
!
!  where l = min(n,m), U and V have orthonormal columns and D is a diagonal matrix
!  with positive real-values, arranged in descending order.
!
!  Input parameters:
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0   indicates successful execution
!      ier = 4   means that the Jacobi iterative procedure failed 
!
!    u - the (n,l) matrix u appearing in the factorization
!    sigma - an array of length l whose entries are the singular values
!    vt - the (l,m) matrix giving the transpose of the matrix V 
!
!

double complex, allocatable :: r(:,:), u2(:,:)
double complex, allocatable :: at(:,:), v(:,:), v2(:,:), ut(:,:)

ier   = 0
ifuv = 1
if (n .ge. m) then
call linalg0_jacqr_c(ier,n,m,a,u,r)
if (ier .ne. 0) return
call linalg0_jacsvd_c(ier,ifuv,m,r,u2,sigma,vt)
if (ier .ne. 0) return
u = matmul(u,u2)

else

allocate(at(m,n))
at = transpose(conjg(a))
call linalg0_jacqr_c(ier,m,n,at,v,r)

if (ier .ne. 0) return
call linalg0_jacsvd_c(ier,ifuv,n,r,v2,sigma,ut)
if (ier .ne. 0) return

v = matmul(v,v2)
allocate(vt(n,m), u(n,n))
u = transpose(conjg(ut))
vt = transpose(conjg(v))

endif


end subroutine



subroutine linalg0_jacsvd_c(ier,ifuv,n,a,u,sigma,vt)
implicit double precision (a-h,o-z)
double complex                                     :: a(:,:)
double complex, allocatable, intent(out)           :: u(:,:), vt(:,:)
double precision , allocatable, intent(out)        :: sigma(:)
!
!  Use Jacobi iterations to compute the SVD of a square matrix.  That is,
!  factor the matrix as
!
!    A(n,n) = U(n,n) D(n,n) V^*(n,n)
!
!  where U and V are unitary and D is diagonal with monotonically decreasing
!  entries.
!
!  If the integer flag ifuv is not set to 1, then only the singular values are
!  calculated.
!
!  Input parameters:
!    ifuv - an integer flag indicating whether to compute u and vt or not
!    eps - the threshold for terminating the Jacobi iterations
!    n - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0   indicates successful execution
!      ier = 4   means that the Jacobi iterative procedure failed 
!
!    u - the (n,n) matrix u appearing in the factorization
!    sigma - an array of length n whose entries are the singular values
!    vt - the (n,n) matrix giving the transpose of the matrix V 
!
!

double complex, allocatable   :: b(:,:), work1(:), work2(:)
double complex, allocatable   :: b2(:,:)


double complex                :: a0(2,2), u0(2,2), vt0(2,2), u00(2,2)
double precision              :: rlams0(2)
double precision, allocatable :: dmaxs(:)
integer, allocatable          :: imaxs(:)

! integer, allocatable          :: ipivs(:), imaxs(:)
! double precision, allocatable :: rnorms(:), dmaxs(:)

ier       = 0
maxsweeps = 100000
thresh    = epsilon(0.0d0)*100

!
!  Handle the pathological 1x1 case
!

if (n .eq. 1) then

allocate(sigma(1))

sigma(1) = abs(a(1,1))

if(ifuv .eq. 1) then
allocate(u(1,1),vt(1,1))

vt(1,1) = 1
u(1,1)  = a(1,1) / abs(a(1,1))

endif

return
endif

 
allocate(b(n,n),u(n,n),vt(n,n))
b = a

! allocate(b2(n,n))
! b2 = b

!
!  Do Jacobi-like iterations to diagonalize B
!

allocate(work1(n),work2(n))
allocate(dmaxs(1:n), imaxs(1:n))

if (ifuv .eq. 1) then
u   = 0
vt  = 0
do i=1,n
u(i,i)   = 1
vt(i,i)  = 1
end do
endif

! do iter=1,1
!thresh =sqrt(epsilon(0.0d0))
!thresh = thresh**2

!
!  Find the maximum off-diagonal element in each row
!

do i=1,n
dmaxs(i) = 0
imaxs(i) = -1

do j=1,n
if (i .eq. j) cycle

dd = abs(b(i,j))
if (dd .gt. dmaxs(i)) then
dmaxs(i) = dd
imaxs(i) = j
endif

end do
end do



do isweep = 1,maxsweeps

!
!  Find the largest nondiagonal entry in the matrix
!

dmax =  0
irow = -1
icol = -1

do i=1,n
dd = dmaxs(i)
if (dd .gt. dmax) then
dmax = dd
irow = i
icol = imaxs(i)
endif
end do


!
!  Terminate the procedure is the remaining entries are
!  below the threshold 
!
if (dmax .lt. thresh) exit


!
!  Find the SVD of the appropriate 2x2 block
!

a0(1,1) = b(irow,irow)
a0(1,2) = b(irow,icol)
a0(2,1) = b(icol,irow)
a0(2,2) = b(icol,icol)
!call linalg0_svd2x2_c(a0,u0,rlams0,vt0)

call vlad_oneuv_c(a0,rlams0,u0,vt0)
u0 = transpose(conjg(u0))
vt0 = transpose(conjg(vt0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check the accuracy of the 2x2 SVD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! u00(:,1) = u0(:,1)*rlams0(1)
! u00(:,2) = u0(:,2)*rlams0(2)
! ! print *,rlams0(1),rlams0(2)
! print *,norm2(abs(matmul(u0,transpose(conjg(u0)))-eye(2)))
! print *,norm2(abs(matmul(transpose(conjg(u0)),u0)-eye(2)))
! print *,norm2(abs(matmul(vt0,transpose(conjg(vt0)))-eye(2)))
! print *,norm2(abs(matmul(transpose(conjg(vt0)),vt0)-eye(2)))
! print *,norm2(abs(a0-matmul(u00,vt0)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if (ifuv .eq. 1) then
! Apply VT to the left side of VT
work1       = vt0(1,1)*vt(irow,:) + vt0(1,2) * vt(icol,:)
work2       = vt0(2,1)*vt(irow,:) + vt0(2,2) * vt(icol,:)
vt(irow,:) = work1
vt(icol,:) = work2

! Apply u0 to the right side of u
work1       = u0(1,1)*u(:,irow) + u0(2,1) * u(:,icol)
work2       = u0(1,2)*u(:,irow) + u0(2,2) * u(:,icol)
u(:,irow)  = work1
u(:,icol)  = work2
endif


! Apply u0^* to the left side of b and v0 to the right side of b

u0  = conjg(transpose(u0))
vt0 = conjg(transpose(vt0))

work1      = u0(1,1)*b(irow,:) + u0(1,2) * b(icol,:)
work2      = u0(2,1)*b(irow,:) + u0(2,2) * b(icol,:)
b(irow,:)  = work1
b(icol,:)  = work2

work1      = vt0(1,1)*b(:,irow) + vt0(2,1) * b(:,icol)
work2      = vt0(1,2)*b(:,irow) + vt0(2,2) * b(:,icol)
b(:,irow)  = work1
b(:,icol)  = work2





!print *,irow,icol,dmax,norm2(abs(matmul(matmul(u,b),vt)-b2))



!
!  Update the dmaxs array
!

dmaxs(irow) = 0
imaxs(irow) = -1

dmaxs(icol) = 0
imaxs(icol) = -1


do j=1,n

if (j .ne. irow) then
dd = abs(b(irow,j))
if (dd .gt. dmaxs(irow)) then
dmaxs(irow) = dd
imaxs(irow) = j
endif
endif
if (j .ne. icol) then
dd = abs(b(icol,j))
if (dd .gt. dmaxs(icol)) then
dmaxs(icol) = dd
imaxs(icol) = j
endif
endif

end do

do i=1,n
if (i .ne. irow ) then
dd = abs(b(i,irow))
if (dd .gt. dmaxs(i)) then
dmaxs(i)  = dd
imaxs(i)  = irow
endif
endif


if (i .ne. icol) then
dd = abs(b(i,icol))
if (dd .ge. dmaxs(i)) then
dmaxs(i)  = dd
imaxs(i)  = icol
endif
endif
end do

end do
! end do


if (isweep .gt. maxsweeps) then
ier = 4
return
endif

!
!  Extract the singular values from b
!

allocate(sigma(n))
do i=1,n
sigma(i)=b(i,i)

if (sigma(i) .lt. 0) then
sigma(i) = -sigma(i)
vt(i,:)  = - vt(i,:)
endif

end do


!
!  Reorder the matrices so as to arrange the singular values in descending order
!  using a simple insertion sort.
!


do i=2,n
val       = sigma(i)

if (ifuv .eq. 1) then
work1     = u(:,i)
work2     = vt(i,:)
endif

j         = i-1
do while (j .ge. 1 .AND. sigma(j) .lt. val) 
sigma(j+1) = sigma(j)
if (ifuv .eq. 1) then
u(:,j+1)   = u(:,j)
vt(j+1,:)  = vt(j,:)
endif
j          = j-1
end do
sigma(j+1) = val
if (ifuv .eq. 1) then
u(:,j+1)   = work1
vt(j+1,:)  = work2
endif

end do


end subroutine




subroutine linalg0_jacqr_c(ier,n,m,a,q,r)
implicit double precision (a-h,o-z)
double complex                                     :: a(:,:)
double complex, allocatable, intent(out)           :: q(:,:), r(:,:)
!
!  Use Jacobi iterations to compute the QR decomposition of a rectangular
!  matrix.  That is, factor the input matrix as
!
!    A(n,m) = Q(n,n) R(n,m)
!
!  where l = min(n,m), R is upper trapezoidal and Q is unitary.
!
!  Input parameters:
!    eps - the threshold for terminating the Jacobi iterations
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0   indicates successful execution
!      ier = 4   means that the Jacobi iterative procedure failed 
!
!  q - the matrix q in the factorization
!  r - the matrix r in the factorization
!
!

double complex, allocatable   :: b(:,:), work1(:), work2(:), work3(:), work4(:)
double complex, allocatable   :: b2(:,:), q0(:,:), r0(:,:)
double complex                :: a0(2,2), u0(2,2)


allocate(r0(n,m),q0(n,n))
allocate(work1(m),work2(m))
allocate(work3(n),work4(n))


r0 = a

q0 = 0
do i=1,n
q0(i,i) = 1
end do


!
!  Eliminate the (i,j) entry from r0 and update q0
!
do j=1,m
do i=n,j+1,-1
call linalg0_givens_c(r0(i-1,j),r0(i,j),u0)
work1      = u0(1,1)*r0(i-1,:) + u0(1,2) * r0(i,:)
work2      = u0(2,1)*r0(i-1,:) + u0(2,2) * r0(i,:)
r0(i-1,:)  = work1
r0(i,:)    = work2

u0 = transpose(conjg(u0))
work3      = u0(1,1)*q0(:,i-1) + u0(2,1) * q0(:,i)
work4      = u0(1,2)*q0(:,i-1) + u0(2,2) * q0(:,i)
q0(:,i-1)  = work3
q0(:,i)    = work4
end do
end do


! print *,norm2(abs(a-matmul(q0,r0)))
! print *,norm2(abs(matmul(transpose(conjg(q0)),q0)-eye(n)))


!
!  Neuter the matrices
!

l = min(n,m)
allocate(q(n,l),r(l,m))
q = q0(:,1:l)
r = r0(1:l,:)


! print *,norm2(abs(a-matmul(q,r)))
! print *,norm2(abs(matmul(transpose(conjg(q)),q)-eye(l)))

end subroutine


subroutine linalg0_jaceig_c(ier,eps,n,a,u,sigma)
implicit double precision (a-h,o-z)
double complex                 :: a(:,:), u(:,:)
double precision               :: sigma(:)
!
!  Use the complex Jacobi method to factor a Hermitian matrix A as 
!
!    A(n,n) = U(n,n) D (n,n) U(n,n)^*
!
!  with U unitary and D diagonal using Jacobi iterations.  The columns of
!  U are ordered so that the entries of D are in descendning order.
!
!  Input parameters:
!    eps - the threshold for terminating the procedure
!    n - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ier - an error return code;
!      ier  = 0  indicates successful execution
!      ier  = 4  means that maximum number of sweeps was exceeded
!
!    u - the (n,n) orthogonal matrix
!    sigma - an array giving the n entries of the diagonal matrix
!

double precision, allocatable :: dmaxs(:)
integer, allocatable          :: imaxs(:)
double complex, allocatable   :: v1(:), v2(:), b(:,:)
double complex                :: a0(2,2), u0(2,2)


maxsweeps = 100000
dnorm     = maxval(abs(a))
thresh    = epsilon(0.0d0)*1000

allocate(v1(n), v2(n))

ier  =  0

allocate(b(n,n))
b = a(1:n,1:n)

!
!  Initialize u
!
u = 0
do i=1,n
u(i,i) = 1
end do

!
!  Find the maximum off-diagonal element in each row, we are ignoring symmetry
!  considerations here for ease of programming
!

allocate(dmaxs(1:n), imaxs(1:n))
do i=1,n
dmaxs(i) =0
imaxs(i) = -1

do j=1,n
if (i .eq. j) cycle

dd = abs(b(i,j))
if (dd .gt. dmaxs(i)) then
dmaxs(i) = dd
imaxs(i) = j
endif

end do
end do


do isweep = 1,maxsweeps

!
!  Find the largest entry in the matrix
!

dmax =  0
irow = -1
icol = -1

do i=1,n
dd = dmaxs(i)
if (dd .gt. dmax) then
dmax = dd
irow = i
icol = imaxs(i)
endif
end do


!print *,isweep,irow,icol,dmax

!
!  Terminate the procedure is the remaining entries are
!  below the threshold 
!
if (dmax .lt. thresh) exit

!
!  Find the Jacobi rotation to apply
!

a0(1,1) = b(irow,irow)
a0(1,2) = b(irow,icol)
a0(2,1) = b(icol,irow)
a0(2,2) = b(icol,icol)
call linalg0_jacobi2x2_c(a0,u0)

!
!  Apply U0^* to the left side of V
!
v1        = u0(1,1)*u(:,irow) + u0(2,1) * u(:,icol)
v2        = u0(1,2)*u(:,irow) + u0(2,2) * u(:,icol)
u(:,irow) = v1
u(:,icol) = v2

!
!  Apply the transform A = U_0^* A U_0
!

v1        = u0(1,1)*b(:,irow) + u0(2,1) * b(:,icol)
v2        = u0(1,2)*b(:,irow) + u0(2,2) * b(:,icol)
b(:,irow) = v1
b(:,icol) = v2

u0 = conjg(transpose(u0))
v1        = u0(1,1)*b(irow,:) + u0(1,2) * b(icol,:)
v2        = u0(2,1)*b(irow,:) + u0(2,2) * b(icol,:)
b(irow,:) = v1
b(icol,:) = v2



!
!  Update the dmaxs array
!

dmaxs(irow) = 0
imaxs(irow) = -1

dmaxs(icol) = 0
imaxs(icol) = -1


do j=1,n

if (j .ne. irow) then
dd = abs(b(irow,j))
if (dd .gt. dmaxs(irow)) then
dmaxs(irow) = dd
imaxs(irow) = j
endif
endif
if (j .ne. icol) then
dd = abs(b(icol,j))
if (dd .gt. dmaxs(icol)) then
dmaxs(icol) = dd
imaxs(icol) = j
endif
endif

end do


do i=1,n

if (i .ne. irow ) then
dd = abs(b(i,irow))
if (dd .gt. dmaxs(i)) then
dmaxs(i)  = dd
imaxs(i)  = irow
endif
endif


if (i .ne. icol) then
dd = abs(b(i,icol))
if (dd .ge. dmaxs(i)) then
dmaxs(i)  = dd
imaxs(i)  = icol
endif
endif
end do

end do

!
!  Copy out the diagonal entries and make sure they are positive
!

do i=1,n
sigma(i) = b(i,i)

! The sign might be negative if the value is close to 0
if (sigma(i) .lt. 0) then
sigma(i) = -sigma(i)
u(:,i)   = -u(:,i)
endif

end do

!
!  Order the entries using a simple insertion sort
!


do i=2,n
val       = sigma(i)
v1        = u(:,i)
j         = i-1
do while (j .ge. 1 .AND. sigma(j) .lt. val) 
sigma(j+1) = sigma(j)
u(:,j+1)   = u(:,j)
j          = j-1
end do
sigma(j+1) = val
u(:,j+1)   = v1
end do


end subroutine



subroutine linalg0_svd2x2_c(a0,u0,rlams,vt0)
implicit double precision (a-h,o-z)
double complex      :: a0(2,2), u0(2,2), vt0(2,2)
double precision    :: rlams(2)
!
!  Factor a 2x2 complex-valued matrix A as
!
!     U0 D V0^*
!
!  with U0, V0 unitary and D diagonal.
!
!
double complex    :: a12
double complex    :: b0(2,2), j0(2,2)

end subroutine


subroutine linalg0_jacobi2x2_c(a0,u0)
implicit double precision (a-h,o-z)
double complex      :: a0(2,2), u0(2,2)
!
!  Given a (2,2) Hermitian matrix A0, find a (2,2) unitary matrix U0 such
!  that 
!
!    U0^* A0 U0  
!
!  is diagonal with real-valued entries.  
!
double complex    :: a12
eps = 0.0d0

a11 = real(a0(1,1))
a12 = a0(1,2)
a22 = real(a0(2,2))

if (abs(a12) .le. eps) then
u0(1,1)  = 1.0d0
u0(2,2)  = 1.0d0
u0(1,2)  = 0.0d0
u0(2,1)  = 0.0d0
return
endif

t = 2*abs(a12) * sign(1.0d0,a11-a22) / (abs(a11-a22) + sqrt( abs(a11-a22)**2 + 4*abs(a12)**2) )
c = 1/sqrt(1+t**2)
s = t*c

u0(1,1) = c
u0(1,2) = -s * a12/abs(a12)
u0(2,1) =  s  * abs(a12)/a12
u0(2,2) = c

end subroutine



subroutine linalg0_givens_c(a,b,u0)
implicit double precision (a-h,o-z)
double complex      :: u0(2,2), a, b, s, c
!
!  Return a (2,2) unitary matrix U0 such that
!
!      [ a ]      [ * ]
!   u0 [ b ]   =  [ 0 ]
!

if (b .eq. 0) then
c = 1
s = 0
elseif (a .eq. 0) then
c = 0
s = b/abs(b)
else
c = abs(a)/sqrt(abs(a)**2+abs(b)**2)
s = a/abs(a) *conjg(b) /sqrt(abs(a)**2+abs(b)**2)
endif 

u0(1,1) = c
u0(1,2) = s
u0(2,1) = -conjg(s)
u0(2,2) = c

end subroutine




subroutine linalg0_svd_r(ier,n,m,a,u,sigma,vt)
implicit double precision (a-h,o-z)
double precision                                     :: a(:,:)
double precision, allocatable, intent(out)           :: u(:,:), vt(:,:)
double precision , allocatable, intent(out)          :: sigma(:)
!
!  Form the SVD of an input matrix A; that is, factor A as 
!
!    A(n,m) = U(n,l) D(l,l) V(l,m)^t
!
!  where l = min(n,m), U and V have orthonormal columns and D is a diagonal matrix
!  with positive real-values, arranged in descending order.
!
!  Input parameters:
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0   indicates successful execution
!      ier = 4   means that the Jacobi iterative procedure failed 
!
!    u - the (n,l) matrix u appearing in the factorization
!    sigma - an array of length l whose entries are the singular values
!    vt - the (l,m) matrix giving the transpose of the matrix V 
!
!

double precision, allocatable :: r(:,:), u2(:,:)
double precision, allocatable :: at(:,:), v(:,:), v2(:,:), ut(:,:)


ier   = 0
ifuv = 1
if (n .ge. m) then

call linalg0_jacqr_r(ier,n,m,a,u,r)
if (ier .ne. 0) return
call linalg0_jacsvd_r(ier,ifuv,m,r,u2,sigma,vt)
if (ier .ne. 0) return
u = matmul(u,u2)

else

allocate(at(m,n))
at = transpose(a)
call linalg0_jacqr_r(ier,m,n,at,v,r)
if (ier .ne. 0) return
call linalg0_jacsvd_r(ier,ifuv,n,r,v2,sigma,ut)
if (ier .ne. 0) return

v = matmul(v,v2)
allocate(vt(n,m), u(n,n))
u = transpose(ut)
vt = transpose(v)

endif


end subroutine



subroutine linalg0_jacsvd_r(ier,ifuv,n,a,u,sigma,vt)
implicit double precision (a-h,o-z)
double precision                                     :: a(:,:)
double precision, allocatable, intent(out)           :: u(:,:), vt(:,:)
double precision , allocatable, intent(out)          :: sigma(:)
!
!  Use Jacobi iterations to compute the SVD of a square matrix.  That is,
!  factor the matrix as
!
!    A(n,n) = U(n,n) D(n,n) V^t(n,n)
!
!  where U and V are orthgonal and D is diagonal with monotonically decreasing
!  entries.
!
!  If the integer flag ifuv is not set to 1, then only the singular values are
!  calculated.
!
!  Input parameters:
!    ifuv - an integer flag indicating whether to compute u and vt or not
!    eps - the threshold for terminating the Jacobi iterations
!    n - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0   indicates successful execution
!      ier = 4   means that the Jacobi iterative procedure failed 
!
!    u - the (n,n) matrix u appearing in the factorization
!    sigma - an array of length n whose entries are the singular values
!    vt - the (n,n) matrix giving the transpose of the matrix V 
!
!

double precision, allocatable   :: b(:,:), work1(:), work2(:)
double precision, allocatable   :: b2(:,:)


double precision              :: a0(2,2), u0(2,2), vt0(2,2), u00(2,2)
double precision              :: rlams0(2)
double precision, allocatable :: dmaxs(:)
integer, allocatable          :: imaxs(:)

! integer, allocatable          :: ipivs(:), imaxs(:)
! double precision, allocatable :: rnorms(:), dmaxs(:)

ier       = 0
maxsweeps = 100000
thresh    = epsilon(0.0d0)*1000

!
!  Handle the pathological 1x1 case
!

if (n .eq. 1) then

allocate(sigma(1))

sigma(1) = abs(a(1,1))

if(ifuv .eq. 1) then
allocate(u(1,1),vt(1,1))

vt(1,1) = 1
u(1,1)  = a(1,1) / abs(a(1,1))

endif

return
endif

 
allocate(b(n,n),u(n,n),vt(n,n))
b = a

! allocate(b2(n,n))
! b2 = b

!
!  Do Jacobi-like iterations to diagonalize B
!

allocate(work1(n),work2(n))
allocate(dmaxs(1:n), imaxs(1:n))

if (ifuv .eq. 1) then
u   = 0
vt  = 0
do i=1,n
u(i,i)   = 1
vt(i,i)  = 1
end do
endif

! do iter=1,1
! thresh =sqrt(epsilon(0.0d0))
!thresh = thresh**2

!
!  Find the maximum off-diagonal element in each row
!

do i=1,n
dmaxs(i) = 0
imaxs(i) = -1

do j=1,n
if (i .eq. j) cycle

dd = abs(b(i,j))
if (dd .gt. dmaxs(i)) then
dmaxs(i) = dd
imaxs(i) = j
endif

end do
end do



do isweep = 1,maxsweeps

!
!  Find the largest nondiagonal entry in the matrix
!

dmax =  0
irow = -1
icol = -1

do i=1,n
dd = dmaxs(i)
if (dd .gt. dmax) then
dmax = dd
irow = i
icol = imaxs(i)
endif
end do


!
!  Terminate the procedure is the remaining entries are
!  below the threshold 
!
if (dmax .lt. thresh) exit



!
!  Find the SVD of the appropriate 2x2 block
!

a0(1,1) = b(irow,irow)
a0(1,2) = b(irow,icol)
a0(2,1) = b(icol,irow)
a0(2,2) = b(icol,icol)
!call linalg0_svd2x2_c(a0,u0,rlams0,vt0)

call vlad_oneuv_r(a0,rlams0,u0,vt0)
!u0 = transpose(u0)
!vt0 = transpose(vt0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check the accuracy of the 2x2 SVD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! u00(:,1) = u0(:,1)*rlams0(1)
! u00(:,2) = u0(:,2)*rlams0(2)
! print *,norm2(abs(matmul(u0,transpose(u0))-eye(2)))
! print *,norm2(abs(matmul(transpose(u0),u0)-eye(2)))
! print *,norm2(abs(matmul(vt0,transpose(vt0))-eye(2)))
! print *,norm2(abs(matmul(transpose(vt0),vt0)-eye(2)))
! print *,norm2(abs(a0-matmul(u00,vt0)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if (ifuv .eq. 1) then
! Apply VT to the left side of VT
work1       = vt0(1,1)*vt(irow,:) + vt0(1,2) * vt(icol,:)
work2       = vt0(2,1)*vt(irow,:) + vt0(2,2) * vt(icol,:)
vt(irow,:) = work1
vt(icol,:) = work2

! Apply u0 to the right side of u
work1       = u0(1,1)*u(:,irow) + u0(2,1) * u(:,icol)
work2       = u0(1,2)*u(:,irow) + u0(2,2) * u(:,icol)
u(:,irow)  = work1
u(:,icol)  = work2
endif


! Apply u0^* to the left side of b and v0 to the right side of b

u0  = transpose(u0)
vt0 = transpose(vt0)

work1      = u0(1,1)*b(irow,:) + u0(1,2) * b(icol,:)
work2      = u0(2,1)*b(irow,:) + u0(2,2) * b(icol,:)
b(irow,:)  = work1
b(icol,:)  = work2

work1      = vt0(1,1)*b(:,irow) + vt0(2,1) * b(:,icol)
work2      = vt0(1,2)*b(:,irow) + vt0(2,2) * b(:,icol)
b(:,irow)  = work1
b(:,icol)  = work2





!print *,irow,icol,dmax,thresh
!norm2(abs(matmul(matmul(u,b),vt)-b2))



!
!  Update the dmaxs array
!

dmaxs(irow) = 0
imaxs(irow) = -1

dmaxs(icol) = 0
imaxs(icol) = -1


do j=1,n

if (j .ne. irow) then
dd = abs(b(irow,j))
if (dd .gt. dmaxs(irow)) then
dmaxs(irow) = dd
imaxs(irow) = j
endif
endif
if (j .ne. icol) then
dd = abs(b(icol,j))
if (dd .gt. dmaxs(icol)) then
dmaxs(icol) = dd
imaxs(icol) = j
endif
endif

end do

do i=1,n
if (i .ne. irow ) then
dd = abs(b(i,irow))
if (dd .gt. dmaxs(i)) then
dmaxs(i)  = dd
imaxs(i)  = irow
endif
endif


if (i .ne. icol) then
dd = abs(b(i,icol))
if (dd .ge. dmaxs(i)) then
dmaxs(i)  = dd
imaxs(i)  = icol
endif
endif
end do

end do
! end do


if (isweep .gt. maxsweeps) then
ier = 4
return
endif

!
!  Extract the singular values from b
!

allocate(sigma(n))
do i=1,n
sigma(i)=b(i,i)

if (sigma(i) .lt. 0) then
sigma(i) = -sigma(i)
vt(i,:)  = - vt(i,:)
endif

end do


!
!  Reorder the matrices so as to arrange the singular values in descending order
!  using a simple insertion sort.
!


do i=2,n
val       = sigma(i)

if (ifuv .eq. 1) then
work1     = u(:,i)
work2     = vt(i,:)
endif

j         = i-1
do while (j .ge. 1 .AND. sigma(j) .lt. val) 
sigma(j+1) = sigma(j)
if (ifuv .eq. 1) then
u(:,j+1)   = u(:,j)
vt(j+1,:)  = vt(j,:)
endif
j          = j-1
end do
sigma(j+1) = val
if (ifuv .eq. 1) then
u(:,j+1)   = work1
vt(j+1,:)  = work2
endif

end do


end subroutine


subroutine linalg0_jacqr_r(ier,n,m,a,q,r)
implicit double precision (a-h,o-z)
double precision                                   :: a(:,:)
double precision, allocatable, intent(out)         :: q(:,:), r(:,:)
!
!  Use Jacobi iterations to compute the QR decomposition of a rectangular
!  matrix.  That is, factor the input matrix as
!
!    A(n,m) = Q(n,n) R(n,m)
!
!  where l = min(n,m), R is upper trapezoidal and Q is orthogonal.
!
!  Input parameters:
!    eps - the threshold for terminating the Jacobi iterations
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0   indicates successful execution
!      ier = 4   means that the Jacobi iterative procedure failed 
!
!  q - the matrix q in the factorization
!  r - the matrix r in the factorization
!
!

double precision, allocatable   :: b(:,:), work1(:), work2(:), work3(:), work4(:)
double precision, allocatable   :: b2(:,:), q0(:,:), r0(:,:)
double precision                :: a0(2,2), u0(2,2)


allocate(r0(n,m),q0(n,n))
allocate(work1(m),work2(m))
allocate(work3(n),work4(n))


r0 = a
q0 = 0
do i=1,n
q0(i,i) = 1
end do


!
!  Eliminate the (i,j) entry from r0 and update q0
!
do j=1,m
do i=n,j+1,-1
call linalg0_givens_r(r0(i-1,j),r0(i,j),u0)
work1      = u0(1,1)*r0(i-1,:) + u0(1,2) * r0(i,:)
work2      = u0(2,1)*r0(i-1,:) + u0(2,2) * r0(i,:)
r0(i-1,:)  = work1
r0(i,:)    = work2

u0 = transpose(u0)
work3      = u0(1,1)*q0(:,i-1) + u0(2,1) * q0(:,i)
work4      = u0(1,2)*q0(:,i-1) + u0(2,2) * q0(:,i)
q0(:,i-1)  = work3
q0(:,i)    = work4
end do
end do



!
!  Neuter the matrices
!

l = min(n,m)
allocate(q(n,l),r(l,m))
q = q0(:,1:l)
r = r0(1:l,:)


end subroutine



subroutine linalg0_jaceig_r(ier,eps,n,a,u,sigma)
implicit double precision (a-h,o-z)
double precision               :: a(:,:), u(:,:)
double precision               :: sigma(:)
!
!  Use the real Jacobi method to factor a symmetric matrix A as 
!
!    A(n,n) = U(n,n) D (n,n) U(n,n)^t
!
!  with U unitary and D diagonal with positive real entries using Jacobi iterations. 
!  The columns of U are arranged so that the entries of D are in descending order.
!
!  Input parameters:
!    eps - the threshold for terminating the procedure
!    n - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    ier - an error return code;
!      ier  = 0  indicates successful execution
!      ier  = 4  means that maximum number of sweeps was exceeded
!
!    u - the (n,n) orthogonal matrix
!    sigma - an array giving the n entries of the diagonal matrix
!

double precision, allocatable :: dmaxs(:), b(:,:)
integer, allocatable          :: imaxs(:)
double precision, allocatable :: v1(:), v2(:)
double precision              :: a0(2,2), u0(2,2)
double precision              :: rlams(2)


maxsweeps = 100000
dnorm     = maxval(abs(a))
thresh    = epsilon(0.0d0)*1000

allocate(b(n,n))
b = a

allocate(v1(n), v2(n))

ier  =  0

!
!  Initialize u
!
u = 0
do i=1,n
u(i,i) = 1
end do

!
!  Find the maximum off-diagonal element in each row, we are ignoring symmetry
!  considerations here for ease of programming
!

allocate(dmaxs(1:n), imaxs(1:n))
do i=1,n
dmaxs(i) =0
imaxs(i) = -1

do j=1,n
if (i .eq. j) cycle

dd = abs(b(i,j))
if (dd .gt. dmaxs(i)) then
dmaxs(i) = dd
imaxs(i) = j
endif

end do
end do


do isweep = 1,maxsweeps

!
!  Find the largest entry in the matrix
!

dmax =  0
irow = -1
icol = -1

do i=1,n
dd = dmaxs(i)
if (dd .gt. dmax) then
dmax = dd
irow = i
icol = imaxs(i)
endif
end do


!print *,isweep,irow,icol,dmax

!
!  Terminate the procedure is the remaining entries are
!  below the threshold 
!
if (dmax .lt. thresh) exit

!
!  Find the Jacobi rotation to apply
!

a0(1,1) = b(irow,irow)
a0(1,2) = b(irow,icol)
a0(2,1) = b(icol,irow)
a0(2,2) = b(icol,icol)
call linalg0_jacobi2x2_r(a0,u0)

!
!  Apply U0^* to the left side of V
!
v1        = u0(1,1)*u(:,irow) + u0(2,1) * u(:,icol)
v2        = u0(1,2)*u(:,irow) + u0(2,2) * u(:,icol)
u(:,irow) = v1
u(:,icol) = v2

!
!  Apply the transform A = U_0^* A U_0
!

v1        = u0(1,1)*b(:,irow) + u0(2,1) * b(:,icol)
v2        = u0(1,2)*b(:,irow) + u0(2,2) * b(:,icol)
b(:,irow) = v1
b(:,icol) = v2

u0 = transpose(u0)
v1        = u0(1,1)*b(irow,:) + u0(1,2) * b(icol,:)
v2        = u0(2,1)*b(irow,:) + u0(2,2) * b(icol,:)
b(irow,:) = v1
b(icol,:) = v2


!
!  Update the dmaxs array
!

dmaxs(irow) = 0
imaxs(irow) = -1

dmaxs(icol) = 0
imaxs(icol) = -1


do j=1,n

if (j .ne. irow) then
dd = abs(b(irow,j))
if (dd .gt. dmaxs(irow)) then
dmaxs(irow) = dd
imaxs(irow) = j
endif
endif
if (j .ne. icol) then
dd = abs(b(icol,j))
if (dd .gt. dmaxs(icol)) then
dmaxs(icol) = dd
imaxs(icol) = j
endif
endif

end do


do i=1,n

if (i .ne. irow ) then
dd = abs(b(i,irow))
if (dd .gt. dmaxs(i)) then
dmaxs(i)  = dd
imaxs(i)  = irow
endif
endif


if (i .ne. icol) then
dd = abs(b(i,icol))
if (dd .ge. dmaxs(i)) then
dmaxs(i)  = dd
imaxs(i)  = icol
endif
endif
end do

end do

!
!  Copy out the diagonal entries and make sure they are positive
!

do i=1,n
sigma(i) = b(i,i)

! The sign might be negative if the value is close to 0
if (sigma(i) .lt. 0) then
sigma(i) = -sigma(i)
u(:,i)   = -u(:,i)
endif

end do


!
!  Order the entries using a single insertion sort
!

do i=2,n
val       = sigma(i)
v1        = u(:,i)
j         = i-1
do while (j .ge. 1 .AND. sigma(j) .lt. val) 
sigma(j+1) = sigma(j)
u(:,j+1)   = u(:,j)
j          = j-1
end do
sigma(j+1) = val
u(:,j+1)   = v1
end do

end subroutine


subroutine linalg0_givens_r(a,b,u0)
implicit double precision (a-h,o-z)
double precision      :: u0(2,2)
!
!  Return a (2,2) orthogonal matrix U0 such that
!
!      [ a ]      [ * ]
!   u0 [ b ]   =  [ 0 ]
!

if (b .eq. 0) then
c = 1
s = 0
elseif (a .eq. 0) then
c = 0
s = b/abs(b)
else
c = abs(a)/sqrt(abs(a)**2+abs(b)**2)
s = a/abs(a) *b /sqrt(abs(a)**2+abs(b)**2)
endif 

u0(1,1) = c
u0(1,2) = s
u0(2,1) = -s
u0(2,2) = c

end subroutine


subroutine linalg0_svd2x2_r(a0,u0,rlams,vt0)
implicit double precision (a-h,o-z)
double precision      :: a0(2,2), u0(2,2), vt0(2,2)
double precision      :: rlams(2)
!
!  Factor a 2x2 real-valued matrix A as
!
!     U0 D V0^T
!
!  with U0, V0 orthogonal and D diagonal.
!
!
double precision    :: b0(2,2), j0(2,2)

end subroutine


subroutine linalg0_jacobi2x2_r(a0,u0)
implicit double precision (a-h,o-z)
double precision    :: a0(2,2), u0(2,2)
!
!  Given a (2,2) symmetric matrix A0, find a (2,2) orthogonal matrix U0 such
!  that 
!
!    U0^t A0 U0  
!
!  is diagonal with real-valued entries.  
!
eps = epsilon(0.0d0)**2

a11 = real(a0(1,1))
a12 = a0(1,2)
a22 = real(a0(2,2))

if (abs(a12) .le. eps) then
u0(1,1)  = 1
u0(1,2)  = 0
u0(2,1)  = 0
u0(2,2)  = 1
return
endif


beta  = (a11 - a22) / (2*a12)
t     = sign(1.0d0,beta) / ( abs(beta)**2 + sqrt(beta**2+1) )
c     = 1/sqrt(1+t**2)
s     = t*c

u0(1,1) = c
u0(1,2) = s
u0(2,1) = -s
u0(2,2) = c


end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Low-level pivoting /QR code 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg0_qr3_r(n,m,a,rnorms,ipivs,q,r)
implicit double precision (a-h,o-z)
double precision                             :: a(:,:)
integer, allocatable, intent(out)            :: ipivs(:)
double precision, allocatable, intent(out)   :: rnorms(:)
double precision, allocatable, intent(out)   :: q(:,:)
double precision, allocatable, intent(out)   :: r(:,:)
!
!  Form a QR decomposition of an input matrix A. That is, factor
!  A as
! 
!     A(n,m) = Q(n,l) * R(l,m)
!
!  where l = min(n,m), Q has orthonormal columns and R(:,ipivs) 
!  is upper trapezoidal.
!
!  This routine does NOT destroy the input matrix.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    rnorms - an array of length m specifying the normalizing factors
!      obtained from the GS procedure
!    ipivs - an integer array of length m specifying the list of pivots
!      from the GS procedure
!    q - the matrix q in the factorization
!    r - the matrix r in the factorization
!
!

double precision, allocatable   :: b(:,:)

l = min(n,m)

! Make a copy of the input array.
allocate(b(n,m))
b = a

allocate(ipivs(m),rnorms(m))
call linalg0_gspiv_r(0.0d0,n,m,b,krank,rnorms,ipivs)

allocate(q(n,l))
q   = b(:,1:l)
r   = matmul(transpose(q),a)

end subroutine


subroutine linalg0_qr3_c(n,m,a,rnorms,ipivs,q,r)
implicit double precision (a-h,o-z)
double complex                               :: a(:,:)
integer, allocatable, intent(out)            :: ipivs(:)
double precision, allocatable, intent(out)   :: rnorms(:)
double complex, allocatable, intent(out)     :: q(:,:)
double complex, allocatable, intent(out)     :: r(:,:)
!
!  Form a QR decomposition of an input matrix A. That is, factor
!  A as
! 
!     A(n,m) = Q(n,l) * R(l,m)
!
!  where l = min(n,m), Q has orthonormal columns and R(:,ipivs) 
!  is upper trapezoidal.
!
!  This routine does NOT destroy the input matrix.
!
!  Input parameters:
!    (n,m) - the dimension of the input matrix
!    a - the input matrix
!
!  Output parameters:
!    rnorms - an array of length m specifying the normalizing factors
!      obtained from the GS procedure
!    ipivs - an integer array of length m specifying the list of pivots
!      from the GS procedure
!    q - the matrix q in the factorization
!    r - the matrix r in the factorization
!
!

double complex, allocatable     :: b(:,:)

l = min(n,m)

! Make a copy of the input array.
allocate(b(n,m))
b = a

allocate(ipivs(m),rnorms(m))
call linalg0_gspiv_c(0.0d0,n,m,b,krank,rnorms,ipivs)

allocate(q(n,l))
q   = b(:,1:l)
r   = matmul(conjg(transpose(q)),a)

end subroutine



subroutine linalg0_gspiv_r(eps,n,m,b,krank,rnorms,ipivs)
implicit double precision (a-h,o-z)
double precision :: rnorms(:)
integer          :: ipivs(:)
double precision :: b(:,:)
!
!  Perform Gram-Schmidt with pivoting and reorthogonalization on the input
!  matrix.
!
double precision :: work(n)

done = 1
dtot = 0
do i=1,m
ipivs(i) = i
d          = dot_product(b(:,i),b(:,i))
dtot       = dtot+d
rnorms(i)  = sqrt(d)
end do


thresh=dtot*eps**2 
thresh=sqrt(thresh)

do i=1,min(n,m)

! choose the next column
ipivot = i
rn     = rnorms(i)
do j=i+1,m
if (rnorms(j) .gt. rn) then
rn     = rnorms(j)
ipivot = j
endif
end do

! swap the chosen column with the ith column
work        = b(:,i)
b(:,i)      = b(:,ipivot)
b(:,ipivot) = work

ii              = ipivs(i)
ipivs(i)      = ipivs(ipivot)
ipivs(ipivot) = ii

d               = rnorms(i)
rnorms(i)       = rnorms(ipivot)
rnorms(ipivot)  = d

! orthogonalize the chosen column to all previous columns
do j=1,i-1
cd = dot_product(b(:,j),b(:,i))
b(:,i) = b(:,i) - cd*b(:,j)
end do

! recompute its norm
d  = dot_product(b(:,i),b(:,i))

if(d .lt. thresh**2 ) return
krank  = i
b(:,i) = b(:,i) / sqrt(d)

do j=i+1,m
if( rnorms(j) .lt. thresh) continue
cd  = dot_product(b(:,i),b(:,j))
b(:,j) = b(:,j) - cd*b(:,i)
d      = dot_product(b(:,j),b(:,j))
rnorms(j) = sqrt(d)
end do

end do
return
end subroutine



subroutine linalg0_gspiv_c(eps,n,m,b,krank,rnorms,ipivs)
implicit double precision (a-h,o-z)
double precision :: rnorms(:)
integer          :: ipivs(:)
double complex   :: b(:,:)
!
!  Perform Gram-Schmidt with pivoting and reorthogonalization on the input
!  matrix.
!
double complex              :: cd
double complex              :: work2(n)

done = 1
dtot = 0
do i=1,m
ipivs(i) = i
d          = dot_product(b(:,i),b(:,i))
dtot       = dtot+d
rnorms(i)  = sqrt(d)
end do


thresh=dtot*eps**2 
thresh=sqrt(thresh)
krank = 0

do i=1,min(n,m)

! choose the next column
ipivot = i
rn     = rnorms(i)
do j=i+1,m
if (rnorms(j) .gt. rn) then
rn     = rnorms(j)
ipivot = j
endif
end do

! swap the chosen column with the ith column
work2         = b(1:n,i)
b(1:n,i)      = b(1:n,ipivot)
b(1:n,ipivot) = work2

ii            = ipivs(i)
ipivs(i)      = ipivs(ipivot)
ipivs(ipivot) = ii

d               = rnorms(i)
rnorms(i)       = rnorms(ipivot)
rnorms(ipivot)  = d

! orthogonalize the chosen column to all previous columns
do j=1,i-1
cd = dot_product(b(:,j),b(:,i))
b(:,i) = b(:,i) - cd*b(:,j)
end do

! recompute its norm
d  = dot_product(b(:,i),b(:,i))

if(d .lt. thresh**2 ) return
krank  = i
b(:,i) = b(:,i) / sqrt(d)

do j=i+1,m
if( rnorms(j) .lt. thresh) continue
cd  = dot_product(b(:,i),b(:,j))
b(:,j) = b(:,j) - cd*b(:,i)
d      = dot_product(b(:,j),b(:,j))
rnorms(j) = sqrt(d)
end do

end do
return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Vladimir Rokhlin's routines for computing 2x2 SVDs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine vlad_oneuv_r(a7,rlam,u,v)
        implicit double precision (a-h,o-z)
        dimension u(2,2),v(2,2),u0(2,2),v0(2,2),wl(2,2), &
           uk(2,2),vk(2,2),u2(2,2),v2(2,2),a(2,2),a7(2,2),&
           rlam(2),a2(2,2),wr(2,2)
!        data eps/1.0d-15/
        eps = epsilon(0.0d0)*10
! 
!        initialize the matrices u,v
! 
        uk(1,1)=1
        uk(2,2)=1
        uk(1,2)=0
        uk(2,1)=0
! 
        call vlad_d2mcopy2(uk,vk)
        call vlad_d2mcopy2(a7,a)
! 
        do 2000 i=1,10
! 
!       attempt to construct the SVD of a
! 
        call vlad_d2muv0(ier,a,rlam,u0,v0,eps)
        if(ier .eq. 0) goto 2200
! 
!        the attempt failed. rotate the matrix a a little from both sides
! 
!        . . . construct the perturbing matrices
! 
        if(i .ne. 1) goto 1600
        d=1
        d=d/10
        alpha=dcos(d)
        beta=dsin(d)
! 
        wl(1,1)=alpha
        wl(2,2)=alpha
        wl(1,2)=beta
        wl(2,1)=-beta
! 
        d=1
        d=d/7
        alpha=dcos(d)
        beta=dsin(d)
! 
        wr(1,1)=alpha
        wr(2,2)=alpha
        wr(1,2)=beta
        wr(2,1)=-beta
 1600 continue
! 
!        . . . perturb the matrix from the left
! 
        call vlad_d2mprod2(wl,a,a2)
        wl(1,2)=-wl(1,2)
        wl(2,1)=-wl(2,1)
        call vlad_d2mprod2(uk,wl,u2)
        wl(1,2)=-wl(1,2)
        wl(2,1)=-wl(2,1)
! 
!        . . . perturb the matrix from the right
! 
        call vlad_d2mprod2(a2,wr,a)
        wr(1,2)=-wr(1,2)
        wr(2,1)=-wr(2,1)
        call vlad_d2mprod2(wr,vk,v2)
        wr(1,2)=-wr(1,2)
        wr(2,1)=-wr(2,1)
! 
        call vlad_d2mcopy2(u2,uk)
        call vlad_d2mcopy2(v2,vk)
 2000 continue
 2200 continue
        call vlad_d2mprod2(u0,uk,u)
        call vlad_d2mprod2(vk,v0,v)
  
        return
        end
! 
! 
! 
! 
        subroutine vlad_d2mcopy2(a,b)
        implicit double precision (a-h,o-z)
        dimension a(4),b(4)
        b(1)=a(1)
        b(2)=a(2)
        b(3)=a(3)
        b(4)=a(4)
        return
        end
  
  
! 
! 
! 
! 
! 
        subroutine vlad_d2muv0(ier,a,rlam,u,v,eps2)
        implicit double precision (a-h,o-z)
        dimension a(2,2),u(2,2),v(2,2),rlam(2), &
           b(2,2),w(2,2),vstar(2,2),z(2,2),rlams(2,2) 
!        data eps/1.0d-10/
        eps = epsilon(0.0d0)*100
! 
!       this subroutine produces the singular value decomposition
!       of a 2 * 2 real matrix a, so that on exit,
! 
!       a = u d v,
! 
!       with u and v orthogonal, and d a diagonal matrix with the
!       vector rlam=(rlam(1),rlam(2)) on the diagonal.
! 
        ier=0
! 
!       . . . simmetrize a
! 
        den=a(2,2)+a(1,1)
        rnum=a(2,1)-a(1,2)
        dd=dabs(a(1,1))+dabs(a(2,1))+dabs(a(1,2))+dabs(a(2,2))
        if(dabs(rnum) .gt. eps2*dd) goto 1100
        alpha=1
        beta=0
        goto 1400
 1100 continue
! 
!       if denominator is not too small, use exact formulae
! 
        if(dabs(den) .le. eps*dabs(rnum)) goto 1200
        tgphi=rnum/den
        phi=datan(tgphi)
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1400
 1200 continue
! 
!       denominator is too small. bomb.
! 
        ier=4
        return
 1400 continue
! 
!      calculate the simmetrizing matrix and the simmetric one
! 
        w(1,1)=alpha
        w(1,2)=beta
        w(2,1)=-beta
        w(2,2)=alpha
        call vlad_d2mprod2(w,a,b)
! 
!       now, diagonalize the simmetrized matrix
! 
        den=b(2,2)-b(1,1)
        rnum=b(1,2)*2
        if(dabs(rnum) .gt. eps2*dd) goto 1500
        alpha=1
        beta=0
        goto 1800
 1500 continue
! 
!       if denominator is not too small, use exact formulae
! 
        if(dabs(den) .le. eps*dabs(rnum)) goto 1600
        tg2phi=-rnum/den
        phi=datan(tg2phi)/2
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1800
 1600 continue
! 
!       denominator is too small. bomb.
! 
        ier=8
        return
 1800 continue

! 
!       construct the diagonalizing matrix
!       and the resulting diagonal
! 
        v(1,1)=alpha
        v(1,2)=beta
        v(2,1)=-beta
        v(2,2)=alpha
! 
        vstar(1,1)=v(1,1)
        vstar(1,2)=v(2,1)
        vstar(2,1)=v(1,2)
        vstar(2,2)=v(2,2)
! 
        call vlad_d2mprod2(v,w,u)
! 
!       finally, compute the resulting diagonal elements
! 
        call vlad_d2mprod2(u,a,z)
        call vlad_d2mprod2(z,vstar,rlams)
        rlam(1)=rlams(1,1)
        rlam(2)=rlams(2,2)
        u(1,2)=-u(1,2)
        u(2,1)=-u(2,1)
        return
        end

        subroutine vlad_d2mprod2(a,b,c)
        implicit double precision (a-h,o-z)
        dimension a(2,2),b(2,2),c(2,2)
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end







        subroutine vlad_oneuv_c(a0,rlam,u,v)
        implicit double precision (a-h,o-z)
        double complex z,b(2,2),u(2,2),a0(2,2),uu(2,2),v(2,2),w1(2,2),w2(2,2), &
           w3(2,2),cd,clam1,clam2,clam,x(2),y(2)
        dimension rlam(2)
! 
!       multiply a by its adjoint
! 
        b = matmul(transpose(conjg(a0)),a0)
! 
!       find the eigenvalues of a^*a
! 
        two=2
        four=4
        cd=(b(1,1)-b(2,2))**2+four*b(1,2)*b(2,1)
        cd=sqrt(cd)
        clam1= (b(1,1)+b(2,2)+cd )/two
        clam2= (b(1,1)+b(2,2)-cd )/two
! 
!        find the eigenvectors
! 
        clam=clam1
        dc2=clam2*conjg(clam2)
        dc1=clam1*conjg(clam1)
        if( dc2 .gt. dc1 ) clam=clam2
! 
        d11=(b(1,1)-clam)*conjg((b(1,1)-clam))
        d22=(b(2,2)-clam)*conjg((b(2,2)-clam))
        d12=b(1,2)*conjg(b(2,2))
! 
        if( (d11 .lt. d22) .or. (d11 .lt. d12) ) goto 1400
! 
        x(2)=1
        x(1)=-b(1,2)/(b(1,1)-clam)
        goto 1800
! 
 1400 continue
! 
        if(d22 .lt. d12) goto 1600
! 
        x(1)=1
        x(2)=-b(2,1)/(b(2,2)-clam)
        goto 1800
! 
 1600 continue
! 
        x(2)=1
        x(1)=(b(1,1)-clam)/b(1,2)
! 
 1800 continue
! 
        d=x(1)*conjg(x(1))+x(2)*conjg(x(2))
        d=sqrt(d)
        x(1)=x(1)/d
        x(2)=x(2)/d
! 
        call  vlad_d2md2ort(x,y)
! 
        uu(1,1)=x(1)
        uu(2,1)=x(2)
        uu(1,2)=y(1)
        uu(2,2)=y(2)
!
        w2 = transpose(conjg(uu))
        v = matmul(a0,uu)
!        call  vlad_d2mcpcon(uu,w2)
!         call  vlad_d2mprodc(a0,uu,v)

! 
!        now, obtain v
! 
! 
         d1=v(1,1)*conjg(v(1,1))+v(2,1)*conjg(v(2,1))
         d2=v(1,2)*conjg(v(1,2))+v(2,2)*conjg(v(2,2))
! 
         if(d2 .gt. d1) goto 3000
! 
         d=sqrt(d1)
         call  vlad_d2md2ort(v(1,1),v(1,2))
         goto 3200
! 
 3000 continue
! 
         d=sqrt(d2)
         call  vlad_d2md2ort(v(1,2),v(1,1))
 3200 continue
! 
!       obtain the diagonal matrix
! 
!         call  vlad_d2mprodc(a0,uu,w1)
         w1 = matmul(a0,uu)
!         call  vlad_d2mcpcon(v,b)
         b = transpose(conjg(v))

!         call  vlad_d2mprodc(b,w1,w3)
         w3 = matmul(b,w1)
! 
!        call  vlad_d2mcopy2(w2,u)
         u = w2
! 
!       finally, make sure that the diagonal elements are real
! 
        z=1
        dzero=0
        d=w3(1,1)*conjg(w3(1,1))
        if(d .ne. dzero) z=w3(1,1)/abs(w3(1,1))
        rlam(1)=abs(w3(1,1))
        v(1,1)=v(1,1)*z
        v(2,1)=v(2,1)*z
! 
        z=1
        dzero=0
        d=w3(2,2)*conjg(w3(2,2))
        if(d .ne. dzero) z=w3(2,2)/abs(w3(2,2))
        rlam(2)=abs(w3(2,2))
        v(1,2)=v(1,2)*z
        v(2,2)=v(2,2)*z
! 
!        invert u, v
! 
!        call vlad_d2mcopy2(v,w3)
        w3 = v
!        call vlad_d2mcpcon(u,v)
        v = conjg(transpose(u))
!        call vlad_d2mcpcon(w3,u)
        u = conjg(transpose(w3))
        return
        end

        subroutine vlad_d2md2ort(x,y)
        implicit double precision (a-h,o-z)
        double complex x(2),y(2)
! 
!       normalize x
! 
        d1=x(1)*conjg(x(1))
        d2=x(2)*conjg(x(2))
        d=d1+d2
        d=sqrt(d)
        x(1)=x(1)/d
        x(2)=x(2)/d
! 
!       make y into a unit vector orthogonal to x
! 
!       . . . when |x(1)| > |x(2)|
! 
        if(d1 .lt. d2) goto 1200
        y(2)=x(1)
        y(1)=-conjg(x(2))*x(1)/conjg(x(1))
        return
 1200 continue
! 
!       . . . when |x(2)| > |x(1)|
! 
        y(1)=x(2)
        y(2)=-conjg(x(1))*x(2)/conjg(x(2))
        return
        end


end module

