!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for solving the Riccati equation
!
!    r'(t) + (r(t))^2 + q(t) = 0,  a <= t <= b,                                           (1)
!
!  satisfied by the logarithmic derivatives of solutions of the normal form
!
!    y''(t) + q(t) y(t) = 0,   a <= t < = b,                                              (2)
!
!  of the Helmoholtz equation. The advantage of solving (1) over (2) is that
!  the solutions of the former are easier to represent than those of the latter
!  in many cases of interest.  This is particularly true in the oscillatory regime, 
!  where the coefficient q(t) >=0, and the code here is intended for use in that
!  regime.  
!
!  The relationship between the solutions of (2) and (1) can be expressed easily
!  in terms of phase functions.  If alpha'(t) = Im(r(t)) then
!
!               cos ( alpha(t) )                  sin ( alpha(t ) ) 
!    u(t) =     ------------------  and   v(t) =  -----------------                       (3)
!               sqrt( alpha'(t) )                 sqrt( alpha'(t) )
!
!
!  form a basis in the space of solutions of (2).   The function alpha(t) is
!  known as a phase function for (2).  Note that the real part of r(t) is 
!   
!           -  alpha''(t)
!            -------------.
!             2 alpha'(t)
!
!  It should be noted that this code can be used to construct phase functions
!  which represent solutions of (2) when q(t) is negative, but numerical difficulties 
!  are encountered when alpha'(t) becomes on the order of machine precision.  
!  Appell's equation (see appell_nop.f90) is much better in this regard.
!
!  The coefficient in the equation can be specified as piecewise Chebyshev expansions
!  or via an external subroutine.
!
!  The following subroutine should be regarded as publicly callable:
!
!    riccati_nop_init - initialize the solver code -- this entails populating a 
!     structure which is then passed to the other codes in this file
!
!    riccati_nop_window - find the value of the first two derivatives of a 
!      nonoscillatory phase function alpha(t) at the left or right endpoint
!      of interval [a,b] on which Q>>0
!
!    riccati_nop_phase -  given the values of alpha'(t) and alpha''(t) at a point c
!      in an interval [a,b] and the value of alpha(t) at a point d in [a,b],
!      solve Riccati's equation in order to construct piecewise Chebyshev expansions
!      of alpha(t) and its first two derivatives
!
!    riccati_nop_extend -  given piecewise Chebyshev expansions of the phase function
!      alpha(t) and its first two derivatives on an interval [a,b], try to extend
!      alpha(t) on one or both sides of the the interval [a,b] an interval
!
!    riccati_nop_solve - construct a complex-valued solution of (1) with a specified 
!      value at a point c in the interval [a,b]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Observations / To do:
!
!    - When  the p(t) y'(t) term is not present in the equation, Appell's equation
!      should generally be used to construct phase functions rather than the solver
!      here.
!
!      When p(t) is present, however, Appell's eqution runs into difficulties because
!      the Wronskian is not constant (see appell_nop.f90 for more discussion of this).
!      In that case, solving Riccati's equation is worthwhile and the code in
!      riccati.f90 can be used to do so.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module riccati_nop

use utils
use linalg0
use chebyshev
use chebpw
use ieee_arithmetic


type      riccati_nop_vars_t
integer                                 :: kcheb
integer                                 :: ntest
integer                                 :: maxsteps  ! max refinement steps for Newton
integer                                 :: maxiters  ! max iters for trapezoidal rule
integer                                 :: maxints   ! max intervals for solutions
double precision                        :: epsint    ! smallest allowed interval size
                                        
double precision, allocatable           :: xscheb(:)
double precision, allocatable           :: whtscheb(:)
double precision, allocatable           :: acoefs(:,:)
double precision, allocatable           :: aintl(:,:)
double precision, allocatable           :: adiff(:,:)
double precision, allocatable           :: aintr(:,:)
end type  riccati_nop_vars_t



interface
subroutine riccati_nop_fun1(n,ts,qs,par1,par2,par3)
implicit double precision (a-h,o-z)
double precision :: ts(:), qs(:)
end subroutine
end interface

interface     riccati_nop_solve
module procedure riccati_nop_solve1
module procedure riccati_nop_solve2
end interface riccati_nop_solve

interface     riccati_nop_phase
module procedure riccati_nop_phase1
module procedure riccati_nop_phase2
end interface riccati_nop_phase

interface     riccati_nop_extend
module procedure riccati_nop_extend1
module procedure riccati_nop_extend2
end interface riccati_nop_extend

interface     riccati_nop_window
module procedure riccati_nop_window1
module procedure riccati_nop_window2
end interface riccati_nop_window

contains


subroutine riccati_nop_init(vars,kcheb0,ntest0,maxints0,epsint0,maxsteps0,maxiters0)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t), intent(out)      :: vars
integer, optional                          :: kcheb0, ntest0, maxints0, maxsteps0, maxiters0
double precision, optional                 :: epsint0
!
!  Initialize the structure containing any data needed by the other
!  procedures in this module.
!
!  If this subroutine is called with no other arguments than vars, then
!  reasonable defaults are chosen.
!
!  Input parameters:
!    kcheb0 - the number of terms in the Chebyshev expansions used to represent solutions
!    ntest0 - the number of trailing Chebyshev coefficients which must be small in order
!      to consider a Chebyshev expansion accurate
!    maxints0 - the maximum number of intervals used to represent solutions
!    epsint0 - the smallest allowable interval when adaptively discretizating
!    maxsteps0 - the maximum number of refinement steps
!    maxiters0 - the maximum number of iterations for the trapezoidal method
!
!  Output parameters:
!    N/A
!

if (.not. present(kcheb0)) then
kcheb    = 16
ntest    = 4
maxiters = 4
maxsteps = 4
maxints  = 1000
epsint   = 1.0d-7
else
ntest    = ntest0
kcheb    = kcheb0
maxiters = maxiters0
maxsteps = maxsteps0
maxints  = maxints0
epsint   = epsint0
endif

vars%kcheb    = kcheb
vars%ntest    = ntest
vars%maxiters = maxiters
vars%maxints  = maxints
vars%maxsteps = maxsteps
vars%epsint   = epsint

call chebyshev_quad(kcheb,vars%xscheb,vars%whtscheb)
call chebyshev_coefsmatrix(kcheb,vars%acoefs)
call chebyshev_intlmatrix(kcheb,vars%aintl)
call chebyshev_intrmatrix(kcheb,vars%aintr)
call chebyshev_diffmatrix(kcheb,vars%adiff)

end subroutine


subroutine riccati_nop_fitr(vars,vals,dcoefs)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)             :: vars
double precision                     :: vals(:), coefs(vars%kcheb)
kcheb  = vars%kcheb
ntest  = vars%ntest
coefs  = matmul(vars%acoefs,vals)
dd1    = norm2(abs(coefs))
dd2    = norm2(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1=1
dcoefs = dd2/dd1
end subroutine

subroutine riccati_nop_fitc(vars,vals,dcoefs)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)             :: vars
double complex                       :: vals(:), coefs(vars%kcheb)
kcheb  = vars%kcheb
ntest  = vars%ntest
coefs  = matmul(vars%acoefs,vals)
dd1    = norm2(abs(coefs))
dd2    = norm2(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1=1
dcoefs = dd2/dd1
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Extension code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine riccati_nop_extend1(vars,ier,eps,amax,bmax,chebcoefs,qs,aout,bout,chebphase, &
  alpha,alphader,alphader2)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)                     :: vars
type(chebpw_scheme)                          :: chebcoefs
double precision                             :: qs(:)

type(chebpw_scheme)                          :: chebphase
double precision, allocatable, intent(inout) :: alpha(:), alphader(:), alphader2(:)
!
!  This routine attempts to extend an existing phase function alpha(t) given on an interval [a,b].
!  If the parameter amax < a, then it attempts to find a point aout in [amax,a] such that
!
!    alpha(aout) = k*pi/16
!
!  with k an integer and, assuming this is achieved, the phase function alpha(t) is extended to
!  an interval containing aout.
!
!  Similarly, if bmax > b, then it attempts to find a point bout in [b,bmax] such that
!
!    alpha(bout) = k*pi/16
!
!  with k an integer and, assuming this is achieved, the phase function alpha(t) is extended
!  to an interval containing bout.
!
!  The coefficient q(t) is specified via a piecewise Chebyshev expansion.
!
!  Input parameters:
!    eps - the precision for the calculations
!    [amax,bmax] - the largest possible interval
!    chebcoefs - a piecewise Chebyshev discretization scheme giving the coefficients on an
!      interval containing [amax,bmax]
!   qs - the values of the coefficient at the discretization nodes
!
!  chebphase - a  structure describing a discretization of the phase function alpha(t)
!  alpha, alphader, alphader2 - the values of alpha(t) and its first two derivatives
!    at the discretization nodes
!
!  Output parameters:
!   ier - an error return code;
!     ier = 0    means successful execution
!     ier = 64   means that the code failed to extend the solution to the left
!     ier = 128  means that the bisection method failed to find the point aout
!     ier = 256   means that the code failed to extend the solution to the right
!     ier = 512  means that the bisection method failed to find the point bout
!   
!
!  (aout, bout) - the values of aout and bout
!  chebphase - this structure is updated to reflect the extension of the phase function
!  alpha, alphader, alphader2 - the values of the extended alpha(t) and its first 
!    two derivatives at the discretization nodes
!
!
double precision, allocatable              :: stack(:,:), ts0(:)
double precision, allocatable              :: qs0(:)
double complex, allocatable                :: rs0(:), rders0(:)

double precision, allocatable              :: ab0(:,:)
double precision, allocatable              :: alphapp0(:,:),alphap0(:,:), alpha0(:,:)

double precision, allocatable              :: abnew(:,:)
double precision, allocatable              :: alphanew(:), alphadernew(:), alphader2new(:)

double complex                             :: rval, ima

data pi / 3.14159265358979323846264338327950288d0 /

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
maxiters = vars%maxiters
maxsteps = vars%maxsteps


ier      = 0
eps0     = epsilon(0.0d0)
ima      = (0.0d0,1.0d0)
maxstack = maxints*2
dmul     = pi/16
call chebpw_interval(chebphase,a,b)

aout     = a
bout     = b

allocate(stack(2,maxstack))

allocate(ts0(kcheb), qs0(kcheb) )
allocate(rs0(kcheb),rders0(kcheb))
allocate(ab0(2,-maxints:maxints))

allocate(alpha0(kcheb,-maxints:maxints))
allocate(alphap0(kcheb,-maxints:maxints))
allocate(alphapp0(kcheb,-maxints:maxints))

nints1   = 0
nints2   = 0

!
!  Try to extend the interval to the right if bmax > b
!


if (bmax .gt. b) then

   call chebpw_interp(chebphase,alpha,alphader,alphader2,b,aval,apval,appval)
   rval = -0.5d0 * appval / apval + ima* apval
   dtarget  = ceiling(aval/dmul)*dmul

   a0       = b
   ind      = 0
   dstep    = (bmax-b)/32
   icount   = 0             ! used to count the number of consecutive successful steps
   ifdone   = 0
    
   do while (a0 .lt. bmax)
   b0    = min(bmax,a0+dstep)

   !
   !  Use a spectral integration method to take a step
   !
    
   ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
   call chebpw_interp(chebcoefs,qs,kcheb,ts0,qs0)
   ifaccept = 0
   dq       = 1d300
   dcoefs   = 1d300
   jer      = -1

   call riccati_nop_fitr(vars,qs0,dq)

   if (dq .lt. eps) then

   call riccati_nop_ivp_specint_c(vars,jer,eps,a0,b0,kcheb,ts0,qs0,rval,rs0,rders0)

   if (jer .eq. 0) then
   alphapp0(:,nints2+1) = imag(rders0)
   alphap0(:,nints2+1)  = imag(rs0)
   alpha0(:,nints2+1)   = aval + (b0-a0)/2*matmul(vars%aintl,alphap0(:,nints2+1))
   
   call riccati_nop_fitr(vars,alphap0(:,nints2+1),dcoefs)

   if (dcoefs .lt. eps) ifaccept = 1
   endif

   endif
     
   ! ind = ind+1
   ! write(*,"(2(I6,1X),5(D15.8,1X),I4)") ind,jer,a0,b0,dstep,dq,dcoefs,ifaccept
   ! write(13,"(2(I6,1X),5(D15.8,1X),I4)") ind,jer,a0,b0,dstep,dq,dcoefs,ifaccept

   !
   !  If the interval is acceptable, take a step
   !
    
   if (ifaccept .eq. 1) then
    
     if (nints2+1 .ge. maxints) then
     ier = 256
     return
     endif
    
     ! add the interval to the list of intervals and update the inital value for
     ! the next step
     nints2           = nints2+1
     ab0(1,nints2)    = a0
     ab0(2,nints2)    = b0     
     rval             = rs0(kcheb)
     aval             = alpha0(kcheb,nints2)
    
     a0               = b0
     icount           = icount+1
    
     ! adjust the step size
     if (icount .ge. 2) then
     dstep         = dstep*2
     endif

     ! check if we have achieved the target value
     if (aval .gt. dtarget) then
     ifdone = 1
     exit
     endif
    
   endif
    
   !
   !  If the interval is not accepted, reset the "success" count to zero and
   !  half the step size
   !
   if (ifaccept .eq. 0) then
     dstep         = dstep/2
     icount        = 0
   endif
    
   end do

  if (ifdone .eq. 0) then
  ier = 256
  return
  endif

endif



!
!  Try to extend the interval to the left if amax < a and a != d
!

if (amax .lt. a) then

   call chebpw_interp(chebphase,alpha,alphader,alphader2,a,aval,apval,appval)
   rval = -0.5d0 * appval/apval + ima*apval
   dtarget  = floor(aval/dmul)*dmul

   b0       = a
   ind      = 0
   dstep    = (a-amax)/32
   icount   = 0             ! used to count the number of consecutive successful steps
   ifdone   = 0   
   
   do while (b0 .gt. amax)
    
   a0    = max(amax,b0-dstep)

   ifaccept = 0
   dcoefs   = 1d300
   dq       = 1d300
   jer      = -1
    
   !
   !  Use a spectral integration method to take a step
   !
    
   ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
   call chebpw_interp(chebcoefs,qs,kcheb,ts0,qs0)

   call riccati_nop_fitr(vars,qs0,dq)

   if (dq .lt. eps) then

   call riccati_nop_tvp_specint_c(vars,jer,eps,a0,b0,kcheb,ts0,qs0,rval,rs0,rders0)
   if (jer .eq. 0) then
   alphapp0(:,nints1-1) = imag(rders0)
   alphap0(:,nints1-1)  = imag(rs0)
   alpha0(:,nints1-1)   = aval + (b0-a0)/2*matmul(vars%aintr,alphap0(:,nints1-1))


   call riccati_nop_fitr(vars,alphap0(:,nints1-1),dcoefs)
   if (dcoefs .lt. eps) ifaccept = 1
   endif
   endif
    
   ! ind = ind+1
   ! write(*,"(2(I6,1X),5(D15.8,1X),I4)") ind,jer,a0,b0,dstep,dq,dcoefs,ifaccept
   ! write(13,"(2(I6,1X),5(D15.8,1X),I4)") ind,jer,a0,b0,dstep,dq,dcoefs,ifaccept

   !
   !  If the interval is acceptable, take a step
   !
    
   if (ifaccept .eq. 1) then
    
     if (-nints1+1 .ge. maxints) then
     ier = 64
     return
     endif
    
     ! add the interval to the list of intervals and update the inital value for
     ! the next step
     nints1           = nints1-1
     ab0(1,nints1)    = a0
     ab0(2,nints1)    = b0

     rval             = rs0(1)
     aval             = alpha0(1,nints1)
    
     b0               = a0
     icount           = icount+1
    
     ! adjust the step size
     if (icount .ge. 2) then
     dstep         = dstep*2
     endif
     
     if (aval .lt. dtarget) then
     ifdone = 1
     exit
     endif

   endif
    
   !
   !  If the interval is not accepted, reset the "success" count to zero and
   !  half the step size
   !
   if (ifaccept .eq. 0) then
     dstep         = dstep/2
     icount        = 0
   endif
   end do


if (ifdone .eq. 0) then
ier = 64
return
endif

endif

!
!  Update the discretization structure and the solutions
!


nints    = chebphase%nints
nintsnew = nints-nints1+nints2

allocate(abnew(2,nintsnew))
allocate(alphanew(kcheb*nintsnew))
allocate(alphadernew(kcheb*nintsnew))
allocate(alphader2new(kcheb*nintsnew))

idx = 0
i1  = 1
i2  = kcheb

do int=nints1,-1
a0 = ab0(1,int)
b0 = ab0(2,int)
idx = idx +1
abnew(1,idx) = a0
abnew(2,idx) = b0

alphanew(i1:i2)     = alpha0(:,int)
alphadernew(i1:i2)  = alphap0(:,int)
alphader2new(i1:i2) = alphapp0(:,int)

i1              = i1+kcheb
i2              = i2+kcheb

end do

j1 = 1
j2 = kcheb

do int=1,nints
a0 = chebphase%ab(1,int)
b0 = chebphase%ab(2,int)
idx = idx +1
abnew(1,idx) = a0
abnew(2,idx) = b0

alphanew(i1:i2)     = alpha(j1:j2)
alphadernew(i1:i2)  = alphader(j1:j2)
alphader2new(i1:i2) = alphader2(j1:j2)

i1              = i1+kcheb
i2              = i2+kcheb
j1              = j1+kcheb
j2              = j2+kcheb

end do

j1 = 1

do int=1,nints2
a0 = ab0(1,int)
b0 = ab0(2,int)
idx = idx +1
abnew(1,idx) = a0
abnew(2,idx) = b0

alphanew(i1:i2)     = alpha0(:,int)
alphadernew(i1:i2)  = alphap0(:,int)
alphader2new(i1:i2) = alphapp0(:,int)

i1              = i1+kcheb
i2              = i2+kcheb

end do


call chebpw_specified(chebphase,kcheb,nintsnew,abnew)
nints =nintsnew

deallocate(alpha,alphader,alphader2)

n = nints*kcheb

allocate(alpha(n),alphader(n),alphader2(n))
alpha     = alphanew
alphader  = alphadernew
alphader2 = alphader2new



!
!  Use bisection to find aout and bout
!

!
! Use bisection to try to find aout and bout
!

maxiters = 200

if (bmax .gt. b) then

call chebpw_interp(chebphase,alpha,b,val)
valout = ceiling(val/dmul)*dmul

xleft  = b
xright = chebphase%ab(2,nints)


do iter=1,maxiters

bout = (xleft+xright)/2
call chebpw_interp(chebphase,alpha,bout,val)

if (val .lt. valout) then
xleft  = bout
else
xright = bout
endif

!print *,iter,xleft,bout,xright,val-valout

if (abs(xleft-xright) .lt. eps0) exit

end do

if (iter .gt. maxiters) then
ier = 512
return
endif

endif



if (a .gt. amax) then

call chebpw_interp(chebphase,alpha,a,val)
valout = floor(val/dmul)*dmul

xleft  = chebphase%ab(1,1)
xright = a


do iter=1,maxiters

aout = (xleft+xright)/2
call chebpw_interp(chebphase,alpha,aout,val)

if (val .lt. valout) then
xleft  = aout
else
xright = aout
endif

!print *,iter,xleft,aout,xright,val-valout

if (abs(xleft-xright) .lt. eps0) exit

!if ( abs(val-valout) .lt. eps*abs(valout)) exit

end do


if (iter .gt. maxiters) then
ier = 128
return
endif

endif


end subroutine



subroutine riccati_nop_extend2(vars,ier,eps,amax,bmax,fun,par1,par2,par3,aout,bout,chebphase, &
  alpha,alphader,alphader2)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)                     :: vars
procedure(riccati_nop_fun1)                  :: fun
type(chebpw_scheme)                          :: chebphase
double precision, allocatable, intent(inout) :: alpha(:), alphader(:), alphader2(:)
!
!  This routine attempts to extend an existing phase function alpha(t) given on an interval [a,b].
!  If the parameter amax < a, then it attempts to find a point aout in [amax,a] such that
!
!    alpha(aout) = k*pi/16
!
!  with k an integer and, assuming this is achieved, the phase function alpha(t) is extended to
!  an interval containing aout.
!
!  Similarly, if bmax > b, then it attempts to find a point bout in [b,bmax] such that
!
!    alpha(bout) = k*pi/16
!
!  with k an integer and, assuming this is achieved, the phase function alpha(t) is extended
!  to an interval containing bout.
!
!  The coefficient q(t) is specified via an external subroutine.
!
!  Input parameters:
!    eps - the precision for the calculations
!    [amax,bmax] - the largest possible interval
!
!   chebphase - a  structure describing a discretization of the phase function alpha(t)
!   alpha, alphader, alphader2 - the values of alpha(t) and its first two derivatives
!     at the discretization nodes
!   fun - an exernal subroutine supplying the values of q(t)
!   par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!   ier - an error return code;
!     ier = 0    means successful execution
!     ier = 64   means that the code failed to extend the solution to the left
!     ier = 128  means that the bisection method failed to find the point aout
!     ier = 256   means that the code failed to extend the solution to the right
!     ier = 512  means that the bisection method failed to find the point bout
!   
!
!  (aout, bout) - the values of aout and bout
!  chebphase - this structure is updated to reflect the extension of the phase function
!  alpha, alphader, alphader2 - the values of the extended alpha(t) and its first 
!    two derivatives at the discretization nodes
!
!
double precision, allocatable              :: stack(:,:), ts0(:)
double precision, allocatable              :: qs0(:)
double complex, allocatable                :: rs0(:), rders0(:)

double precision, allocatable              :: ab0(:,:)
double precision, allocatable              :: alphapp0(:,:),alphap0(:,:), alpha0(:,:)

double precision, allocatable              :: abnew(:,:)
double precision, allocatable              :: alphanew(:), alphadernew(:), alphader2new(:)

double complex                             :: rval, ima

data pi / 3.14159265358979323846264338327950288d0 /

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
maxiters = vars%maxiters
maxsteps = vars%maxsteps

ier      = 0
eps0     = epsilon(0.0d0)
ima      = (0.0d0,1.0d0)
maxstack = maxints*2
dmul     = pi/16
call chebpw_interval(chebphase,a,b)

aout     = a
bout     = b

allocate(stack(2,maxstack))

allocate(ts0(kcheb), qs0(kcheb) )
allocate(rs0(kcheb),rders0(kcheb))
allocate(ab0(2,-maxints:maxints))

allocate(alpha0(kcheb,-maxints:maxints))
allocate(alphap0(kcheb,-maxints:maxints))
allocate(alphapp0(kcheb,-maxints:maxints))

nints1   = 0
nints2   = 0

!
!  Try to extend the interval to the right if bmax > b
!


if (bmax .gt. b) then

   call chebpw_interp(chebphase,alpha,alphader,alphader2,b,aval,apval,appval)
   rval = -0.5d0 * appval / apval + ima* apval
   dtarget  = ceiling(aval/dmul)*dmul

   a0       = b
   ind      = 0
   dstep    = (bmax-b)/32
   icount   = 0             ! used to count the number of consecutive successful steps
   ifdone   = 0
    
   do while (a0 .lt. bmax)
   b0    = min(bmax,a0+dstep)

   !
   !  Use a spectral integration method to take a step
   !
    
   ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
   call fun(kcheb,ts0,qs0,par1,par2,par3)
!   call chebpw_interp(chebcoefs,qs,kcheb,ts0,qs0)
   ifaccept = 0
   dq       = 1d300
   dcoefs   = 1d300
   jer      = -1

   call riccati_nop_fitr(vars,qs0,dq)

   if (dq .lt. eps) then

   call riccati_nop_ivp_specint_c(vars,jer,eps,a0,b0,kcheb,ts0,qs0,rval,rs0,rders0)

   if (jer .eq. 0) then
   alphapp0(:,nints2+1) = imag(rders0)
   alphap0(:,nints2+1)  = imag(rs0)
   alpha0(:,nints2+1)   = aval + (b0-a0)/2*matmul(vars%aintl,alphap0(:,nints2+1))
   
   call riccati_nop_fitr(vars,alphap0(:,nints2+1),dcoefs)

   if (dcoefs .lt. eps) ifaccept = 1
   endif

   endif
     
   ! ind = ind+1
   ! write(*,"(2(I6,1X),5(D15.8,1X),I4)") ind,jer,a0,b0,dstep,dq,dcoefs,ifaccept
   ! write(13,"(2(I6,1X),5(D15.8,1X),I4)") ind,jer,a0,b0,dstep,dq,dcoefs,ifaccept

   !
   !  If the interval is acceptable, take a step
   !
    
   if (ifaccept .eq. 1) then
    
     if (nints2+1 .ge. maxints) then
     ier = 256
     return
     endif
    
     ! add the interval to the list of intervals and update the inital value for
     ! the next step
     nints2           = nints2+1
     ab0(1,nints2)    = a0
     ab0(2,nints2)    = b0     
     rval             = rs0(kcheb)
     aval             = alpha0(kcheb,nints2)
    
     a0               = b0
     icount           = icount+1
    
     ! adjust the step size
     if (icount .ge. 2) then
     dstep         = dstep*2
     endif

     ! check if we have achieved the target value
     if (aval .gt. dtarget) then
     ifdone = 1
     exit
     endif
    
   endif
    
   !
   !  If the interval is not accepted, reset the "success" count to zero and
   !  half the step size
   !
   if (ifaccept .eq. 0) then
     dstep         = dstep/2
     icount        = 0
   endif
    
   end do

  if (ifdone .eq. 0) then
  ier = 256
  return
  endif

endif



!
!  Try to extend the interval to the left if amax < a and a != d
!

if (amax .lt. a) then

   call chebpw_interp(chebphase,alpha,alphader,alphader2,a,aval,apval,appval)
   rval = -0.5d0 * appval/apval + ima*apval
   dtarget  = floor(aval/dmul)*dmul

   b0       = a
   ind      = 0
   dstep    = (a-amax)/32
   icount   = 0             ! used to count the number of consecutive successful steps
   ifdone   = 0   
   
   do while (b0 .gt. amax)
    
   a0    = max(amax,b0-dstep)

   ifaccept = 0
   dcoefs   = 1d300
   dq       = 1d300
   jer      = -1
    
   !
   !  Use a spectral integration method to take a step
   !
    
   ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
   call fun(kcheb,ts0,qs0,par1,par2,par3)

!   call chebpw_interp(chebcoefs,qs,kcheb,ts0,qs0)

   call riccati_nop_fitr(vars,qs0,dq)

   if (dq .lt. eps) then

   call riccati_nop_tvp_specint_c(vars,jer,eps,a0,b0,kcheb,ts0,qs0,rval,rs0,rders0)
   if (jer .eq. 0) then
   alphapp0(:,nints1-1) = imag(rders0)
   alphap0(:,nints1-1)  = imag(rs0)
   alpha0(:,nints1-1)   = aval + (b0-a0)/2*matmul(vars%aintr,alphap0(:,nints1-1))


   call riccati_nop_fitr(vars,alphap0(:,nints1-1),dcoefs)
   if (dcoefs .lt. eps) ifaccept = 1
   endif
   endif
    
   ! ind = ind+1
   ! write(*,"(2(I6,1X),5(D15.8,1X),I4)") ind,jer,a0,b0,dstep,dq,dcoefs,ifaccept
   ! write(13,"(2(I6,1X),5(D15.8,1X),I4)") ind,jer,a0,b0,dstep,dq,dcoefs,ifaccept

   !
   !  If the interval is acceptable, take a step
   !
    
   if (ifaccept .eq. 1) then
    
     if (-nints1+1 .ge. maxints) then
     ier = 64
     return
     endif
    
     ! add the interval to the list of intervals and update the inital value for
     ! the next step
     nints1           = nints1-1
     ab0(1,nints1)    = a0
     ab0(2,nints1)    = b0

     rval             = rs0(1)
     aval             = alpha0(1,nints1)
    
     b0               = a0
     icount           = icount+1
    
     ! adjust the step size
     if (icount .ge. 2) then
     dstep         = dstep*2
     endif
     
     if (aval .lt. dtarget) then
     ifdone = 1
     exit
     endif

   endif
    
   !
   !  If the interval is not accepted, reset the "success" count to zero and
   !  half the step size
   !
   if (ifaccept .eq. 0) then
     dstep         = dstep/2
     icount        = 0
   endif
   end do


if (ifdone .eq. 0) then
ier = 64
return
endif

endif

!
!  Update the discretization structure and the solutions
!


nints    = chebphase%nints
nintsnew = nints-nints1+nints2

allocate(abnew(2,nintsnew))
allocate(alphanew(kcheb*nintsnew))
allocate(alphadernew(kcheb*nintsnew))
allocate(alphader2new(kcheb*nintsnew))

idx = 0
i1  = 1
i2  = kcheb

do int=nints1,-1
a0 = ab0(1,int)
b0 = ab0(2,int)
idx = idx +1
abnew(1,idx) = a0
abnew(2,idx) = b0

alphanew(i1:i2)     = alpha0(:,int)
alphadernew(i1:i2)  = alphap0(:,int)
alphader2new(i1:i2) = alphapp0(:,int)

i1              = i1+kcheb
i2              = i2+kcheb

end do

j1 = 1
j2 = kcheb

do int=1,nints
a0 = chebphase%ab(1,int)
b0 = chebphase%ab(2,int)
idx = idx +1
abnew(1,idx) = a0
abnew(2,idx) = b0

alphanew(i1:i2)     = alpha(j1:j2)
alphadernew(i1:i2)  = alphader(j1:j2)
alphader2new(i1:i2) = alphader2(j1:j2)

i1              = i1+kcheb
i2              = i2+kcheb
j1              = j1+kcheb
j2              = j2+kcheb

end do

j1 = 1

do int=1,nints2
a0 = ab0(1,int)
b0 = ab0(2,int)
idx = idx +1
abnew(1,idx) = a0
abnew(2,idx) = b0

alphanew(i1:i2)     = alpha0(:,int)
alphadernew(i1:i2)  = alphap0(:,int)
alphader2new(i1:i2) = alphapp0(:,int)

i1              = i1+kcheb
i2              = i2+kcheb

end do


call chebpw_specified(chebphase,kcheb,nintsnew,abnew)
nints =nintsnew

deallocate(alpha,alphader,alphader2)

n = nints*kcheb

allocate(alpha(n),alphader(n),alphader2(n))
alpha     = alphanew
alphader  = alphadernew
alphader2 = alphader2new



!
!  Use bisection to find aout and bout
!

!
! Use bisection to try to find aout and bout
!

maxiters = 200

if (bmax .gt. b) then

call chebpw_interp(chebphase,alpha,b,val)
valout = ceiling(val/dmul)*dmul

xleft  = b
xright = chebphase%ab(2,nints)


do iter=1,maxiters

bout = (xleft+xright)/2
call chebpw_interp(chebphase,alpha,bout,val)

if (val .lt. valout) then
xleft  = bout
else
xright = bout
endif

!print *,iter,xleft,bout,xright,val-valout

if (abs(xleft-xright) .lt. eps0) exit

end do

if (iter .gt. maxiters) then
ier = 512
return
endif

endif



if (a .gt. amax) then

call chebpw_interp(chebphase,alpha,a,val)
valout = floor(val/dmul)*dmul

xleft  = chebphase%ab(1,1)
xright = a


do iter=1,maxiters

aout = (xleft+xright)/2
call chebpw_interp(chebphase,alpha,aout,val)

if (val .lt. valout) then
xleft  = aout
else
xright = aout
endif

!print *,iter,xleft,aout,xright,val-valout

if (abs(xleft-xright) .lt. eps0) exit

!if ( abs(val-valout) .lt. eps*abs(valout)) exit

end do


if (iter .gt. maxiters) then
ier = 128
return
endif

endif


end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Phase and windowing code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine riccati_nop_window1(vars,ier,ifleft,eps,a,b,chebcoefs,qs,apval,appval)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)                           :: vars
type(chebpw_scheme)                                :: chebcoefs
double precision                                   :: qs(:)
!
!   Use the "windowing algorithm" to compute the value of the first two derivatives
!   of the nonoscillatory phase function alpha(t) at the left or right endpoint of an
!   interval in which Q>>0.
!
!  Input parameters:
!   
!    ifleft - an integer parameter indicating whether the value should be computed
!      at the left or right endpoints of the interval;
!
!      ifleft = 1   means compute the value of the nonoscillatory solution at a
!      ifleft = 0   means compute the value of the nonoscillatory solution at b
!
!    eps - the precision used by the adaptive discretization procedure
!    (a,b) - the interval over which the problem is given
!    chebcoefs - structure describing the piecewise discretization scheme used to
!      represent the coefficients
!    qs - the vectors of values of the coefficients
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 4    means that the maximum number of intervals was exceeded
!                  while trying to construct the solution
!       ier = 8    means the internal stack overflowed
!       ier = 16   means that the adaptive subdivision procedure produced an
!                  interval which was too small
!
!  apval, appval - the values of the first and second derivatives of the
!    nonoscillatory phase function at the left or right endpoint of the
!    interval
!

double precision, allocatable            :: ts(:)
double precision, allocatable            :: stack(:,:), phi(:)

double precision, allocatable            :: qs0(:)
double precision, allocatable            :: qs1(:)
double precision, allocatable            :: qsw(:)
type(chebpw_scheme)                      :: chebout

double complex, allocatable              :: rs0(:), rders0(:)
double complex                           :: ima

double complex                           :: rval, rder

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
maxiters = vars%maxiters
maxsteps = vars%maxsteps


ier      = 0
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
maxstack = maxints
nints0   = 0
ind      = 0

!
!  Determine parameters for the window function
!

if (ifleft .eq. 1) then
x = b
call chebpw_interp(chebcoefs,qs,x,dnu)
else
x = a
call chebpw_interp(chebcoefs,qs,x,dnu)
endif

! ! nn = 30
! ! allocate(ts(nn),qsw(nn))
! ! do i=1,nn
! ! ts(i) = a + (i-1.0d0)/(nn-1.0d0)*(b-a)
! ! end do
! ! call chebpw_interp(chebcoefs,qs,nn,ts,qsw)
! ! dnu  = maxval(qsw)
! ! deallocate(qsw,ts)

c1   = (b+a)/2
c2   = 11.0d0/(b-a)
if (eps0 .lt. 1.0d-16) c2   = 12.0d0/(b-a)
if (eps0 .lt. 1.0d-30) c2   = 16.0d0/(b-a)

dnu  = sqrt(dnu)
rval = ima*dnu


allocate(stack(2,maxstack))
allocate(rs0(kcheb),phi(kcheb),ts(kcheb), rders0(kcheb))
allocate(qs0(kcheb))
allocate(qs1(kcheb))
allocate(qsw(kcheb))

nints0     = 0
nstack     = 1
stack(1,1) = a
stack(2,1) = b


if (ifleft .eq. 0) then

do while (nstack > 0 )
a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

if (b0-a0 .lt. epsint) then
ier = 16
return
endif


ifaccept = 0
dq       = 1d300
dcoefs   = 1d300
jer      = -1

ts = (b0-a0)/2 * vars%xscheb + (b0+a0)/2
call chebpw_interp(chebcoefs,qs,kcheb,ts,qs1)

phi = (erf(c2*(ts-c1))+1)/2.0d0
qs0  = dnu**2
qsw  = (1-phi)*qs0 + phi*qs1

call riccati_nop_fitr(vars,qsw,dq)

if (dq .lt. eps) then

   call riccati_nop_ivp_specint_c(vars,jer,eps,a0,b0,kcheb,ts,qsw,rval,rs0,rders0)
   if (jer .eq. 0 ) then
   call riccati_nop_fitc(vars,rs0,dcoefs)
   if (dcoefs .lt. eps) ifaccept = 1
   endif

endif

! ind = ind+1
! write (*,"(2(I6,1X),4(D20.10,2X),1X,I2)")   ind,jer,a0,b0,dq,dcoefs,ifaccept
! write (13,"(2(I6,1X),4(D20.10,2X),1X,I2)")  ind,jer,a0,b0,dq,dcoefs,ifaccept

if (ifaccept .eq. 0) then

if (nstack+2 .ge. maxstack) then
ier = 8
return
endif

nstack = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

rval   = rs0(kcheb)
rder   = rders0(kcheb)
nints0 = nints0+1
if (nints0 .gt. maxints) then
ier = 4
return
endif

endif

end do

else


do while (nstack > 0 )
a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

if (b0-a0 .lt. epsint) then
ier = 16
return
endif

ifaccept = 0
dq       = 1d300
dcoefs   = 1d300
jer      = -1

ts = (b0-a0)/2 * vars%xscheb + (b0+a0)/2
call chebpw_interp(chebcoefs,qs,kcheb,ts,qs1)
phi = (erf(c2*(ts-c1))+1)/2.0d0
qs0  = dnu**2
qsw  = (phi)*qs0 + (1-phi)*qs1

call riccati_nop_fitr(vars,qsw,dq)

if (dq .le. eps) then

  call riccati_nop_tvp_specint_c(vars,jer,eps,a0,b0,kcheb,ts,qsw,rval,rs0,rders0)
  if (jer .eq. 0 ) then
  call riccati_nop_fitc(vars,rs0,dcoefs)
  if (dcoefs .lt. eps) ifaccept = 1
  endif

endif

! ind = ind+1
! write  (*,"(2(I6,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept
! write  (13,"(2(I6,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept

if (ifaccept .eq. 0) then

if (nstack+2 .ge. maxstack) then
ier = 8
return
endif

nstack = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

nstack = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

else

rval   = rs0(1)
rder   = rders0(1)
nints0 = nints0+1
if (nints0 .gt. maxints) then
ier = 4
return
endif

endif

end do

endif

apval  = imag(rval)
appval = imag(rder)

end subroutine




subroutine riccati_nop_window2(vars,ier,ifleft,eps,a,b,fun,par1,par2,par3,apval,appval)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)                           :: vars
procedure(riccati_nop_fun1)                        :: fun
!
!   Use the "windowing algorithm" to compute the value of the first two derivatives
!   of the nonoscillatory phase function alpha(t) at the left or right endpoint of an
!   interval in which Q>>0.
!
!  Input parameters:
!   
!    ifleft - an integer parameter indicating whether the value should be computed
!      at the left or right endpoints of the interval;
!
!      ifleft = 1   means compute the value of the nonoscillatory solution at a
!      ifleft = 0   means compute the value of the nonoscillatory solution at b
!
!    eps - the precision used by the adaptive discretization procedure
!    (a,b) - the interval over which the problem is given
!    fun - an external subroutine supplying the values of q(t)
!    par? - parameters which are passed to the subroutine fun
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 4    means that the maximum number of intervals was exceeded
!                  while trying to construct the solution
!       ier = 8    means the internal stack overflowed
!       ier = 16   means that the adaptive subdivision procedure produced an
!                  interval which was too small
!
!  apval, appval - the values of the first and second derivatives of the
!    nonoscillatory phase function at the left or right endpoint of the
!    interval
!

double precision, allocatable            :: ts(:)
double precision, allocatable            :: stack(:,:), phi(:)

double precision, allocatable            :: qs0(:)
double precision, allocatable            :: qs1(:)
double precision, allocatable            :: qsw(:)
type(chebpw_scheme)                      :: chebout

double complex, allocatable              :: rs0(:), rders0(:)
double complex                           :: ima

double complex                           :: rval, rder

double precision :: ts00(1), qs00(1)

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
maxiters = vars%maxiters
maxsteps = vars%maxsteps


ier      = 0
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
maxstack = maxints
nints0   = 0
ind      = 0

!
!  Determine parameters for the window function
!

if (ifleft .eq. 1) then
ts00(1) = b
call fun(1,ts00,qs00,par1,par2,par3)
dnu = qs00(1)
else
x = a
ts00(1) = a
call fun(1,ts00,qs00,par1,par2,par3)
dnu = qs00(1)
endif

! ! nn = 30
! ! allocate(ts(nn),qsw(nn))
! ! do i=1,nn
! ! ts(i) = a + (i-1.0d0)/(nn-1.0d0)*(b-a)
! ! end do
! ! call chebpw_interp(chebcoefs,qs,nn,ts,qsw)
! ! dnu  = maxval(qsw)
! ! deallocate(qsw,ts)

c1   = (b+a)/2
c2   = 11.0d0/(b-a)
if (eps0 .lt. 1.0d-16) c2   = 12.0d0/(b-a)
if (eps0 .lt. 1.0d-30) c2   = 16.0d0/(b-a)

dnu  = sqrt(dnu)
rval = ima*dnu


allocate(stack(2,maxstack))
allocate(rs0(kcheb),phi(kcheb),ts(kcheb), rders0(kcheb))
allocate(qs0(kcheb))
allocate(qs1(kcheb))
allocate(qsw(kcheb))


nints0     = 0
nstack     = 1
stack(1,1) = a
stack(2,1) = b


if (ifleft .eq. 0) then

do while (nstack > 0 )
a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

if (b0-a0 .lt. epsint) then
ier = 16
return
endif


ifaccept = 0
dq       = 1d300
dcoefs   = 1d300
jer      = -1

ts = (b0-a0)/2 * vars%xscheb + (b0+a0)/2
!call chebpw_interp(chebcoefs,qs,kcheb,ts,qs1)
call fun(kcheb,ts,qs1,par1,par2,par3)

phi = (erf(c2*(ts-c1))+1)/2.0d0
qs0  = dnu**2
qsw  = (1-phi)*qs0 + phi*qs1

call riccati_nop_fitr(vars,qsw,dq)

if (dq .lt. eps) then

   call riccati_nop_ivp_specint_c(vars,jer,eps,a0,b0,kcheb,ts,qsw,rval,rs0,rders0)
   if (jer .eq. 0 ) then
   call riccati_nop_fitc(vars,rs0,dcoefs)
   if (dcoefs .lt. eps) ifaccept = 1
   endif

endif

! ind = ind+1
! write (*,"(2(I6,1X),4(D20.10,2X),1X,I2)")   ind,jer,a0,b0,dq,dcoefs,ifaccept
! write (13,"(2(I6,1X),4(D20.10,2X),1X,I2)")  ind,jer,a0,b0,dq,dcoefs,ifaccept

if (ifaccept .eq. 0) then

if (nstack+2 .ge. maxstack) then
ier = 8
return
endif

nstack = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

rval   = rs0(kcheb)
rder   = rders0(kcheb)
nints0 = nints0+1
if (nints0 .gt. maxints) then
ier = 4
return
endif

endif

end do

else


do while (nstack > 0 )
a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

if (b0-a0 .lt. epsint) then
ier = 16
return
endif

ifaccept = 0
dq       = 1d300
dcoefs   = 1d300
jer      = -1

ts = (b0-a0)/2 * vars%xscheb + (b0+a0)/2
!call chebpw_interp(chebcoefs,qs,kcheb,ts,qs1)
call fun(kcheb,ts,qs1,par1,par2,par3)

phi = (erf(c2*(ts-c1))+1)/2.0d0
qs0  = dnu**2
qsw  = (phi)*qs0 + (1-phi)*qs1

call riccati_nop_fitr(vars,qsw,dq)

if (dq .le. eps) then

  call riccati_nop_tvp_specint_c(vars,jer,eps,a0,b0,kcheb,ts,qsw,rval,rs0,rders0)
  if (jer .eq. 0 ) then
  call riccati_nop_fitc(vars,rs0,dcoefs)
  if (dcoefs .lt. eps) ifaccept = 1
  endif

endif

! ind = ind+1
! write  (*,"(2(I6,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept
! write  (13,"(2(I6,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept

if (ifaccept .eq. 0) then

if (nstack+2 .ge. maxstack) then
ier = 8
return
endif

nstack = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

nstack = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

else

rval   = rs0(1)
rder   = rders0(1)
nints0 = nints0+1
if (nints0 .gt. maxints) then
ier = 4
return
endif

endif

end do

endif

apval  = imag(rval)
appval = imag(rder)

end subroutine




subroutine riccati_nop_phase1(vars,ier,eps,a,b,c,apval,appval,d,aval,chebcoefs,qs, &
  chebphase,alpha,alphader,alphader2)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)                           :: vars
type(chebpw_scheme)                                :: chebcoefs
double precision                                   :: qs(:)
type(chebpw_scheme), intent(out)                   :: chebphase
double precision, allocatable, intent(out)         :: alpha(:), alphader(:), alphader2(:)
!
!  Given the values of the first two derivatives of the phase function alpha(t) at
!  a point c in the interval [a,b] and the value of alpha(t) at the point d in [a,b],
!  construct piecewise Chebyshev expansions of alpha(t) and its first two derivatives
!  on the interval [a,b] by solving Riccati's equation.
!
!  The coefficient is specified via a piecewise Chebyshev expansion.
!
!  Input parameters:
!    eps - the desired precision for the solution
!    c - the point at which the values of alpha'(t) and alpha''(t) are given
!    apval - the value of alpha'(c)
!    appval - the value of alpha''(c)
!    chebcoefs - the piecewise Chebyshev discretization scheme used to
!      represent the coefficients
!    qs the vector of values of values of q(t) at the discretization nodes
!
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 4    indicates that the maximum number of intervals was exceeded
!                  while trying to construct the solution
!       ier = 8    means that this routine's internal stack overflowed
!       ier = 16   means that the adaptive subdivision procedure produced an
!                  interval which was too small
!
!    chebphase - a structure describing the piecewise Chebyshev discretization scheme
!      used to represent the solution
!    alpha, alphader, alphader2 - the values of alpha(t) and its first two
!      derivatives at the discretization nodes
!

double complex, allocatable        :: rs(:), rders(:)
double complex                     :: rval, ima

ier = 0

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
maxiters = vars%maxiters
maxsteps = vars%maxsteps


ima = (0.0d0,1.0d0)

!
!  Call the solver
!


rval = -0.5d0*appval/apval + ima*apval
call riccati_nop_solve(vars,ier,eps,a,b,c,rval,chebcoefs,qs,chebphase,rs,rders)
if (ier .ne. 0) return

call chebpw_info(chebphase,k0,nints0)
ncoefs = k0*nints0

allocate(alphader(ncoefs), alphader2(ncoefs))

alphader  = imag(rs)
alphader2 = imag(rders)

!
!  Integrate alpha'(t)
!

call chebpw_int(chebphase,alphader,d,aval,alpha)

end subroutine


subroutine riccati_nop_phase2(vars,ier,eps,a,b,c,apval,appval,d,aval,fun,par1,par2,par3, &
  chebphase,alpha,alphader,alphader2)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)                           :: vars
procedure(riccati_nop_fun1)                        :: fun
type(chebpw_scheme), intent(out)                   :: chebphase
double precision, allocatable, intent(out)         :: alpha(:), alphader(:), alphader2(:)
!
!  Given the values of the first two derivatives of the phase function alpha(t) at
!  a point c in the interval [a,b] and the value of alpha(t) at the point d in [a,b],
!  construct piecewise Chebyshev expansions of alpha(t) and its first two derivatives
!  on the interval [a,b] by solving Riccati's equation.
!
!  The coefficient is specified via an external subroutine.
!
!  Input parameters:
!    eps - the desired precision for the solution
!    c - the point at which the values of alpha'(t) and alpha''(t) are given
!    apval - the value of alpha'(c)
!    appval - the value of alpha''(c)
!    fun - an external subroutine supplying the values of q(t)
!    par? - parameters which are passed to the subroutine fun
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 4    indicates that the maximum number of intervals was exceeded
!                  while trying to construct the solution
!       ier = 8    means that this routine's internal stack overflowed
!       ier = 16   means that the adaptive subdivision procedure produced an
!                  interval which was too small
!
!    chebphase - a structure describing the piecewise Chebyshev discretization scheme
!      used to represent the solution
!    alpha, alphader, alphader2 - the values of alpha(t) and its first two
!      derivatives at the discretization nodes
!

double complex, allocatable        :: rs(:), rders(:)
double complex                     :: rval, ima

ier = 0

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
maxiters = vars%maxiters
maxsteps = vars%maxsteps


ima = (0.0d0,1.0d0)

!
!  Call the solver
!


rval = -0.5d0*appval/apval + ima*apval
call riccati_nop_solve(vars,ier,eps,a,b,c,rval,fun,par1,par2,par3,chebphase,rs,rders)
if (ier .ne. 0) return

call chebpw_info(chebphase,k0,nints0)
ncoefs = k0*nints0

allocate(alphader(ncoefs), alphader2(ncoefs))

alphader  = imag(rs)
alphader2 = imag(rders)

!
!  Integrate alpha'(t)
!

call chebpw_int(chebphase,alphader,d,aval,alpha)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Solver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine riccati_nop_solve1(vars,ier,eps,a,b,c,rval,chebcoefs,qs,chebout,rs,rders)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)                           :: vars
type(chebpw_scheme)                                :: chebcoefs
double precision                                   :: qs(:)
type(chebpw_scheme), intent(out)                   :: chebout
double complex, allocatable, intent(out)           :: rs(:), rders(:)
double complex                                     :: rval
!
!  Adaptively construct a piecewise Chebyshev expansion which represents the solution
!  of the problem
!
!    { r'(t) + (r(t))^2 + p(t) r(t) + q(t) = 0      a <= t <= b
!    { r(c)                                = rval.
!
!  The coefficient are specified as piecewise Chebyshev expansions.
!
!  Input parameters:
!    eps - the desired precision for the solution
!    (a,b) - the interval over which the problem is given
!    c - the point at which the value of the solution is specified
!    rval - the value of the solution
!    chebcoefs - structure describing the piecewise discretization scheme used to
!      represent the coefficients
!    qs - the vectors of values of the coefficient
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 4    means that the maximum number of intervals was exceeded
!                  while trying to construct the solution
!       ier = 8    means the internal stack overflowed
!       ier = 16   means that the adaptive subdivision procedure produced an
!                  interval which was too small
!
!    chebout - the piecewise Chebyshev discretization scheme used to represent the
!      solution
!    rs - the values of the solution at the piecewise Chebyshev discretization nodes
!    rders - the values of the derivative of the  solution at the piecewise Chebyshev
!       discretization nodes
!

double precision, allocatable            :: stack(:,:), ab(:,:)
double precision, allocatable            :: qs0(:), ts0(:)

double complex, allocatable              :: rs1(:,:), rders1(:,:)
double precision, allocatable            :: ab1(:,:)

double complex, allocatable              :: rs2(:,:), rders2(:,:)
double precision, allocatable            :: ab2(:,:)

double complex                           :: rval0

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
maxiters = vars%maxiters
maxsteps = vars%maxsteps

ier      = 0
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
maxstack = maxints
nints1   = 0
nints2   = 0


allocate(stack(2,maxstack))
allocate(ts0(kcheb),qs0(kcheb))

allocate(ab1(2,maxints))
allocate(rs1(kcheb,maxints+1))
allocate(rders1(kcheb,maxints+1))

allocate(ab2(2,maxints))
allocate(rs2(kcheb,maxints+1))
allocate(rders2(kcheb,maxints+1))

!
!  Solve backward on the interval (a,c)
!


if (a .lt. c) then

nstack          = 1
stack(1,nstack) = a
stack(2,nstack) = c
ind             = 0

do while (nstack > 0) 
a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

if (b0-a0 .lt. epsint) then
ier = 16
return
endif



if (nints1 .eq. 0) then
rval0 = rval
else
rval0 = rs1(1,nints1)
endif

ifaccept = 0
dcoefs   = 1d300
dq       = 13d00
jer      = -1

ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
call chebpw_interp(chebcoefs,qs,kcheb,ts0,qs0)

call riccati_nop_fitr(vars,qs0,dq)

if (dq .lt. eps) then

  call riccati_nop_tvp_specint_c(vars,jer,eps,a0,b0,kcheb,ts0,qs0,rval0, &
    rs1(:,nints1+1),rders1(:,nints1+1))

  if (jer .eq. 0) then
  call riccati_nop_fitc(vars,rs1(:,nints1+1),dcoefs)
  if (dcoefs .lt. eps) ifaccept=1
  endif

endif

! ind = ind+1
! write (*,"(2(I5,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept
! write (13,"(2(I5,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept

if (ifaccept .eq. 0) then

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif

nstack = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

nstack = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

else

if (nints1+1 .gt. maxints) then
ier = 4
return
endif

nints1        = nints1+1
ab1(1,nints1) = a0
ab1(2,nints1) = b0

endif

end do

endif


!
!  Solve forward on the interval (c,b)
!


if (c .lt. b) then


ind    = 0
nstack = 1
stack(1,nstack) = c
stack(2,nstack) = b

do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1


if (b0-a0 .lt. epsint) then
ier = 16
return
endif

if (nints2 .eq. 0) then
rval0 = rval
else
rval0 = rs2(kcheb,nints2)
endif

ifaccept = 0
dq       = 1d300
dcoefs   = 1d300
jer      = -1

ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
call chebpw_interp(chebcoefs,qs,kcheb,ts0,qs0)

call riccati_nop_fitr(vars,qs0,dq)

if (dq .lt. eps) then

  call riccati_nop_ivp_specint_c(vars,jer,eps,a0,b0,kcheb,ts0,qs0,rval0,&
    rs2(:,nints2+1),rders2(:,nints2+1))
  if (jer .eq. 0) then
  call riccati_nop_fitc(vars,rs2(:,nints2+1),dcoefs)
  if (dcoefs .lt. eps) ifaccept = 1
  endif

endif

! ind = ind+1
! write (*,"(2(I5,1X),4(D20.10,2X),I3)")   ind,jer,a0,b0,dq,dcoefs,ifaccept
! write (13,"(2(I5,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept

if (ifaccept .eq. 0) then


nstack = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2


else

if (nints2+1 .ge. maxints) then
ier = 4
return
endif

nints2         = nints2+1
ab2(1,nints2) = a0
ab2(2,nints2) = b0
a0            = b0

endif

end do

endif



!
!  Copy out the solution
!


nints = nints1+nints2

allocate(ab(2,1:nints))
ab(:,1:nints1)       = ab1(:,nints1:1:-1)
ab(:,nints1+1:nints) = ab2(:,1:nints2)

call chebpw_specified(chebout,kcheb,nints,ab)

! nints          = nints1+nints2
! chebout%a      = a
! chebout%b      = b
! chebout%k      = kcheb
! chebout%nints  = nints

! allocate(chebout%xscheb(kcheb), chebout%whtscheb(kcheb) )
! allocate(chebout%acoefs(kcheb,kcheb))
! allocate(chebout%aintl(kcheb,kcheb))
! allocate(chebout%aintr(kcheb,kcheb))

! allocate(chebout%ab(2,nints))

! chebout%xscheb                        = xscheb
! chebout%whtscheb                      = whtscheb
! chebout%acoefs                        = acoefs
! chebout%aintl                         = aintl
! chebout%aintr                         = aintr
! chebout%ab(:,1:nints1)                = ab1(:,nints1:1:-1)
! chebout%ab(:,nints1+1:nints)          = ab2(:,1:nints2)
! chebout%nints                         = nints

allocate(rs(kcheb*nints), rders(kcheb*nints) )

n1 = kcheb*nints1
n2 = kcheb*nints

rs(1:n1)     = reshape(rs1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
rs(n1+1:n2)  = reshape(rs2(1:kcheb,1:nints2), [kcheb*nints2] )

rders(1:n1)     = reshape(rders1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
rders(n1+1:n2)  = reshape(rders2(1:kcheb,1:nints2), [kcheb*nints2] )

end subroutine


subroutine riccati_nop_solve2(vars,ier,eps,a,b,c,rval,fun,par1,par2,par3,chebout,rs,rders)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)                           :: vars
procedure(riccati_nop_fun1)                        :: fun
type(chebpw_scheme), intent(out)                   :: chebout
double complex, allocatable, intent(out)           :: rs(:), rders(:)
double complex                                     :: rval
!
!  Adaptively construct a piecewise Chebyshev expansion which represents the solution
!  of the problem
!
!    { r'(t) + (r(t))^2 + p(t) r(t) + q(t) = 0      a <= t <= b
!    { r(c)                                = rval.
!
!  The coefficient is specified via an external subroutine.
!
!  Input parameters:
!    eps - the desired precision for the solution
!    (a,b) - the interval over which the problem is given
!    c - the point at which the value of the solution is specified
!    rval - the value of the solution
!    fun  - an exteranl subroutine which supplies the values of the coefficient q(t)
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 4    means that the maximum number of intervals was exceeded
!                  while trying to construct the solution
!       ier = 8    means the internal stack overflowed
!       ier = 16   means that the adaptive subdivision procedure produced an
!                  interval which was too small
!
!    chebout - the piecewise Chebyshev discretization scheme used to represent the
!      solution
!    rs - the values of the solution at the piecewise Chebyshev discretization nodes
!    rders - the values of the derivative of the  solution at the piecewise Chebyshev
!       discretization nodes
!

double precision, allocatable            :: stack(:,:), ab(:,:)
double precision, allocatable            :: qs0(:), ts0(:)

double complex, allocatable              :: rs1(:,:), rders1(:,:)
double precision, allocatable            :: ab1(:,:)

double complex, allocatable              :: rs2(:,:), rders2(:,:)
double precision, allocatable            :: ab2(:,:)

double complex                           :: rval0

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
maxiters = vars%maxiters
maxsteps = vars%maxsteps

ier      = 0
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
maxstack = maxints
nints1   = 0
nints2   = 0


allocate(stack(2,maxstack))
allocate(ts0(kcheb),qs0(kcheb))

allocate(ab1(2,maxints))
allocate(rs1(kcheb,maxints+1))
allocate(rders1(kcheb,maxints+1))

allocate(ab2(2,maxints))
allocate(rs2(kcheb,maxints+1))
allocate(rders2(kcheb,maxints+1))

!
!  Solve backward on the interval (a,c)
!


if (a .lt. c) then

nstack          = 1
stack(1,nstack) = a
stack(2,nstack) = c
ind             = 0

do while (nstack > 0) 
a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

if (b0-a0 .lt. epsint) then
ier = 16
return
endif



if (nints1 .eq. 0) then
rval0 = rval
else
rval0 = rs1(1,nints1)
endif


ifaccept = 0
dcoefs   = 1d300
dq       = 13d00
jer      = -1


ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
call fun(kcheb,ts0,qs0,par1,par2,par3)
!call chebpw_interp(chebcoefs,qs,kcheb,ts0,qs0)

call riccati_nop_fitr(vars,qs0,dq)

if (dq .lt. eps) then

  call riccati_nop_tvp_specint_c(vars,jer,eps,a0,b0,kcheb,ts0,qs0,rval0, &
    rs1(:,nints1+1),rders1(:,nints1+1))

  if (jer .eq. 0) then
  call riccati_nop_fitc(vars,rs1(:,nints1+1),dcoefs)
  if (dcoefs .lt. eps) ifaccept=1
  endif

endif

! ind = ind+1
! write (*,"(2(I5,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept
! write (13,"(2(I5,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept

if (ifaccept .eq. 0) then

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif

nstack = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

nstack = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

else

if (nints1+1 .gt. maxints) then
ier = 4
return
endif

nints1        = nints1+1
ab1(1,nints1) = a0
ab1(2,nints1) = b0

endif

end do

endif


!
!  Solve forward on the interval (c,b)
!


if (c .lt. b) then


ind    = 0
nstack = 1
stack(1,nstack) = c
stack(2,nstack) = b

do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1


if (b0-a0 .lt. epsint) then
ier = 16
return
endif

if (nints2 .eq. 0) then
rval0 = rval
else
rval0 = rs2(kcheb,nints2)
endif

ifaccept = 0
dq       = 1d300
dcoefs   = 1d300
jer      = -1

ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
call fun(kcheb,ts0,qs0,par1,par2,par3)

!call chebpw_interp(chebcoefs,qs,kcheb,ts0,qs0)

call riccati_nop_fitr(vars,qs0,dq)

if (dq .lt. eps) then

  call riccati_nop_ivp_specint_c(vars,jer,eps,a0,b0,kcheb,ts0,qs0,rval0,&
    rs2(:,nints2+1),rders2(:,nints2+1))
  if (jer .eq. 0) then
  call riccati_nop_fitc(vars,rs2(:,nints2+1),dcoefs)
  if (dcoefs .lt. eps) ifaccept = 1
  endif

endif

! ind = ind+1
! write (*,"(2(I5,1X),4(D20.10,2X),I3)")   ind,jer,a0,b0,dq,dcoefs,ifaccept
! write (13,"(2(I5,1X),4(D20.10,2X),I3)")  ind,jer,a0,b0,dq,dcoefs,ifaccept

if (ifaccept .eq. 0) then


nstack = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2


else

if (nints2+1 .ge. maxints) then
ier = 4
return
endif

nints2         = nints2+1
ab2(1,nints2) = a0
ab2(2,nints2) = b0
a0            = b0

endif

end do

endif



!
!  Copy out the solution
!


nints = nints1+nints2

allocate(ab(2,1:nints))
ab(:,1:nints1)       = ab1(:,nints1:1:-1)
ab(:,nints1+1:nints) = ab2(:,1:nints2)

call chebpw_specified(chebout,kcheb,nints,ab)

! nints          = nints1+nints2
! chebout%a      = a
! chebout%b      = b
! chebout%k      = kcheb
! chebout%nints  = nints

! allocate(chebout%xscheb(kcheb), chebout%whtscheb(kcheb) )
! allocate(chebout%acoefs(kcheb,kcheb))
! allocate(chebout%aintl(kcheb,kcheb))
! allocate(chebout%aintr(kcheb,kcheb))

! allocate(chebout%ab(2,nints))

! chebout%xscheb                        = xscheb
! chebout%whtscheb                      = whtscheb
! chebout%acoefs                        = acoefs
! chebout%aintl                         = aintl
! chebout%aintr                         = aintr
! chebout%ab(:,1:nints1)                = ab1(:,nints1:1:-1)
! chebout%ab(:,nints1+1:nints)          = ab2(:,1:nints2)
! chebout%nints                         = nints

allocate(rs(kcheb*nints), rders(kcheb*nints) )

n1 = kcheb*nints1
n2 = kcheb*nints

rs(1:n1)     = reshape(rs1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
rs(n1+1:n2)  = reshape(rs2(1:kcheb,1:nints2), [kcheb*nints2] )

rders(1:n1)     = reshape(rders1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
rders(n1+1:n2)  = reshape(rders2(1:kcheb,1:nints2), [kcheb*nints2] )

end subroutine


subroutine riccati_nop_tvp_specint_c(vars,ier,eps,a,b,kcheb,ts,qs,rb,rs,rders)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)                           :: vars
double precision                                   :: ts(:), qs(:)
double complex                                     :: rs(:), rders(:)
double complex                                     :: rb
!
!  Calculate the solution of the terminal value problem
!
!    { r'(t) + (r(t))^2 + p(t) r(t) + q(t) = 0          a < t < b
!    { r(b)                                = alpha
!
!  at the Chebyshev nodes on the interval [a,b].  An initial approximation is
!  formed using the trapezoidal method and then a spectral integration
!  method is used to iteratively refine the solution.
!
!  Input parameters:
!    eps - the desired precision for the calculations
!    (a,b) - the interval on which to solve the Riccati equation
!    ts - the Chebyshev nodes on the the interval (a,b)
!    qs - the values of the coefficient q(t) at the Chebyshev nodes
!    rb - the desired initial value for the solution
!
!  Output parameters: 
!    ier - an error return code;
!      ier = 0     indicates successful execution
!      ier = 1024  means that NaN or Inf was encountered
!      ier = 2048   means that Newton's method failed to converge
!
!    rs - the values of the computed solution
!    rders - the values of the derivative of the computed solution
!
!
double complex             :: ima
double complex             :: res(kcheb)
double complex             :: amatr(kcheb,kcheb)
double complex             :: rs2(kcheb)
double complex             :: rders2(kcheb)
double complex             :: res2(kcheb)
double complex             :: rb0

ier      = 0
nsteps   = 0 
ima      = (0.0d0,1.0d0)

kcheb    = vars%kcheb
maxiters = vars%maxiters
maxsteps = vars%maxsteps

!
! Form an initial approximation 
!
rs(kcheb) = rb
call riccati_nop_tvp_trap_c(ier,maxiters,eps,kcheb,ts,qs,rs)
if (ier .ne. 0 ) return

!
!  d'(t) + 2 r(t) d(t) + p(t) d(t) = -res(t)
!
!                  t 
!  d(t) = rb + int    d'(s) ds
!                  b
!
!                                  t
!  d'(t) + ( 2 r(t) + p(t) )   \int d'(s) ds   = -res(t) - rb ( 2 r(t) + p )
!                                  b 
!

rders  = -rs**2  - qs
rs     = rb + (b-a)/2*matmul(vars%aintr,rders)
dnorm0 = norm2(abs(rs))
if (dnorm0 .eq. 0) dnorm0 = 1

!
! Solve the linearized system
!

res   = rders + rs**2  + qs

do istep=1,maxsteps

rders2 = -res
rb0    = rb - rs(kcheb)


do i=1,kcheb
amatr(i,:) = (2*rs(i))*vars%aintr(i,:)*(b-a)/2
amatr(i,i) = amatr(i,i) + 1.0d0
rders2(i)  = rders2(i) - rb0*(2*rs(i))
end do

call linalg0_solve(kcheb,amatr,rders2)

do i=1,kcheb
ddd = imag(rders2(i))
if (ieee_is_nan(ddd) .or. .not. ieee_is_finite(ddd))  then
ier = 1024
return
endif
end do

rs2    = rb0 + (b-a)/2*matmul(vars%aintr,rders2)
dnorm  = norm2(abs(rs2)) 
rs2    = rs2 + rs

rders2 = rders2 + rders
res2   = rders2 + rs2**2  + qs

rs     = rs2
rders  = rders2
res    = res2
dres1  = dres2
nsteps = istep

if (dnorm .lt. eps*dnorm0) exit

end do

if (istep .gt. maxsteps) then
ier = 2048
return
endif


end subroutine


subroutine riccati_nop_tvp_trap_c(ier,maxiters,eps,k,ts,qs,rs)
implicit double precision (a-h,o-z)
integer                                            :: k
double precision                                   :: ts(k), qs(k)
double complex                                     :: rs(k), rb
!
!  Use the implicit trapezoidal method to solve the terminal value problem
! 
!     r'(t) + (r(t))^2 + q(t) = 0      a < t < b
!     r(b)                    = rb
!
!  The solution is tabulated at a collection of nodes in the interval (a,b)
!  specified by the user the right-most of which must be the right  endpoint
!  b.
!
!  Input parameters:
!    k - the number of nodes at which to approximate the solution
!    ts - an array containing the nodes
!    qs - the values of the function q at the nodes
!    rb - the terminal value for the solution
!
!  Output parameters:
!    ier - an error return code;
!      ier =    0      indicates successful execution
!      ier =    1024   indicates that NaN or Inf was encountered
!
!    rs - the values of the solution at the specified nodes
!    rders - the values of the derivative of the solution at the specified
!      nodes
!    

double complex               :: dd1, dd2, r0, r1, r
double complex               :: val, der, delta

ier      = 0

do j=k-1,1,-1

t0 = ts(j)
t1 = ts(j+1)
h  = t1-t0

q0 = qs(j)
q1 = qs(j+1)

r1 = rs(j+1)

r    = r1
dd2  = -r1-h/2*(q0+q1)-h/2*r1**2

do iter=1,maxiters
val   = -h/2*r**2+r+dd2
der   = -h*r+1
delta = -val/der
r     = r+delta
if (abs(delta) .lt. eps*abs(r)) exit
end do

ddd = imag(r)
if (ieee_is_nan(ddd) .OR. .not. ieee_is_finite(ddd)) then
ier = 1024
return
endif

rs(j) = r

end do


end subroutine


subroutine riccati_nop_ivp_specint_c(vars,ier,eps,a,b,kcheb,ts,qs,ra,rs,rders)
implicit double precision (a-h,o-z)
type(riccati_nop_vars_t)   :: vars
double precision           :: ts(:), qs(:)
double complex             :: rs(:), rders(:)
double complex             :: ra
!
!  Calculate the solution of the initial value problem
!
!    { r'(t) + (r(t))^2 + p(t) r(t) + q(t) = 0          a < t < b
!    { r(a)                                = alpha
!
!  at the Chebyshev nodes on the interval [a,b].  An initial approximation is
!  formed using the trapezoidal method and then a spectral integration
!  method is used to iteratively refine the solution.
!
!  Input parameters:
!    eps - the desired precision for the calculations
!    (a,b) - the interval on which to solve the Riccati equation
!    ts - the Chebyshev nodes on the the interval (a,b)
!    qs - the values of the coefficient q(t) at the Chebyshev nodes
!    ra - the desired initial value for the solution
!
!  Output parameters: 
!    ier - an error return code;
!      ier = 0     indicates successful execution
!      ier = 1024  means that NaN or Inf was encountered
!      ier = 2048  means that Newton's method failed to converge
!
!    rs - the values of the computed solution
!    rders - the values of derivative of the computed solution

!
!
double complex             :: ima
double complex             :: res(kcheb)
double complex             :: amatr(kcheb,kcheb)

double complex             :: rs2(kcheb)
double complex             :: rders2(kcheb)
double complex             :: res2(kcheb)
double complex             :: ra0

ier      = 0
ima      = (0.0d0,1.0d0)
nsteps   = 0

kcheb    = vars%kcheb
maxiters = vars%maxiters
maxsteps = vars%maxsteps

!
! Form an initial approximation 
!
rs(1) = ra
call riccati_nop_ivp_trap_c(ier,maxiters,eps,kcheb,ts,qs,rs)
if (ier .ne. 0) return


!rs = ima * sqrt(qs)

rders  = -rs**2 - qs
rs     = ra + (b-a)/2*matmul(vars%aintl,rders)
dnorm0 = norm2(abs(rs))
if (dnorm0 .eq. 0) dnorm0 = 1

!
! Solve the system
!
! d'(t) + (2 r(t) + p(t)) d(t) + res(t)    = 0
! d(a)                                     = ra0 = ra-rs(1)
!
! by solving the integral equation obtained by letting
!
!                 t
! d(t) = ra0 + int  sigma(s) ds
!                 a
!
!                          t
! sigma(t) + (2 r + p )\int   sigma(s) ds = -res(t) - 
!                          a
!

res   = rders + rs**2  + qs
!dres1 = maxval(abs(res))

do istep=1,maxsteps

rders2 = -res
ra0    = ra - rs(1)

do i=1,kcheb
amatr(i,:) = (2*rs(i))*vars%aintl(i,:)*(b-a)/2
amatr(i,i) = amatr(i,i) + 1.0d0
rders2(i)  = rders2(i) - ra0*(2*rs(i))
end do

call linalg0_solve(kcheb,amatr,rders2)

do i=1,kcheb
ddd = imag(rders2(i))
if (ieee_is_nan(ddd) .OR. .not. ieee_is_finite(ddd)) then
ier = 1024
return
endif
end do

rs2    = ra0 + (b-a)/2*matmul(vars%aintl,rders2)
dnorm  = norm2(abs(rs2))

rs2    = rs2 + rs


rders2 = rders2 + rders
res2   = rders2 + rs2**2  + qs
!dres2  = maxval(abs(res2))

!if (dres1 .lt. 2*dres2) exit

rs     = rs2
rders  = rders2
res    = res2
dres1  = dres2
nsteps = istep


if (dnorm .lt. eps*dnorm0) exit

end do


if (istep .gt. maxsteps) then
ier = 2048
return
endif

end subroutine


subroutine riccati_nop_ivp_trap_c(ier,maxiters,eps,k,ts,qs,rs)
implicit double precision (a-h,o-z)
double precision            :: ts(:), qs(:)
double complex              :: rs(:)
!
!  Approximate the values of the solution of the IVP
!
!     { r'(t) + (r(t))^2 + p(t) * r(t) + q(t) = 0
!     { r(a)                                  = alpha,
!
!  where p and q are real-valued, at a collection of points using the 
!  implicit trapezoidal method.
!
!  Input parameters:
!    ier - an error return code; 
!      ier = 0      indicates successful execution
!      ier = 1024  indicates that NaN or Inf was encountered
!
!    eps - precision for the calculations
!    k - the number of points at which to compute the solution
!    ts - an array of length k specifying the points
!    qs - the values of the coefficient q at the Chebyshev nodes
!    rs(1) - the initial value alpha for the solution
!
!  Output parameters:
!    vals - the values of the solution of the IVP at the specified nodes
!
!
double complex               :: dd, r0, r
double complex               :: val,der,delta

ier      = 0
niters   = 0

do i=2,k
t0 = ts(i-1)
t  = ts(i)
h  = t-t0

r0 = rs(i-1)
q0 = qs(i-1)
q  = qs(i)

dd = -r0+h/2*r0**2+h/2*q0+h/2*q
r  = r0
!-h*(r0**2-q0)

do iter=1,maxiters
val   = r+h/2*r**2+h/2*p*r+dd
der   = 1+h*r+h/2*p
delta = -val/der
r     = r+delta
if (abs(delta) .lt. eps*abs(r)) exit
end do

ddd = imag(r)
if (ieee_is_nan(ddd) .OR. .not. ieee_is_finite(ddd)) then
ier = 1024
return
endif

rs(i) = r
end do

end subroutine

end module
 
