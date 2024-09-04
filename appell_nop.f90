!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for solving the special form 
!
!    w'''(t) + 4 q(t) w'(t) +  2 q'(t) w(t) = 0,  a <= t <= b,                            (1)
!
!  of Appell's equation which is satisfied by the products of solutions of the
!  normal form
!
!    y''(t) + q(t) y(t) = 0,  a< = t <= b.                                                (2)
!
!  of the Helmholtz equation.  Coefficients are specified either via piecewise
!  Chebyshev expansions or external subroutines.
!
!  The equation (1) is of interest because if u,v is a pair of linearly independent
!  solutions of (2) then
!
!                       1
!    alpha'(t) = -----------------                                                        (3)
!                 u(t)^2 + v(t)^2
!
!  is the derivative of a phase function for (2) and
!
!    w(t) = 1 / alpha'(t)
!
!  solves (1).  That alpha(t) is a phase function for (2) means that
!
!
!               cos ( alpha(t) )                  sin ( alpha(t ) ) 
!    u(t) =     ------------------  and   v(t) =  -----------------                       (4)
!               sqrt( alpha'(t) )                 sqrt( alpha'(t) )
!
!
!  form a basis in the space of solutions of (2).
!
!  The follow subroutines should be regarded as publicly callable:
!
!    appell_nop_init - set parameters for the solver and construct the various spectral
!      integration and differentiation matrices it requires
!
!    appell_nop_solve - solve Appell's equation on an interval [a,b] given its
!      initial values at a point c in [a,b]. 
!
!    appell_nop_phase - given the values of alpha'(t) and alpha''(t) at a point c
!      in an interval [a,b] and the value of alpha(t) at a point d in [a,b],
!      solve Appell's equation in order to construct piecewise Chebyshev expansions
!      of alpha(t) and its first two derivatives.  
!
!    appell_nop_extend - given piecewise Chebyshev expansions of the phase function
!      alpha(t) and its first two derivatives on an interval [a,b], try to extend
!      alpha(t) on one or both sides of the the interval [a,b] an interval
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Observations / To do:
!
!    - The windowing algorithm implemented in riccati.f90 and riccati_nop.f90
!      does not work for Appell's equation.
!
!    - Appell's equation is much less useful when the equation is of the form
!      y''(t) + p(t) y'(t) + q(t) y(t) = 0.  In that event, the function
!      one wants is  alpha'(t)(t) = [u,v][t]/ (u(t)^2 + v(t)^2) with [u,v] the 
!      Wronskian of the pair [u,v].  The solution of Appell's equation is merely
!      u(t)^2 + v(t)^2.  If one solves Appell's equation and constructs alpha'(t) 
!      with the assumption [u,v](c) = 1 at some point c, then the resulting procedure 
!      is quite sensitive to errors in the initial values for Appell's equation.
!
!      All this is moot when the term p(t) y'(t) is missing, because the Wronskian
!      of any pair of solutions is constant in that case.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module      appell_nop

use utils
use linalg0
use chebyshev
use chebpw
use ieee_arithmetic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The structure which stores any data needed by the procedures in this modules and
!  which is populated by the initialization routine.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type       appell_nop_vars_t
double precision                                        :: epsint
double precision                                        :: dmax
                                                        
integer                                                 :: ntest
integer                                                 :: kcheb
integer                                                 :: maxints
                                                        
double precision, allocatable                           :: xscheb(:)
double precision, allocatable                           :: whtscheb(:)
double precision, allocatable                           :: acoefs(:,:)
double precision, allocatable                           :: aintl(:,:)
double precision, allocatable                           :: aintl2(:,:)
double precision, allocatable                           :: aintl3(:,:)
                                                        
double precision, allocatable                           :: adiff(:,:)
double precision, allocatable                           :: aintr(:,:)
double precision, allocatable                           :: aintr2(:,:)
double precision, allocatable                           :: aintr3(:,:)
end type   appell_nop_vars_t


interface
subroutine appell_nop_fun1(n,ts,qs,qders,par1,par2,par3)
implicit double precision (a-h,o-z)
double precision :: ts(:), qs(:), qders(:)
end subroutine
end interface

interface           appell_nop_solve
module procedure    appell_nop_solve1
module procedure    appell_nop_solve2
end interface       appell_nop_solve

interface           appell_nop_phase
module procedure    appell_nop_phase1
module procedure    appell_nop_phase2
end interface       appell_nop_phase

interface           appell_nop_extend
module procedure    appell_nop_extend1
module procedure    appell_nop_extend2
end interface       appell_nop_extend

contains

subroutine appell_nop_init(vars,kcheb0,ntest0,maxints0,epsint0,dmax0)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t), intent(out)        :: vars
integer, optional                           :: kcheb0, ntest0, maxints0
double precision, optional                  :: epsint0, dmax0
!
!  Initialize the solver for the Appell solver.
!
!  If this subroutine is called with no other arguments than vars, then
!  reasonable defaults are chosen.
!
!  Input parameters:
!    kcheb0 - the number of terms in the Chebyshev expansions used to represent solutions
!    ntest0 - the number of trailing Chebyshev coefficients which must be small in order
!      to consider a Chebyshev expansion accurate
!    maxints0 - the maximum number of intervals used to represent solutions
!    epsint0 - the smallest allowable interval size
!    dmax0 - at the present time this parameter is unused
!
!  Output parameters:
!    vars - the structure containing all of the data needed by the solver
!

if (.not. present(kcheb0) )  then
kcheb    = 16
ntest    = 2
maxints  = 10000
epsint   = 1.0d-7
dmax     = 1d200
else
kcheb    = kcheb0
ntest    = ntest0
maxints  = maxints0
epsint   = epsint0
dmax     = dmax0
endif


vars%kcheb   = kcheb
vars%ntest   = ntest
vars%maxints = maxints
vars%epsint  = epsint
vars%dmax    = dmax

call chebyshev_quad(kcheb,vars%xscheb,vars%whtscheb)
call chebyshev_coefsmatrix(kcheb,vars%acoefs)
call chebyshev_intlmatrix(kcheb,vars%aintl)
call chebyshev_intrmatrix(kcheb,vars%aintr)

allocate(vars%aintl2(kcheb,kcheb), vars%aintl3(kcheb,kcheb))
allocate(vars%aintr2(kcheb,kcheb), vars%aintr3(kcheb,kcheb))

vars%aintl2 = matmul(vars%aintl,vars%aintl)
vars%aintl3 = matmul(vars%aintl2,vars%aintl)

vars%aintr2 = matmul(vars%aintr,vars%aintr)
vars%aintr3 = matmul(vars%aintr2,vars%aintr)

end subroutine


subroutine appell_nop_fit(vars,vals,dcoefs)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)   :: vars
double precision          :: vals(:), coefs(vars%kcheb)
kcheb  = vars%kcheb
ntest  = vars%ntest
coefs  = matmul(vars%acoefs,vals)
dd1    = norm2(abs(coefs))+1
dd2    = norm2(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1=1
dcoefs = dd2/dd1
end subroutine


subroutine appell_nop_extend1(vars,ier,eps,amax,bmax,chebcoefs,qs,qders,aout,bout,chebphase,&
  alpha,alphader,alphader2)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)                      :: vars
type(chebpw_scheme)                          :: chebcoefs
double precision                             :: qs(:), qders(:)

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
!  The parameters are specified via piecewise Chebyshev expansions.
!
!  Input parameters:
!    eps - the precision for the calculations
!    [amax,bmax] - the largest possible interval
!    chebcoefs - a piecewise Chebyshev discretization scheme giving the coefficients on an
!      interval containing [amax,bmax]
!   qs, qders -the values of the coefficients at the discretization nodes
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
double precision, allocatable              :: qs0(:), qders0(:)
double precision, allocatable              :: ws0(:), wders0(:), wder2s0(:)

double precision, allocatable              :: ab0(:,:)
double precision, allocatable              :: alphapp0(:,:),alphap0(:,:), alpha0(:,:)

double precision, allocatable              :: abnew(:,:)
double precision, allocatable              :: alphanew(:), alphadernew(:), alphader2new(:)


data pi / 3.14159265358979323846264338327950288d0 /

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
dmax     = vars%dmax

ier      = 0
eps0     = epsilon(0.0d0)
maxstack = maxints*2
dmul     = pi/16
call chebpw_interval(chebphase,a,b)

aout     = a
bout     = b

allocate(stack(2,maxstack))

allocate(ts0(kcheb), qs0(kcheb), qders0(kcheb) )
allocate(ws0(kcheb),wders0(kcheb),wder2s0(kcheb))
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
   call chebpw_interp(chebcoefs,qs,b,qval)

   apppval = ( 4*qval*apval**2-4*apval**4+3*appval**2) / (2*apval)

   w2   = 1/apval
   wp2  = -appval/apval**2
   wpp2 = 2*appval**2/apval**3 - apppval/apval**2

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
   call chebpw_interp(chebcoefs,qs,qders,kcheb,ts0,qs0,qders0)
   ifaccept = 0
   dq       = 1d300
   dcoefs   = 1d300
   jer      = -1

   call appell_nop_fit(vars,qs0,dq)

   if (dq .lt. eps) then

   call appell_ivp_specint_nop(vars,jer,a0,b0,kcheb,ts0,qs0,qders0,w2,wp2,wpp2, &
     ws0, wders0, wder2s0)


   if (jer .eq. 0) then
   alphapp0(:,nints2+1) = -wders0/ws0**2
   alphap0(:,nints2+1)  = 1.0d0/ws0
   alpha0(:,nints2+1)   = aval + (b0-a0)/2*matmul(vars%aintl,alphap0(:,nints2+1))
   
   call appell_nop_fit(vars,alphap0(:,nints2+1),dcoefs)

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
     w2               = ws0(kcheb)
     wp2              = wders0(kcheb)
     wpp2             = wder2s0(kcheb)
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
   call chebpw_interp(chebcoefs,qs,a,qval)

   apppval = ( 4*qval*apval**2-4*apval**4+3*appval**2) / (2*apval)

   w1   = 1/apval
   wp1  = -appval/apval**2
   wpp1 = 2*appval**2/apval**3 - apppval/apval**2


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
    
   ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2)  )
   call chebpw_interp(chebcoefs,qs,qders,kcheb,ts0,qs0,qders0)

   call appell_nop_fit(vars,qs0,dq)

   if (dq .lt. eps) then

   call appell_tvp_specint_nop(vars,jer,a0,b0,kcheb,ts0,qs0,qders0,w1,wp1,wpp1, &
     ws0, wders0, wder2s0)
   if (jer .eq. 0) then
   alphapp0(:,nints1-1) = -wders0/ws0**2
   alphap0(:,nints1-1)  = 1.0d0/ws0
   alpha0(:,nints1-1)   = aval + (b0-a0)/2*matmul(vars%aintr,alphap0(:,nints1-1))


   call appell_nop_fit(vars,alphap0(:,nints1-1),dcoefs)
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

     w1               = ws0(1)
     wp1              = wders0(1)
     wpp1             = wder2s0(1)
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
nints = nintsnew

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


subroutine appell_nop_extend2(vars,ier,eps,amax,bmax,fun,par1,par2,par3,aout,bout,chebphase,&
  alpha,alphader,alphader2)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)                      :: vars
procedure(appell_nop_fun1)                   :: fun
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
!  The parameters are specified via an external subroutine
!
!  Input parameters:
!    eps - the precision for the calculations
!    [amax,bmax] - the largest possible interval
!    fun - the external subroutine supplying the values of q(t) and q'(t)
!    par? - user-supplied parameters which are passed to fun
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
double precision, allocatable              :: qs0(:), qders0(:)
double precision, allocatable              :: ws0(:), wders0(:), wder2s0(:)

double precision, allocatable              :: ab0(:,:)
double precision, allocatable              :: alphapp0(:,:),alphap0(:,:), alpha0(:,:)

double precision, allocatable              :: abnew(:,:)
double precision, allocatable              :: alphanew(:), alphadernew(:), alphader2new(:)

double precision                           :: ts00(1), qs00(1), qders00(1)

data pi / 3.14159265358979323846264338327950288d0 /

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
dmax     = vars%dmax

ier      = 0
eps0     = epsilon(0.0d0)
maxstack = maxints*2
dmul     = pi/16
call chebpw_interval(chebphase,a,b)

aout     = a
bout     = b

allocate(stack(2,maxstack))

allocate(ts0(kcheb), qs0(kcheb), qders0(kcheb) )
allocate(ws0(kcheb),wders0(kcheb),wder2s0(kcheb))
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
   ts00(1) = b
   call fun(1,ts00,qs00,qders00,par1,par2,par3)
   qval = qs00(1)
!   call chebpw_interp(chebcoefs,qs,b,qval)


   apppval = ( 4*qval*apval**2-4*apval**4+3*appval**2) / (2*apval)

   w2   = 1/apval
   wp2  = -appval/apval**2
   wpp2 = 2*appval**2/apval**3 - apppval/apval**2

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
    
    ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2)    )
   call fun(kcheb,ts0,qs0,qders0,par1,par2,par3)
!   call chebpw_interp(chebcoefs,qs,qders,kcheb,ts0,qs0,qders0)
   ifaccept = 0
   dq       = 1d300
   dcoefs   = 1d300
   jer      = -1

   call appell_nop_fit(vars,qs0,dq)

   if (dq .lt. eps) then

   call appell_ivp_specint_nop(vars,jer,a0,b0,kcheb,ts0,qs0,qders0,w2,wp2,wpp2, &
     ws0, wders0, wder2s0)


   if (jer .eq. 0) then
   alphapp0(:,nints2+1) = -wders0/ws0**2
   alphap0(:,nints2+1)  = 1.0d0/ws0
   alpha0(:,nints2+1)   = aval + (b0-a0)/2*matmul(vars%aintl,alphap0(:,nints2+1))
   
   call appell_nop_fit(vars,alphap0(:,nints2+1),dcoefs)

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
     w2               = ws0(kcheb)
     wp2              = wders0(kcheb)
     wpp2             = wder2s0(kcheb)
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
   ts00(1) = b
   call fun(1,ts00,qs00,qders00,par1,par2,par3)
   qval = qs00(1)
!   call chebpw_interp(chebcoefs,qs,a,qval)

   apppval = ( 4*qval*apval**2-4*apval**4+3*appval**2) / (2*apval)

   w1   = 1/apval
   wp1  = -appval/apval**2
   wpp1 = 2*appval**2/apval**3 - apppval/apval**2


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
    
    ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2)  )
!   call chebpw_interp(chebcoefs,qs,qders,kcheb,ts0,qs0,qders0)
   call fun(kcheb,ts0,qs0,qders0,par1,par2,par3)

   call appell_nop_fit(vars,qs0,dq)

   if (dq .lt. eps) then

   call appell_tvp_specint_nop(vars,jer,a0,b0,kcheb,ts0,qs0,qders0,w1,wp1,wpp1, &
     ws0, wders0, wder2s0)
   if (jer .eq. 0) then
   alphapp0(:,nints1-1) = -wders0/ws0**2
   alphap0(:,nints1-1)  = 1.0d0/ws0
   alpha0(:,nints1-1)   = aval + (b0-a0)/2*matmul(vars%aintr,alphap0(:,nints1-1))


   call appell_nop_fit(vars,alphap0(:,nints1-1),dcoefs)
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

     w1               = ws0(1)
     wp1              = wders0(1)
     wpp1             = wder2s0(1)
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
nints = nintsnew

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


subroutine appell_nop_solve1(vars,ier,eps,a,b,c,w,wp,wpp,chebcoefs,qs,qders, &
  chebsol,ws,wders,wder2s,nints00,ab00)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)                            :: vars
type(chebpw_scheme)                                :: chebcoefs
double precision                                   :: qs(:), qders(:)
type(chebpw_scheme), intent(out)                   :: chebsol
double precision, allocatable, intent(out)         :: ws(:), wders(:), wder2s(:)

integer, optional                                  :: nints00
double precision, optional                         :: ab00(:,:)

!
!  Find a solution of (1) such given the values of w(t) and its first two derivatives
!  at a point c in [a,b].
!
!  The coefficient and its derivative are specified by piecewise Chebyshev expansions.
!
!  The user has the option of providing an initial list of intervals for the 
!  discretization (this list will be adaptively refine).   In this event, the 
!  parameters a and b specifying the solution interval will be ignored.
!
!  This routine truncates the interval [a,b] in the event that the solution
!  of Appell's equation overflows or underflows. 
!
!  Input parameters:
!    vars - the structure populated by appell_nop_init
!    eps - the desired precision for the solution
!    c - the point at which the values of w and its derivatives is given
!    w - the desired value of w(t)
!    wpa - the desired value of w'(t)
!    wppa - the desired value of w''(t)
!    chebcoefs - the piecewise Chebyshev discretization scheme used to
!      represent the coefficients
!    qs, qders - the vector of values of q(t) and q'(t) at the discretization nodes
!
!    nints00,ab00 - an optional list of initial intervals for the discretization
!     of the solution 
!
!         NOTE: when these are specified, the interval (a,b) is ignored
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
!    chebsol - a structure describing the piecewise Chebyshev discretization scheme
!      used to represent the solution
!    ws - the values of the solution at the piecewise Chebyshev nodes
!    wders - the values of the derivative of the solution at the piecewise 
!      Chebyshev nodes
!    wder2s -  the values of the second derivative of the solution at the 
!      piecewise Chebyshev nodes
!

double precision, allocatable      :: stack(:,:), ab(:,:)

double precision, allocatable      :: ws1(:,:), wders1(:,:), wder2s1(:,:)
double precision, allocatable      :: ws2(:,:), wders2(:,:), wder2s2(:,:)

double precision, allocatable      :: ab1(:,:), ab2(:,:), vals0(:)
double precision, allocatable      :: ts0(:), qders0(:), qs0(:)

ier      = 0

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
dmax     = vars%dmax

allocate(ab1(2,maxints))
allocate(ws1(kcheb,maxints+1))
allocate(wders1(kcheb,maxints+1))
allocate(wder2s1(kcheb,maxints+1))

allocate(ab2(2,maxints))
allocate(ws2(kcheb,maxints+1))
allocate(wders2(kcheb,maxints+1))
allocate(wder2s2(kcheb,maxints+1))

allocate(vals0(kcheb))
allocate(qs0(kcheb),ts0(kcheb),qders0(kcheb))

maxstack = maxints
allocate(stack(2,maxstack))

nints1     = 0
nints2     = 0

!
!  Solve on (a,c)
!

ind             = 0

if (c .gt. a) then


!
!  Form the initial list of intervals
!


if (present(nints00)) then

   nstack     = 0
   do int=1,nints00
   a00 = ab00(1,int)
   b00 = ab00(2,int)
    
   if (a00 .ge. c) exit
    
   if (b00 .le. c) then
   nstack = nstack+1
   stack(1,nstack) = a00
   stack(2,nstack) = b00
   else if (a00 .lt. c) then
   nstack = nstack+1
   stack(1,nstack) = a00
   stack(2,nstack) = c
   endif
   end do

else


  ! nstack           = 1
  ! stack(1,nstack)  = a
  ! stack(2,nstack)  = c

!
!  Use the discretization intervals from q(t)
!
   nstack     = 0
   do int=1,chebcoefs%nints
   a00 = chebcoefs%ab(1,int)
   b00 = chebcoefs%ab(2,int)
    
   if (b00 .le. a) cycle
   if (a00 .ge. c) exit
    
   if (b00 .le. c) then
   nstack = nstack+1
   stack(1,nstack) = max(a,a00)
   stack(2,nstack) = b00
   else if (a00 .lt. c) then
   nstack = nstack+1
   stack(1,nstack) = max(a,a00)
   stack(2,nstack) = c
   endif
   end do

endif


do while ( nstack > 0 )
a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
call chebpw_interp(chebcoefs,qs,qders,kcheb,ts0,qs0,qders0)

if (nints1 .eq. 0) then
wb0   = w
wpb0  = wp
wppb0 = wpp
else
wb0   = ws1(1,nints1)
wpb0  = wders1(1,nints1)
wppb0 = wder2s1(1,nints1)
endif

ifexit   = 0
ifaccept = 0
dcoefs1  = 1d300
dcoefs2  = 1d300
dq       =-1
jer      = 0

!if (present(nints00)) call appell_nop_fit(qs0,dq)

! call appell_nop_fit(vars,qs0,dq)

! if (dq .lt. eps) then

   call appell_tvp_specint_nop(vars,jer,a0,b0,kcheb,ts0,qs0,qders0,wb0,wpb0,wppb0,&
     ws1(:,nints1+1),wders1(:,nints1+1), wder2s1(:,nints1+1))
    
   if (jer .eq. 0) then
    
      vals0 = 1/ws1(:,nints1+1)
      call appell_nop_fit(vars,ws1(:,nints1+1),dcoefs1)
      call appell_nop_fit(vars,vals0,dcoefs2)
       
      if (dcoefs1 .lt. eps .AND. dcoefs2 .lt. eps) then
      ifaccept=1
      
      ! Quit if we exceed the maximum allowable size
!      if (maxval(ws1(:,nints1+1)) .gt. dmax) ifexit=1

      endif

   endif

  ! try to handle overflow or underflow
  if (jer .eq. 1024) ifexit = 1


! endif

! ind = ind+1
! write (*,"(I6,I6,6(D20.10,2X),I2)")  ind,jer,a0,b0,dstep,dq,dcoefs1,dcoefs2,ifaccept
! write (13,"(I6,I6,6(D20.10,2X),I2)")  ind,jer,a0,b0,dstep,dq,dcoefs1,dcoefs2,ifaccept


if (ifaccept .eq. 1) then

if (nints1+1 .gt. maxints) then
ier = 4
return
endif

nints1        = nints1+1
ab1(1,nints1) = a0
ab1(2,nints1) = b0

else


if (nstack+2 .gt. maxstack) then
ier = 8
return
endif

nstack           = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

nstack           = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

endif


if (ifexit .eq. 1) exit

end do

endif


!
!  Solve on the interval (c,b)
!

if (b .gt. c) then
ind        = 0



!
!  From an initial list of intervals
!

if (present(nints00)) then

   nstack     = 0
   do int=nints00,1,-1
   a00 = ab00(1,int)
   b00 = ab00(2,int)
    
   if (a00 .ge. c) then
   nstack = nstack+1
   stack(1,nstack) = a00
   stack(2,nstack) = b00
   else if (c .lt. b00) then
   nstack = nstack+1
   stack(1,nstack) = c
   stack(2,nstack) = b00
   endif
 
   end do

else

  ! nstack          = 1
  ! stack(1,nstack) = c
  ! stack(2,nstack) = b

!
!  Use the discretization intervals from q(t)
!
   nstack     = 0
   do int=chebcoefs%nints,1,-1
   a00 = chebcoefs%ab(1,int)
   b00 = chebcoefs%ab(2,int)
    
   if (a00 .ge. b) exit
    
   if (a00 .ge. c) then
   nstack = nstack+1
   stack(1,nstack) = a00
   stack(2,nstack) = min(b00,b)
   else if (c .lt. b00) then
   nstack = nstack+1
   stack(1,nstack) = c
   stack(2,nstack) = min(b00,b)
   endif
 
   end do

endif


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
call chebpw_interp(chebcoefs,qs,qders,kcheb,ts0,qs0,qders0)

ifaccept = 0
ifexit   = 0
dcoefs1  = 1d300
dcoefs2  = 1d300
dq       = -1
jer      = 0

if (nints2 .eq. 0) then
wa0   = w
wpa0  = wp
wppa0 = wpp
else
wa0   = ws2(kcheb,nints2)
wpa0  = wders2(kcheb,nints2)
wppa0 = wder2s2(kcheb,nints2)
endif


!if (present(nints00)) call appell_nop_fit(qs0,dq)
! call appell_nop_fit(vars,qs0,dq)
! if (dq .lt. eps) then

   call appell_ivp_specint_nop(vars,jer,a0,b0,kcheb,ts0,qs0,qders0,wa0,wpa0,wppa0,&
     ws2(:,nints2+1),wders2(:,nints2+1), wder2s2(:,nints2+1))
    
   if (jer .eq. 0) then
      
      vals0 = 1/ws2(:,nints2+1)
      call appell_nop_fit(vars,ws2(:,nints2+1),dcoefs1)
      call appell_nop_fit(vars,vals0,dcoefs2)
      
      if (dcoefs1 .lt. eps .AND. dcoefs2 .lt. eps) then
      ifaccept=1
!      if (maxval(ws2(:,nints2+1)) .gt. dmax) ifexit = 1
      endif

   endif

  ! try to handle overflow or underflow
  if (jer .eq. 1024) ifexit = 1

! endif

! ind = ind+1
! write (*,"(I6,I6,4(D20.10,2X),I2)")   ind,jer,a0,b0,dcoefs1,dcoefs2,ifaccept
! write (13,"(I6,I6,4(D20.10,2X),I2)")  ind,jer,a0,b0,dcoefs1,dcoefs2,ifaccept

if (ifaccept .eq. 1) then

if (nints2+1 .gt. maxints) then
ier = 4
return
endif

nints2        = nints2+1
ab2(1,nints2) = a0
ab2(2,nints2) = b0


else

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif

nstack           = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack           = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2


endif



if (ifexit .eq. 1) exit

end do

endif

!
! Copy out the solution
!

nints = nints1+nints2
allocate(ab(2,1:nints))
ab(:,1:nints1)       = ab1(:,nints1:1:-1)
ab(:,nints1+1:nints) = ab2(:,1:nints2)
call chebpw_specified(chebsol,kcheb,nints,ab)

allocate(ws(kcheb*nints), wders(kcheb*nints), wder2s(kcheb*nints) )

n1 = kcheb*nints1
n2 = kcheb*nints

ws(1:n1)     = reshape(ws1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
ws(n1+1:n2)  = reshape(ws2(1:kcheb,1:nints2), [kcheb*nints2] )

wders(1:n1)     = reshape(wders1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
wders(n1+1:n2)  = reshape(wders2(1:kcheb,1:nints2), [kcheb*nints2] )

wder2s(1:n1)     = reshape(wder2s1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
wder2s(n1+1:n2)  = reshape(wder2s2(1:kcheb,1:nints2), [kcheb*nints2] )

end subroutine



subroutine appell_nop_solve2(vars,ier,eps,a,b,c,w,wp,wpp,fun,par1,par2,par3, &
  chebsol,ws,wders,wder2s,nints00,ab00)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)                            :: vars
procedure(appell_nop_fun1)                         :: fun
type(chebpw_scheme), intent(out)                   :: chebsol
double precision, allocatable, intent(out)         :: ws(:), wders(:), wder2s(:)

integer, optional                                  :: nints00
double precision, optional                         :: ab00(:,:)

!
!  Find a solution of (1) such given the values of w(t) and its first two derivatives
!  at a point c in [a,b].
!
!  The coefficient and its derivative are specified via an external subroutine.
!
!  The user has the option of providing an initial list of intervals for the 
!  discretization (this list will be adaptively refine).   In this event, the 
!  parameters a and b specifying the solution interval will be ignored.
!
!  In region where the solution of Appell's equation w(t) is greater than
!  the parameter dmax (set by appell_nop_init), the value of alpha'(t) = 1 / w(t) 
!  will be set to 0.
!
!  Input parameters:
!    vars - the structure populated by appell_nop_init
!    eps - the desired precision for the solution
!    c - the point at which the values of w and its derivatives is given
!    w - the desired value of w(t)
!    wpa - the desired value of w'(t)
!    wppa - the desired value of w''(t)
!    fun - an exteranl subroutine supplying the values of the coefficients
!    par? - user-supplied parameters which are passed to the external subroutine
!
!    nints00,ab00 - an optional list of initial intervals for the discretization
!     of the solution 
!
!         NOTE: when these are specified, the interval (a,b) is ignored
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
!    chebsol - a structure describing the piecewise Chebyshev discretization scheme
!      used to represent the solution
!    ws - the values of the solution at the piecewise Chebyshev nodes
!    wders - the values of the derivative of the solution at the piecewise 
!      Chebyshev nodes
!    wder2s -  the values of the second derivative of the solution at the 
!      piecewise Chebyshev nodes
!

double precision, allocatable      :: stack(:,:), ab(:,:)

double precision, allocatable      :: ws1(:,:), wders1(:,:), wder2s1(:,:)
double precision, allocatable      :: ws2(:,:), wders2(:,:), wder2s2(:,:)

double precision, allocatable      :: ab1(:,:), ab2(:,:), vals0(:)
double precision, allocatable      :: ts0(:), qders0(:), qs0(:)

ier      = 0

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
dmax     = vars%dmax

allocate(ab1(2,maxints))
allocate(ws1(kcheb,maxints+1))
allocate(wders1(kcheb,maxints+1))
allocate(wder2s1(kcheb,maxints+1))

allocate(ab2(2,maxints))
allocate(ws2(kcheb,maxints+1))
allocate(wders2(kcheb,maxints+1))
allocate(wder2s2(kcheb,maxints+1))

allocate(vals0(kcheb))
allocate(qs0(kcheb),ts0(kcheb),qders0(kcheb))

maxstack = maxints
allocate(stack(2,maxstack))

nints1     = 0
nints2     = 0


!
!  Solve on (a,c)
!

ind             = 0

if (c .gt. a) then


!
!  Form the initial list of intervals
!


if (present(nints00)) then

   nstack     = 0
   do int=1,nints00
   a00 = ab00(1,int)
   b00 = ab00(2,int)
    
   if (a00 .ge. c) exit
    
   if (b00 .le. c) then
   nstack = nstack+1
   stack(1,nstack) = a00
   stack(2,nstack) = b00
   else if (a00 .lt. c) then
   nstack = nstack+1
   stack(1,nstack) = a00
   stack(2,nstack) = c
   endif
   end do

else


  nstack           = 1
  stack(1,nstack)  = a
  stack(2,nstack)  = c


endif


do while ( nstack > 0 )
a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
call fun(kcheb,ts0,qs0,qders0,par1,par2,par3)


if (nints1 .eq. 0) then
wb0   = w
wpb0  = wp
wppb0 = wpp
else
wb0   = ws1(1,nints1)
wpb0  = wders1(1,nints1)
wppb0 = wder2s1(1,nints1)
endif

ifexit   = 0
ifaccept = 0
dcoefs1  = 1d300
dcoefs2  = 1d300
dq       =-1
jer      = 0

call appell_nop_fit(vars,qs0,dq)

if (dq .lt. eps) then

   call appell_tvp_specint_nop(vars,jer,a0,b0,kcheb,ts0,qs0,qders0,wb0,wpb0,wppb0,&
     ws1(:,nints1+1),wders1(:,nints1+1), wder2s1(:,nints1+1))

    
   if (jer .eq. 0) then
    
      vals0 = 1/ws1(:,nints1+1)
      call appell_nop_fit(vars,ws1(:,nints1+1),dcoefs1)
      call appell_nop_fit(vars,vals0,dcoefs2)
       
      if (dcoefs1 .lt. eps .AND. dcoefs2 .lt. eps) then
      ifaccept=1
      
      ! Quit if we exceed the maximum allowable size
      if (maxval(ws1(:,nints1+1)) .gt. dmax) ifexit=1
      
      endif

   endif

  ! try to handle overflow or underflow
  if (jer .eq. 1024) ifexit = 1

endif

! ind = ind+1
! write (*,"(I6,I6,5(D20.10,2X),I2)")  ind,jer,a0,b0,dq,dcoefs1,dcoefs2,ifaccept
! write (13,"(I6,I6,5(D20.10,2X),I2)")  ind,jer,a0,b0,dq,dcoefs1,dcoefs2,ifaccept



if (ifaccept .eq. 1) then

if (nints1+1 .gt. maxints) then
ier = 4
return
endif

nints1        = nints1+1
ab1(1,nints1) = a0
ab1(2,nints1) = b0

else


if (nstack+2 .gt. maxstack) then
ier = 8
return
endif

nstack           = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

nstack           = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

endif


if (ifexit .eq. 1) exit

end do

endif



!
!  Solve on the interval (c,b)
!

if (b .gt. c) then
ind        = 0



!
!  From an initial list of intervals
!

if (present(nints00)) then

   nstack     = 0
   do int=nints00,1,-1
   a00 = ab00(1,int)
   b00 = ab00(2,int)
    
   if (a00 .ge. c) then
   nstack = nstack+1
   stack(1,nstack) = a00
   stack(2,nstack) = b00
   else if (c .lt. b00) then
   nstack = nstack+1
   stack(1,nstack) = c
   stack(2,nstack) = b00
   endif
 
   end do

else

  nstack          = 1
  stack(1,nstack) = c
  stack(2,nstack) = b

 
endif


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack = nstack-1

ts0 = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
call fun(kcheb,ts0,qs0,qders0,par1,par2,par3)

ifaccept = 0
ifexit   = 0
dcoefs1  = 1d300
dcoefs2  = 1d300
dq       = -1
jer      = 0

if (nints2 .eq. 0) then
wa0   = w
wpa0  = wp
wppa0 = wpp
else
wa0   = ws2(kcheb,nints2)
wpa0  = wders2(kcheb,nints2)
wppa0 = wder2s2(kcheb,nints2)
endif


call appell_nop_fit(vars,qs0,dq)


if (dq .lt. eps) then

   call appell_ivp_specint_nop(vars,jer,a0,b0,kcheb,ts0,qs0,qders0,wa0,wpa0,wppa0,&
     ws2(:,nints2+1),wders2(:,nints2+1), wder2s2(:,nints2+1))
    
   if (jer .eq. 0) then
      
      vals0 = 1/ws2(:,nints2+1)
      call appell_nop_fit(vars,ws2(:,nints2+1),dcoefs1)
      call appell_nop_fit(vars,vals0,dcoefs2)
      
      if (dcoefs1 .lt. eps .AND. dcoefs2 .lt. eps) then
      ifaccept=1
!      if (maxval(ws2(:,nints2+1)) .gt. dmax) ifexit = 1
      endif

   endif


  ! try to handle overflow or underflow
  if (jer .eq. 1024) ifexit = 1

endif

! ind = ind+1
! write (*,"(I6,I6,4(D20.10,2X),I2)")   ind,jer,a0,b0,dcoefs1,dcoefs2,ifaccept
! write (13,"(I6,I6,4(D20.10,2X),I2)")  ind,jer,a0,b0,dcoefs1,dcoefs2,ifaccept

if (ifaccept .eq. 1) then

if (nints2+1 .gt. maxints) then
ier = 4
return
endif

nints2        = nints2+1
ab2(1,nints2) = a0
ab2(2,nints2) = b0


else

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif

nstack           = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack           = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2


endif




if (ifexit .eq. 1) exit

end do

endif


!
! Copy out the solution
!

nints = nints1+nints2
allocate(ab(2,1:nints))
ab(:,1:nints1)       = ab1(:,nints1:1:-1)
ab(:,nints1+1:nints) = ab2(:,1:nints2)
call chebpw_specified(chebsol,kcheb,nints,ab)

allocate(ws(kcheb*nints), wders(kcheb*nints), wder2s(kcheb*nints) )

n1 = kcheb*nints1
n2 = kcheb*nints

ws(1:n1)     = reshape(ws1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
ws(n1+1:n2)  = reshape(ws2(1:kcheb,1:nints2), [kcheb*nints2] )

wders(1:n1)     = reshape(wders1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
wders(n1+1:n2)  = reshape(wders2(1:kcheb,1:nints2), [kcheb*nints2] )

wder2s(1:n1)     = reshape(wder2s1(1:kcheb,nints1:1:-1), [kcheb*nints1] )
wder2s(n1+1:n2)  = reshape(wder2s2(1:kcheb,1:nints2), [kcheb*nints2] )

end subroutine




subroutine appell_nop_phase1(vars,ier,eps,a,b,c,apval,appval,d,aval,chebcoefs,qs,qders, &
  chebphase,alpha,alphader,alphader2,nints00,ab00)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)                            :: vars
type(chebpw_scheme)                                :: chebcoefs
double precision                                   :: qs(:), qders(:)
type(chebpw_scheme), intent(out)                   :: chebphase
double precision, allocatable, intent(out)         :: alpha(:), alphader(:), alphader2(:)
integer, optional                                  :: nints00
double precision, optional                         :: ab00(:,:)
!
!  Given the values of the first two derivatives of the phase function alpha(t) at
!  a point c in the interval [a,b] and the value of alpha(t) at the point d in [a,b],
!  construct piecewise Chebyshev expansions of alpha(t) and its first two derivatives
!  on the interval [a,b] by solving Appell's equation.
!
!  The coefficient and its derivative are specified by piecewise Chebyshev expansions.
!
!  The user has the option of providing an initial list of intervals for the 
!  discretization (this list will be adaptively refine).   In this event, the 
!  parameters a and b specifying the solution interval will be ignored.
! 
!  This routine will truncate the solution domain if w(t) grows larger than
!  the parameter dmax which is set by appell_nop_init.  In this case, if d 
!  is outside the truncated solution domain, it will be set to the closest endpoint.
!
!  Input parameters:
!    vars - the structure populated by appell_nop_init
!    eps - the desired precision for the solution
!    c - the point at which the values of alpha'(t) and alpha''(t) are given
!    apval - the value of alpha'(c)
!    appval - the value of alpha''(c)
!    chebcoefs - the piecewise Chebyshev discretization scheme used to
!      represent the coefficients
!    qs, qders - the vector of values of q(t) and q'(t) at the discretization nodes
!
!    nints00,ab00 - an optional list of initial intervals for the discretization
!     of the solution 
!
!         NOTE: when these are specified, the interval (a,b) is ignored
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

double precision, allocatable        :: ws(:), wders(:), wder2s(:)

ier = 0

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
dmax     = vars%dmax

!  
!  Calculate the initial values for w(t) at c ... the value of alpha''(t) is 
!  determined via the Riccati/Kummer equation
! 

call chebpw_interp(chebcoefs,qs,c,qval)
apppval = ( 4*qval*apval**2-4*apval**4+3*appval**2) / (2*apval)

w      = 1.0d0/apval
wp     = -appval/apval**2
wpp    = 2*appval**2/apval**3 - apppval/apval**2

!
!  Call the solver
!

call appell_nop_solve(vars,ier,eps,a,b,c,w,wp,wpp,chebcoefs,qs,qders, &
  chebphase,ws,wders,wder2s,nints00,ab00)
if (ier .ne. 0) return


call chebpw_interval(chebphase,a0,b0)
call chebpw_info(chebphase,k0,nints0)
ncoefs = k0*nints0

allocate(alphader(ncoefs), alphader2(ncoefS))

alphader  = 1/ws
alphader2 = -wders/ws**2

!
!  Integrate alpha'(t)
!

d0 = max(d,a0)
d0 = min(d,b0)

call chebpw_int(chebphase,alphader,d0,aval,alpha)

end subroutine


subroutine appell_nop_phase2(vars,ier,eps,a,b,c,apval,appval,d,aval,fun,par1,par2,par3, &
  chebphase,alpha,alphader,alphader2,nints00,ab00)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)                            :: vars
procedure(appell_nop_fun1)                         :: fun
type(chebpw_scheme), intent(out)                   :: chebphase
double precision, allocatable, intent(out)         :: alpha(:), alphader(:), alphader2(:)
integer, optional                                  :: nints00
double precision, optional                         :: ab00(:,:)
!
!  Given the values of the first two derivatives of the phase function alpha(t) at
!  a point c in the interval [a,b] and the value of alpha(t) at the point d in [a,b],
!  construct piecewise Chebyshev expansions of alpha(t) and its first two derivatives
!  on the interval [a,b] by solving Appell's equation.
!
!  The coefficient and its derivative are specified via an external subroutine.
!
!  The user has the option of providing an initial list of intervals for the 
!  discretization (this list will be adaptively refine).   In this event, the 
!  parameters a and b specifying the solution interval will be ignored.
! 
!  This routine will truncate the solution domain if w(t) grows larger than
!  the parameter dmax which is set by appell_nop_init.  In this case, if d 
!  is outside the truncated solution domain, it will be set to the closest endpoint.
!
!  Input parameters:
!    vars - the structure populated by appell_nop_init
!    eps - the desired precision for the solution
!    c - the point at which the values of alpha'(t) and alpha''(t) are given
!    apval - the value of alpha'(c)
!    appval - the value of alpha''(c)
!    fun - an external subroutine which supplies  the values of the coefficients
!    par? - user-supplied parameters which are passed to the external subroutine
!
!    nints00,ab00 - an optional list of initial intervals for the discretization
!     of the solution 
!
!         NOTE: when these are specified, the interval (a,b) is ignored
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

double precision, allocatable        :: ws(:), wders(:), wder2s(:)
double precision                     :: ts0(1), qs0(1), qders0(1)

ier = 0

kcheb    = vars%kcheb
ntest    = vars%ntest
maxints  = vars%maxints
epsint   = vars%epsint
dmax     = vars%dmax

!  
!  Calculate the initial values for w(t) at c ... the value of alpha''(t) is 
!  determined via the Riccati/Kummer equation
! 

!call chebpw_interp(chebcoefs,qs,c,qval)
ts0(1) = c
call fun(1,ts0,qs0,qders0,par1,par2,par3)
qval = qs0(1)
apppval = ( 4*qval*apval**2-4*apval**4+3*appval**2) / (2*apval)


w      = 1.0d0/apval
wp     = -appval/apval**2
wpp    = 2*appval**2/apval**3 - apppval/apval**2

!
!  Call the solver
!

call appell_nop_solve(vars,ier,eps,a,b,c,w,wp,wpp,fun,par1,par2,par3, &
  chebphase,ws,wders,wder2s,nints00,ab00)
if (ier .ne. 0) return


call chebpw_interval(chebphase,a0,b0)
call chebpw_info(chebphase,k0,nints0)
ncoefs = k0*nints0

allocate(alphader(ncoefs), alphader2(ncoefS))

alphader  = 1/ws
alphader2 = -wders/ws**2

!
!  Integrate alpha'(t)
!

d0 = max(d,a0)
d0 = min(d,b0)

call chebpw_int(chebphase,alphader,d0,aval,alpha)

end subroutine



subroutine appell_ivp_specint_nop(vars,ier,a,b,kcheb,ts,qs,qders,wa,wpa,wppa,&
  ws,wders,wder2s)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)             :: vars
double precision                    :: ts(:), qs(:), qders(:)
double precision                    :: ws(:), wders(:),wder2s(:)
!
!  Approximate the solution of an initial value problem for (1)
!  at the practical Chebyshev nodes on the interval [a,b].
!
!  Input parameters:
!    (a,b) - the interval on which the IVP is given
!    ts - the practical Chebyshev nodes on [a,b]
!    qs - the values of the coefficient q(t) at the Chebyshev nodes
!    qders - the values of the derivative q'(t) at the Chebyshev nodes
!    wa - the initial value for the solution
!    wpa - the initial value for the derivative of the solution
!    wppa - the initial value for the second derivative of the solution
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0     indicates successful execution
!       ier = 1024  means that NaN or INF was encountered
!

double precision                    :: rs0(kcheb), ps0(kcheb), qs0(kcheb)

ier    = 0

ps0    = 4*qs
qs0    = 2*qders

call appell_ivp_specint_nop0(vars,a,b,kcheb,ts,ps0,qs0,wa,wpa,wppa,&
 ws,wders,wder2s)

do i=1,kcheb
if (.not. ieee_is_finite(ws(i)) .or.  ieee_is_nan(ws(i))) then
ier = 1024
return
endif
end do


end subroutine


subroutine appell_ivp_specint_nop0(vars,a,b,kcheb,ts,ps,qs,wa,wpa,wppa,&
  ws,wders,wder2s)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)             :: vars
double precision                    :: ts(:)
double precision                    :: ps(:), qs(:)
double precision                    :: ws(:), wders(:),wder2s(:)
!
!  Approximate the solution of the initial value problem
!
!    { w'''(t) + p(t) w'(t) +  q(t) w(t) = 0     for a < t < b 
!    { w(a)   = wa
!    { w'(a)  = wpa
!    { w''(a) = wppa
!
!  at the practical Chebyshev nodes on the interval [a,b] using a spectral 
!  integration method.
!
!  Input parameters:
!    (a,b) - the interval on which the IVP is given
!    ts - the practical Chebyshev nodes on [a,b]
!    ps, qs - the values of the coefficients at the Chebyshev nodes
!    wa - the initial value for the solution
!    wpa - the initial value for the derivative of the solution
!    wppa - the initial value for the second derivative of the solution
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0  indicates successful execution
!       ier = 4  means that NaN or Inf was encountered while solving the
!                equation
!
!    ws, wders, wder2s - the values of the solution of the differential
!       equation and its first two derivatives at the relevant points
!
!
double precision                    :: amatr(kcheb,kcheb), bmatr(kcheb,kcheb)
double precision                    :: rhs(kcheb), sigma(kcheb)

!
!  Use the representation of the solution
! 
!                                                  t
!    w(t) = wa + wpa (t-a) + wppa/2 (t-a)^2 + \int   (t-s)^2/2 sigma(s) ds
!                                                  a
!
!  to solve and initial value problem for (1).  This representation leads to the
!  integral equation
!
!                        t                         t
!    sigma(t) + r(t) \int   sigma(s) ds + p(t) \int  (t-s) sigma(s) ds +
!                        a                         a
!
!                        t
!              q(t)  \int   (t-s)^2/2 sigma(s) ds = rhs(t)
!                        a
!
!  where
!
!    rhs(t) = - wppa r(t) - ( wpa + wppa(t-a) )  p(t) - (wa + wpa (t-a) + wppa (t-a)^2/2) q(t)
!
!


!
!  Allocate memory for the procedure and setup some parameters.
!

amatr  = 0
sigma  = 0

dd     = (b-a)/2 
dd2    = dd*dd
dd3    = dd*dd2

do i=1,kcheb
amatr(i,i) = 1.0d0
end do


!
! Handle the p(t) * w'(t) term.
!
do i=1,kcheb
amatr(i,:) = amatr(i,:) + ps(i) * vars%aintl2(i,:)*dd2
sigma(i)   = sigma(i) - ps(i)*(wpa + wppa * (ts(i)-a))
end do

!
!  Handle the q(t) w(t) term
!

do i=1,kcheb
amatr(i,:) = amatr(i,:) + qs(i) * vars%aintl3(i,:)*dd3
sigma(i)   = sigma(i) - qs(i) * (wa + wpa*(ts(i)-a) + wppa * (ts(i)-a)**2/2)
end do

! Solve the system
call linalg0_solve(kcheb,amatr,sigma)

! Construct the solutions
wder2s = wppa + dd*matmul(vars%aintl,sigma)
wders  = wpa  + dd*matmul(vars%aintl,wder2s)
ws     = wa   + dd*matmul(vars%aintl,wders)

end subroutine


subroutine appell_tvp_specint_nop(vars,ier,a,b,kcheb,ts,qs,qders,wb,wpb,wppb,&
  ws,wders,wder2s)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)             :: vars
double precision                    :: ts(:), qs(:), qders(:)
double precision                    :: ws(:), wders(:),wder2s(:)
!
!  Approximate the solution of a terminal value problem for (1)
!  at the practical Chebyshev nodes on the interval [a,b].
!
!  Input parameters:
!    (a,b) - the interval on which the IVP is given
!    ts - the practical Chebyshev nodes on [a,b]
!    qs - the values of the coefficient q(t) at the Chebyshev nodes
!    qders - the values of the derivative q'(t) at the Chebyshev nodes
!    wa - the initial value for the solution
!    wpa - the initial value for the derivative of the solution
!    wppa - the initial value for the second derivative of the solution
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0     indicates successful execution
!       ier = 1024  means that NaN or Inf was encountered while solving the
!                   equation
!
!

double precision                    :: rs0(kcheb), ps0(kcheb), qs0(kcheb)

ier    = 0

ps0    = 4*qs
qs0    = 2*qders

call appell_tvp_specint_nop0(vars,a,b,kcheb,ts,ps0,qs0,wb,wpb,wppb,&
 ws,wders,wder2s)

do i=1,kcheb
if (.not. ieee_is_finite(ws(i)) .or. ieee_is_nan(ws(i))) then
ier = 1024
return
endif
end do

end subroutine


subroutine appell_tvp_specint_nop0(vars,a,b,kcheb,ts,ps,qs,wb,wpb,wppb,&
  ws,wders,wder2s)
implicit double precision (a-h,o-z)
type(appell_nop_vars_t)             :: vars
double precision                    :: ts(:)
double precision                    :: ps(:), qs(:)
double precision                    :: ws(:), wders(:),wder2s(:)
!
!  Approximate the solution of the terminal value problem
!
!    { w'''(t) + p(t) w'(t) +  q(t) w(t) = 0     for a < t < b 
!    { w(b)   = wb
!    { w'(b)  = wpb
!    { w''(b) = wppb
!
!  at the practical Chebyshev nodes on the interval [a,b] using a spectral 
!  integration method.
!
!  Input parameters:
!    (a,b) - the interval on which the IVP is given
!    ts - the practical Chebyshev nodes on [a,b]
!    rs, ps, qs - the values of the coefficients at the Chebyshev nodes
!    wb - the termainl value for the solution
!    wpb - the terminal value for the derivative of the solution
!    wppb - the terminal value for the second derivative of the solution
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0  indicates successful execution
!       ier = 4  means that NaN or Inf was encountered while solving the
!                equation
!
!    ws, wders, wder2s - the values of the solution of the terminal value
!       problem and its first two derivatives at the relevant points
!
!
double precision                    :: amatr(kcheb,kcheb), bmatr(kcheb,kcheb)
double precision                    :: rhs(kcheb), sigma(kcheb)

!
!  Use the representation of the solution
! 
!                                                  t
!    w(t) = wb + wpb (t-b) + wppb/2 (t-b)^2 + \int   (t-s)^2/2 sigma(s) ds
!                                                  b
!
!  to solve and initial value problem for (1).  This representation leads to the
!  integral equation
!
!                        t                         t
!    sigma(t) + r(t) \int   sigma(s) ds + p(t) \int  (t-s) sigma(s) ds +
!                        b                         b
!
!                        t
!              q(t)  \int   (t-s)^2/2 sigma(s) ds = rhs(t)
!                        b
!
!  where
!
!    rhs(t) = - wppb r(t) - ( wpb + wppb(t-b) )  p(t) - (wb + wpb (t-b) + wppb (t-b)^2/2) q(t)
!
!


!
!  Allocate memory for the procedure and setup some parameters.
!

amatr  = 0
sigma  = 0

dd     = (b-a)/2 
dd2    = dd*dd
dd3    = dd*dd2

do i=1,kcheb
amatr(i,i) = 1.0d0
end do


!
! Handle the p(t) * w'(t) term.
!
do i=1,kcheb
amatr(i,:) = amatr(i,:) + ps(i) * vars%aintr2(i,:)*dd2
sigma(i)   = sigma(i) - ps(i)*(wpb + wppb * (ts(i)-b))
end do

!
!  Handle the q(t) w(t) term
!

do i=1,kcheb
amatr(i,:) = amatr(i,:) + qs(i) * vars%aintr3(i,:)*dd3
sigma(i)   = sigma(i) - qs(i) * (wb + wpb*(ts(i)-b) + wppb * (ts(i)-b)**2/2)
end do

! Solve the system
call linalg0_solve(kcheb,amatr,sigma)

! Construct the solutions
wder2s = wppb + dd*matmul(vars%aintr,sigma)
wders  = wpb  + dd*matmul(vars%aintr,wder2s)
ws     = wb   + dd*matmul(vars%aintr,wders)

end subroutine


end module
