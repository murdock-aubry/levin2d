!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for efficiently evaluating integrals of the forms
!
!        b
!    \int   f(x) exp ( g(x) ) dx                                                          (1)
!        a
!
!        b
!    \int  f(x) cos( g(x) ) dx                                                            (2)
!        a
!
!  and
!
!        b
!    \int  f(x) sin( g(x) ) dx                                                            (3)
!        a
!
!  via the "adaptive Levin method."  The functions f(x) and g(x) are 
!  specified by the user via an external subroutine or piecewise Chebyshev
!  expansions.
!
!  In main versions of the routine in this module, the user specifies only
!  f(x) and g(x).  The function g'(x), which is required by the method, is
!  computed via spectral differentiation.  However, versions of the routines
!  which take g'(x) as input are provided.
!
!  The following subroutines should be regarded as publicly callable:
!
!    levin_init - initialize the code -- this entails populating a structure which
!      is then passed to the other codes in this file
!
!
!    levin_adapexp - adaptively compute an integral of the form (1) with f(x)
!      and g(x) comple-valued functions specified by the user
!
!    levin_adapcos - adaptively compute an integral of the form (2) with f(x)
!      and g(x) real-valued functions specified by the user
!
!    levin_adapsin - adaptively compute an integral of the form (3) with f(x)
!      and g(x) real-valued functions specified by the user
!
!
!
!    levin_adapexp_diff - adaptively compute an integral of the form (1) with f(x),
!      g(x) and g'(x) complex-valued functions specified by the user
!
!    levin_adapcos_diff - adaptively compute an integral of the form (2) with f(x),
!      g(x) and g'(x) real-valued functions specified by the user
!
!    levin_adapsin_diff - adaptively compute an integral of the form (3) with f(x),
!      g(x) and g'(x) real-valued functions specified by the user
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Observations / To Do:
!
!    - In the case g(x) = I dlambda, the versions of the routines which perform spectral
!      differentiation lose accuracy on the order of dlambda.
!
!    - Write the versions of these routines in which the input functions are
!      specified via Chebyshev exapnsions
!
!    - The adapcos and adapsin routines use essentially the same solve procedure as 
!      adapexp, and so are not much faster.  However, in some cases f(x) cos(g(x)) is
!      well-behaved whereas f(x) sin(g(x)) is badly behaved, or vice-versa.  In such
!      cases the adaptive integration procedure can be much faster when evaluating
!      only one of these two integrands. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module levin

use utils
use linalg0
use chebyshev
use chebpw
use iso_c_binding

interface

subroutine levin_fun1(n,ts,gs,fs,par1,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double complex              :: gs(:), fs(:)
end subroutine

subroutine levin_fun2(n,ts,gs,gders,fs,par1,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double complex              :: gs(:), gders(:), fs(:)
end subroutine

subroutine levin_fun3(n,ts,gs,fs,par1,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double precision             :: gs(:), fs(:)
end subroutine

subroutine levin_fun4(n,ts,gs,gders,fs,par1,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
doubleprecision             :: gs(:), gders(:), fs(:)
end subroutine

subroutine levin_fun5(n,ts,gs,fs,userptr)
use iso_c_binding
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double complex              :: gs(:), fs(:)
type(c_ptr)                 :: userptr
end subroutine

subroutine levin_fun6(n,ts,gs,fs,userptr)
use iso_c_binding
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double precision            :: gs(:), fs(:)
type(c_ptr)                 :: userptr
end subroutine


end interface


interface        levin_adapexp
module procedure levin_adapexp1
module procedure levin_adapexp2
module procedure levin_adapexp3
end interface    levin_adapexp

interface        levin_adapcos
module procedure levin_adapcos1
module procedure levin_adapcos2
module procedure levin_adapcos3
end interface    levin_adapcos

interface        levin_adapsin
module procedure levin_adapsin1
module procedure levin_adapsin2
module procedure levin_adapsin3
end interface    levin_adapsin


interface        levin_adapexp_diff
module procedure levin_adapexp_diff1
end interface    levin_adapexp_diff

interface        levin_adapcos_diff
module procedure levin_adapcos_diff1
end interface    levin_adapcos_diff

interface        levin_adapsin_diff
module procedure levin_adapsin_diff1
end interface    levin_adapsin_diff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The structure which stores any data needed by the procedures in this modules and
!  which is populated by the initialization routine.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type       levin_vars_t
integer                                    :: kcheb
integer                                    :: maxstack
double precision, allocatable              :: xscheb(:)
double precision, allocatable              :: whtscheb(:)
double precision, allocatable              :: adiff(:,:)
double precision, allocatable              :: acoefs(:,:)
end type   levin_vars_t

contains

subroutine levin_init(vars,k0,maxstack0)
implicit double precision (a-h,o-z)
type(levin_vars_t)                          :: vars
integer, optional                           :: k0, maxstack0
!
!  Initialize the structure containing any data needed by the other procedures 
!  in this module.
!
!  If vars is the only parameter passed to this routine, then reasonable defaults
!  are chosen.
!
!  Input parameters:
!    k0 - the number of points per interval for the Gauss-Legendre rule
!    maxstack0 - the maximum number of entries in the stack used by the adaptive
!     integration procedure -- this is roughly equal to the maximum number !     intervals into which the integration domain can be decomposed
!
!  Output parameters:
!    vars - the structure containing all of the data needed by the other procedures in
!      this module
!

if (.not. present(k0) ) then
kcheb    = 12
maxstack = 10000
else
kcheb    = k0
maxstack = maxstack0
endif

vars%kcheb    = kcheb
vars%maxstack = maxstack

call chebyshev_quad(kcheb,vars%xscheb,vars%whtscheb)
call chebyshev_diffmatrix(kcheb,vars%adiff)
call chebyshev_coefsmatrix(kcheb,vars%acoefs)

end subroutine


subroutine levin_adapexp1(vars,ier,eps,a,b,fun,par1,par2,par3,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun1)           :: fun
double complex                  :: val
!
!  Evaluate an integral of the form (1) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) nad conforms to the levin_fun interface
!    par? - parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double complex, allocatable      :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(4,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))



ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)
epssq      = eps**2

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_exp(vars,a,b,kcheb,ts0,fs0,gs0,valr,vali)


nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = valr
stack(4,1) = vali

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
valr0  = stack(3,nstack)
vali0  = stack(4,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_exp(vars,a0,c0,kcheb,ts0,fs0,gs0,valr1,vali1)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_exp(vars,c0,b0,kcheb,ts0,fs0,gs0,valr2,vali2)

dnorm = (valr0-valr1-valr2)**2 + (vali0-vali1-vali2)**2


if (dnorm .lt. epssq) then
val = val + valr0+ima*vali0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr1
stack(4,nstack) = vali1

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = valr2
stack(4,nstack) = vali2

endif

end do

end subroutine


subroutine levin_adapexp11(vars,ier,eps,a,b,fun,par1,par2,par3,val,nints)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun1)           :: fun
double complex                  :: val
!
!  Evaluate an integral of the form (1) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) nad conforms to the levin_fun interface
!    par? - parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double complex, allocatable      :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(4,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)
epssq      = eps**2

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_exp(vars,a,b,kcheb,ts0,fs0,gs0,valr,vali)


nints = 0

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = valr
stack(4,1) = vali

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
valr0  = stack(3,nstack)
vali0  = stack(4,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_exp(vars,a0,c0,kcheb,ts0,fs0,gs0,valr1,vali1)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_exp(vars,c0,b0,kcheb,ts0,fs0,gs0,valr2,vali2)

dnorm = (valr0-valr1-valr2)**2 + (vali0-vali1-vali2)**2


if (dnorm .lt. epssq ) then
val = val + valr0+ima*vali0
nints = nints+1
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr1
stack(4,nstack) = vali1

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = valr2
stack(4,nstack) = vali2

endif

end do

end subroutine


subroutine levin_adapexp2(vars,ier,eps,a,b,fun,userptr,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun5)           :: fun
type(c_ptr)                     :: userptr
double complex                  :: val
!
!  Evaluate an integral of the form (1) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) 
!    userptr - a "void *" pointer which is passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double complex, allocatable      :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(4,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)
epssq      = eps**2

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
call fun(kcheb,ts0,gs0,fs0,userptr)
call levin_oneint_exp(vars,a,b,kcheb,ts0,fs0,gs0,valr,vali)

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = valr
stack(4,1) = vali

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
valr0  = stack(3,nstack)
vali0  = stack(4,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
call fun(kcheb,ts0,gs0,fs0,userptr)
call levin_oneint_exp(vars,a0,c0,kcheb,ts0,fs0,gs0,valr1,vali1)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
call fun(kcheb,ts0,gs0,fs0,userptr)
call levin_oneint_exp(vars,c0,b0,kcheb,ts0,fs0,gs0,valr2,vali2)

dnorm = (valr0-valr1-valr2)**2 + (vali0-vali1-vali2)**2


if (dnorm .lt. epssq) then
val = val + valr0+ima*vali0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr1
stack(4,nstack) = vali1

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = valr2
stack(4,nstack) = vali2

endif

end do

end subroutine


subroutine levin_adapexp3(vars,ier,eps,a,b,chebin,gs,fs,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
type(chebpw_scheme)             :: chebin
double complex                  :: gs(:), fs(:)
type(c_ptr)                     :: userptr
double complex                  :: val
!
!  Evaluate an integral of the form (1) given piecewise Chebyshev expansions
!  of functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    chebin - a structure describing the piecewise Chebyshev scheme
!     used to represent f(x) and g(x)
!    gs, fs - the values of the functions f(x) and g(x) at the
!     discretization nodes
!  
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double complex, allocatable      :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(4,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)
epssq      = eps**2

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)


! call fun(kcheb,ts0,gs0,fs0,userptr)
call chebpw_interp(chebin,gs,fs,kcheb,ts0,gs0,fs0)
call levin_oneint_exp(vars,a,b,kcheb,ts0,fs0,gs0,valr,vali)

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = valr
stack(4,1) = vali

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
valr0  = stack(3,nstack)
vali0  = stack(4,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
call chebpw_interp(chebin,gs,fs,kcheb,ts0,gs0,fs0)
call levin_oneint_exp(vars,a0,c0,kcheb,ts0,fs0,gs0,valr1,vali1)



ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
call chebpw_interp(chebin,gs,fs,kcheb,ts0,gs0,fs0)
call levin_oneint_exp(vars,c0,b0,kcheb,ts0,fs0,gs0,valr2,vali2)

dnorm = (valr0-valr1-valr2)**2 + (vali0-vali1-vali2)**2


if (dnorm .lt. epssq) then
val = val + valr0+ima*vali0
else


if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr1
stack(4,nstack) = vali1

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = valr2
stack(4,nstack) = vali2

endif

end do

end subroutine



subroutine levin_adapexp4(vars,ier,eps,a,b,gs,fs,val)
    implicit double precision (a-h,o-z)
    type(levin_vars_t)              :: vars
    double complex                  :: gs(:), fs(:)
    type(c_ptr)                     :: userptr
    double complex                  :: val
    !
    !  Evaluate an integral of the form (1) given piecewise Chebyshev expansions
    !  of functions f(x) and g(x).
    !
    !  Input parameters:
    !    vars - the structure populated by the initialization code
    !    eps - the desired precision for the calculations
    !    chebin - a structure describing the piecewise Chebyshev scheme
    !     used to represent f(x) and g(x)
    !    gs, fs - the values of the functions f(x) and g(x) at the
    !     discretization nodes
    !  
    !
    !  Output parameters:
    !    ier - an error return code;
    !      ier = 0  indicates successful execution
    !      ier = 4  means that the maximum number of intervals was exceeded while
    !        attempting to compute the integral
    !
    !    val - the value of the integral as calculated by this routine
    !
    
    double precision, allocatable    :: stack(:,:)
    double precision, allocatable    :: ts0(:), ts1(:)
    double complex, allocatable      :: gs0(:), fs0(:), gs01(:), fs01(:)
    double complex                   :: ima
    
    kcheb    = vars%kcheb
    maxstack = 10000!vars%maxstack
    
    
    allocate(stack(4,maxstack) )
    allocate(ts0(kcheb), fs0(kcheb), gs0(kcheb))
    allocate(ts1(kcheb))
    
    ier        = 0
    val        = 0
    ind        = 0
    ima        = (0.0d0,1.0d0)
    epssq      = eps ** 2
    

    ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)

    call chebyshev_interp5(kcheb, ts0, fs, gs, ts0, fs0, gs0)
    call levin_oneint_exp(vars,a,b,kcheb,ts0,fs,gs,valr,vali)
    
    nstack     = 1
    stack(1,1) = a
    stack(2,1) = b
    stack(3,1) = valr
    stack(4,1) = vali
    
    do while (nstack > 0)
    
    a0     = stack(1,nstack)
    b0     = stack(2,nstack)
    valr0  = stack(3,nstack)
    vali0  = stack(4,nstack)
    nstack = nstack-1
    
    
    c0     = (a0+b0)/2

    ts1    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
    call chebyshev_interp5(kcheb, ts0, fs, gs, ts1, fs0, gs0)
    call levin_oneint_exp(vars,a0,c0,kcheb,ts1,fs0,gs0,valr1,vali1)
    
    ts1    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
    call chebyshev_interp5(kcheb, ts0, fs, gs, ts1, fs0, gs0)
    call levin_oneint_exp(vars,c0,b0,kcheb,ts1,fs0,gs0,valr2,vali2)

    dnorm = (valr0-valr1-valr2)**2 + (vali0-vali1-vali2)**2

    if (dnorm .lt. epssq) then
        val = val + valr0+ima*vali0
    else
        
        if (nstack+2 .gt. maxstack) then
            ier = 4
            return
        endif
    
        nstack          = nstack+1
        stack(1,nstack) = a0
        stack(2,nstack) = c0
        stack(3,nstack) = valr1
        stack(4,nstack) = vali1

        
        nstack          = nstack+1
        stack(1,nstack) = c0
        stack(2,nstack) = b0
        stack(3,nstack) = valr2
        stack(4,nstack) = vali2
    endif
    
    end do
    
    end subroutine


subroutine levin_oneint_exp(vars,a,b,kcheb,ts,fs,gs,valr,vali)
implicit double precision (a-h,o-z)
type(levin_vars_t)        :: vars
double complex            :: gs(kcheb), fs(kcheb)
double complex            :: gders(kcheb), ps(kcheb)
double complex            :: amatr(kcheb,kcheb), coefs(kcheb)
double precision          :: ts(kcheb)
double complex            :: val

amatr = 2.0d0/(b-a)*vars%adiff
gders = matmul(amatr,gs)

do i=1,kcheb
amatr(i,i) = amatr(i,i) + gders(i)
end do


! call linalg0_solve(kcheb,amatr,fs)
! ps = fs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! call prin2("ts    = ",ts)
! call prinz("gs    = ",gs)
! call prinz("gders = ",gders)

!eps = epsilon(0.0d0)*10*norm2(abs(amatr))
eps = epsilon(0.0d0)*10*min(norm2(abs(amatr)),1.0d0)
call linalg0_qrtsolve(eps,kcheb,kcheb,amatr,ps,fs)

! coefs = matmul(vars%acoefs,ps)
! call prinz("coefs = ",coefs)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

val  = ps(kcheb)*exp(gs(kcheb)) - ps(1)*exp(gs(1))


valr = real(val)
vali = imag(val)

end subroutine


subroutine levin_adapcos1(vars,ier,eps,a,b,fun,par1,par2,par3,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun3)           :: fun
double precision                :: val
!
!  Evaluate an integral of the form (2) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) 
!    par? - parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double precision, allocatable    :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(3,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_cos(vars,a,b,kcheb,ts0,fs0,gs0,val0)

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = val0

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
val0   = stack(3,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_cos(vars,a0,c0,kcheb,ts0,fs0,gs0,valr)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_cos(vars,c0,b0,kcheb,ts0,fs0,gs0,vall)


dnorm = abs(val0-valr-vall)


if (dnorm .lt. eps) then
val = val + val0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif


nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = vall

endif

end do

end subroutine

subroutine levin_adapcos2(vars,ier,eps,a,b,fun,userptr,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun6)           :: fun
type(c_ptr)                     :: userptr
double precision                :: val
!
!  Evaluate an integral of the form (2) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) 
!    userptr - a "void *" pointer which is passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double precision, allocatable    :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(3,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
call fun(kcheb,ts0,gs0,fs0,userptr)
call levin_oneint_cos(vars,a,b,kcheb,ts0,fs0,gs0,val0)

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = val0

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
val0   = stack(3,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
call fun(kcheb,ts0,gs0,fs0,userptr)
call levin_oneint_cos(vars,a0,c0,kcheb,ts0,fs0,gs0,valr)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
call fun(kcheb,ts0,gs0,fs0,userptr)
call levin_oneint_cos(vars,c0,b0,kcheb,ts0,fs0,gs0,vall)


dnorm = abs(val0-valr-vall)


if (dnorm .lt. eps) then
val = val + val0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = vall

endif

end do

end subroutine

subroutine levin_adapcos3(vars,ier,eps,a,b,chebin,gs,fs,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
type(chebpw_scheme)             :: chebin
double precision                :: gs(:),fs(:)
type(c_ptr)                     :: userptr
double precision                :: val
!
!  Evaluate an integral of the form (2) with f(x) and g(x) specified
!  via piecewise Chebyshev expansions.
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    chebin - a structure describing the piecewise Chebyshev expansion
!      used to represent f(x) and g(x)
!    gs, fs - the values of g(x) and f(x) at the discretization nodes
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double precision, allocatable    :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(3,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
!call fun(kcheb,ts0,gs0,fs0,userptr)
call chebpw_interp(chebin,gs,fs,kcheb,ts0,gs0,fs0)
call levin_oneint_cos(vars,a,b,kcheb,ts0,fs0,gs0,val0)

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = val0

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
val0   = stack(3,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
!call fun(kcheb,ts0,gs0,fs0,userptr)
call chebpw_interp(chebin,gs,fs,kcheb,ts0,gs0,fs0)
call levin_oneint_cos(vars,a0,c0,kcheb,ts0,fs0,gs0,valr)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
!call fun(kcheb,ts0,gs0,fs0,userptr)
call chebpw_interp(chebin,gs,fs,kcheb,ts0,gs0,fs0)
call levin_oneint_cos(vars,c0,b0,kcheb,ts0,fs0,gs0,vall)


dnorm = abs(val0-valr-vall)


if (dnorm .lt. eps) then
val = val + val0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = vall

endif

end do

end subroutine


subroutine levin_oneint_cos(vars,a,b,kcheb,ts,fs,gs,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)        :: vars
double precision          :: gs(kcheb), fs(kcheb)
double complex            :: ps(kcheb), amatr(kcheb,kcheb), gders(kcheb), rhs(kcheb)
double precision          :: ts(kcheb)
double complex            :: ima, valz

ima = (0.0d0,1.0d0)

amatr = 2.0d0/(b-a)*vars%adiff
gders = ima*2.0d0/(b-a)*matmul(vars%adiff,gs)

do i=1,kcheb
amatr(i,i) = amatr(i,i) + gders(i)
end do

! rhs = fs
! call linalg0_solve(kcheb,amatr,rhs)

rhs = fs
eps = epsilon(0.0d0)*10
call linalg0_qrtsolve(eps,kcheb,kcheb,amatr,ps,rhs)


valz  = ps(kcheb)*exp(ima*gs(kcheb)) - ps(1)*exp(ima*gs(1))
val   = real(valz)

end subroutine




subroutine levin_adapsin11(vars,ier,eps,a,b,fun,par1,par2,par3,val,nints)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun3)           :: fun
double precision                :: val
!
!  Evaluate an integral of the form (2) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) 
!    par? - parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double precision, allocatable    :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(3,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_sin(vars,a,b,kcheb,ts0,fs0,gs0,val0)

nints =0 

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = val0

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
val0   = stack(3,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_sin(vars,a0,c0,kcheb,ts0,fs0,gs0,valr)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_sin(vars,c0,b0,kcheb,ts0,fs0,gs0,vall)


dnorm = abs(val0-valr-vall)


if (dnorm .lt. eps) then
val   = val + val0
nints = nints+1
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = vall

endif

end do

end subroutine



subroutine levin_adapsin1(vars,ier,eps,a,b,fun,par1,par2,par3,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun3)           :: fun
double precision                :: val
!
!  Evaluate an integral of the form (2) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) 
!    par? - parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double precision, allocatable    :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(3,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_sin(vars,a,b,kcheb,ts0,fs0,gs0,val0)


nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = val0

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
val0   = stack(3,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_sin(vars,a0,c0,kcheb,ts0,fs0,gs0,valr)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
call fun(kcheb,ts0,gs0,fs0,par1,par2,par3)
call levin_oneint_sin(vars,c0,b0,kcheb,ts0,fs0,gs0,vall)


dnorm = abs(val0-valr-vall)


if (dnorm .lt. eps) then
val = val + val0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = vall

endif

end do

end subroutine


subroutine levin_adapsin2(vars,ier,eps,a,b,fun,userptr,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun6)           :: fun
type(c_ptr)                     :: userptr
double precision                :: val
!
!  Evaluate an integral of the form (3) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) 
!    userptr - a "void *" pointer which is passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double precision, allocatable    :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(3,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)


ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
call fun(kcheb,ts0,gs0,fs0,userptr)
call levin_oneint_sin(vars,a,b,kcheb,ts0,fs0,gs0,val0)

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = val0

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
val0   = stack(3,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
call fun(kcheb,ts0,gs0,fs0,userptr)
call levin_oneint_sin(vars,a0,c0,kcheb,ts0,fs0,gs0,valr)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
call fun(kcheb,ts0,gs0,fs0,userptr)
call levin_oneint_sin(vars,c0,b0,kcheb,ts0,fs0,gs0,vall)


dnorm = abs(val0-valr-vall)



if (dnorm .lt. eps) then
val = val + val0

else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = vall

endif

end do

end subroutine


subroutine levin_adapsin3(vars,ier,eps,a,b,chebin,gs,fs,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
type(chebpw_scheme)             :: chebin
double precision                :: gs(:),fs(:)
type(c_ptr)                     :: userptr
double precision                :: val
!
!  Evaluate an integral of the form (3) with f(x) and g(x) specified
!  via piecewise Chebyshev expansions.
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    chebin - a structure describing the piecewise Chebyshev expansion
!      used to represent f(x) and g(x)
!    gs, fs - the values of g(x) and f(x) at the discretization nodes
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double precision, allocatable    :: gs0(:), fs0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(3,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)

ts0 = min(max(a,(b-a)/2 * vars%xscheb + (b+a)/2),b)
!call fun(kcheb,ts0,gs0,fs0,userptr)
call chebpw_interp(chebin,gs,fs,kcheb,ts0,gs0,fs0)
call levin_oneint_sin(vars,a,b,kcheb,ts0,fs0,gs0,val0)

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = val0

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
val0   = stack(3,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = min(c0,max(a0,(c0-a0)/2 * vars%xscheb + (c0+a0)/2))
!call fun(kcheb,ts0,gs0,fs0,userptr)
call chebpw_interp(chebin,gs,fs,kcheb,ts0,gs0,fs0)
call levin_oneint_sin(vars,a0,c0,kcheb,ts0,fs0,gs0,valr)

ts0    = min(b0,max(c0,(b0-c0)/2 * vars%xscheb + (b0+c0)/2))
!call fun(kcheb,ts0,gs0,fs0,userptr)
call chebpw_interp(chebin,gs,fs,kcheb,ts0,gs0,fs0)
call levin_oneint_sin(vars,c0,b0,kcheb,ts0,fs0,gs0,vall)

dnorm = abs(val0-valr-vall)


if (dnorm .lt. eps) then
val = val + val0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = vall

endif

end do

end subroutine

subroutine levin_oneint_sin(vars,a,b,kcheb,ts,fs,gs,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)        :: vars
double precision          :: gs(kcheb), fs(kcheb)
double complex            :: ps(kcheb), amatr(kcheb,kcheb), gders(kcheb), rhs(kcheb)
double precision          :: ts(kcheb)
double complex            :: ima, valz


ima = (0.0d0,1.0d0)

amatr = 2.0d0/(b-a)*vars%adiff
gders = ima*2.0d0/(b-a)*matmul(vars%adiff,gs)

do i=1,kcheb
amatr(i,i) = amatr(i,i) + gders(i)
end do

! ps = fs
! call linalg0_solve(kcheb,amatr,ps)

rhs = fs
eps = epsilon(0.0d0)*10
call linalg0_qrtsolve(eps,kcheb,kcheb,amatr,ps,rhs)

valz  = ps(kcheb)*exp(ima*gs(kcheb)) - ps(1)*exp(ima*gs(1))
val   = imag(valz)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The versions of these routines which take g'(x) as a input paramter rather than
!  computing it.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine levin_adapexp_diff11(vars,ier,eps,a,b,fun,par1,par2,par3,val,nints)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun2)           :: fun
double complex                  :: val
!
!  Evaluate an integral of the form (1) given an external subroutine for
!  evaluating the functions f(x), g(x) and g'(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) nad conforms to the levin_fun interface
!    par? - parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double complex, allocatable      :: gs0(:), fs0(:), gders0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack

nints = 0

allocate(stack(4,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb), gders0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)
epssq      = eps**2

ts0 = (b-a)/2 * vars%xscheb + (b+a)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_expdiff(vars,a,b,kcheb,ts0,fs0,gs0,gders0,valr,vali)


nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = valr
stack(4,1) = vali

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
valr0  = stack(3,nstack)
vali0  = stack(4,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = (c0-a0)/2 * vars%xscheb + (c0+a0)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_expdiff(vars,a0,c0,kcheb,ts0,fs0,gs0,gders0,valr1,vali1)

ts0    = (b0-c0)/2 * vars%xscheb + (b0+c0)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_expdiff(vars,c0,b0,kcheb,ts0,fs0,gs0,gders0,valr2,vali2)

dnorm =  max(  (valr0-valr1-valr2)**2 , (vali0-vali1-vali2)**2)
if (dnorm .lt. epssq ) then
val = val + valr0+ima*vali0
nints = nints+1
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr1
stack(4,nstack) = vali1

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = valr2
stack(4,nstack) = vali2

endif

end do

end subroutine


subroutine levin_adapexp_diff1(vars,ier,eps,a,b,fun,par1,par2,par3,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun2)           :: fun
double complex                  :: val
!
!  Evaluate an integral of the form (1) given an external subroutine for
!  evaluating the functions f(x), g(x) and g'(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) nad conforms to the levin_fun interface
!    par? - parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double complex, allocatable      :: gs0(:), fs0(:), gders0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(4,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb), gders0(kcheb))

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)
epssq      = eps**2

ts0 = (b-a)/2 * vars%xscheb + (b+a)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_expdiff(vars,a,b,kcheb,ts0,fs0,gs0,gders0,valr,vali)


nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = valr
stack(4,1) = vali

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
valr0  = stack(3,nstack)
vali0  = stack(4,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = (c0-a0)/2 * vars%xscheb + (c0+a0)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)

call levin_oneint_expdiff(vars,a0,c0,kcheb,ts0,fs0,gs0,gders0,valr1,vali1)

ts0    = (b0-c0)/2 * vars%xscheb + (b0+c0)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_expdiff(vars,c0,b0,kcheb,ts0,fs0,gs0,gders0,valr2,vali2)

dnorm = (valr0-valr1-valr2)**2 + (vali0-vali1-vali2)**2


if (dnorm .lt. epssq) then
val = val + valr0+ima*vali0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr1
stack(4,nstack) = vali1

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = valr2
stack(4,nstack) = vali2

endif

end do

end subroutine




subroutine levin_oneint_expdiff(vars,a,b,kcheb,ts,fs,gs,gders,valr,vali)
implicit double precision (a-h,o-z)
type(levin_vars_t)        :: vars
double complex            :: gs(kcheb), fs(kcheb), gders(kcheb)
double complex            :: amatr(kcheb,kcheb), ps(kcheb)
double precision          :: ts(kcheb)
double complex            :: val

amatr = 2.0d0/(b-a)*vars%adiff

do i=1,kcheb
amatr(i,i) = amatr(i,i) + gders(i)
end do

! call linalg0_solve(kcheb,amatr,fs)

! eps = epsilon(0.0d0)*norm2(abs(amatr))*10
eps = epsilon(0.0d0)*10*min(norm2(abs(amatr)),1.0d0)
call linalg0_qrtsolve(eps,kcheb,kcheb,amatr,ps,fs)

! eps = 1.0d-7 * norm2(abs(amatr))
! call linalg0_svdsolve_c1(ier,eps,kcheb,kcheb,amatr,ps,fs)

 
val  = ps(kcheb)*exp(gs(kcheb)) - ps(1)*exp(gs(1))
valr = real(val)
vali = imag(val)

end subroutine


subroutine levin_adapcos_diff1(vars,ier,eps,a,b,fun,par1,par2,par3,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun4)           :: fun
double precision                :: val
!
!  Evaluate an integral of the form (2) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) 
!    par? - parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double precision, allocatable    :: gs0(:), fs0(:), gders0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(3,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb), gders0(kcheb) )

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)

ts0 = (b-a)/2 * vars%xscheb + (b+a)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_cosdiff(vars,a,b,kcheb,ts0,fs0,gs0,gders0,val0)

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = val0

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
val0   = stack(3,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = (c0-a0)/2 * vars%xscheb + (c0+a0)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_cosdiff(vars,a0,c0,kcheb,ts0,fs0,gs0,gders0,valr)

ts0    = (b0-c0)/2 * vars%xscheb + (b0+c0)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_cosdiff(vars,c0,b0,kcheb,ts0,fs0,gs0,gders0,vall)


dnorm = abs(val0-valr-vall)


if (dnorm .lt. eps) then
val = val + val0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = vall

endif

end do

end subroutine



subroutine levin_oneint_cosdiff(vars,a,b,kcheb,ts,fs,gs,gders,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)        :: vars
double precision          :: gs(kcheb), fs(kcheb), gders(kcheb)
double complex            :: ps(kcheb), amatr(kcheb,kcheb), rhs(kcheb)
double precision          :: ts(kcheb)
double complex            :: ima, valz

ima = (0.0d0,1.0d0)

amatr = 2.0d0/(b-a)*vars%adiff
!gders = ima*2.0d0/(b-a)*matmul(vars%adiff,gs)

do i=1,kcheb
amatr(i,i) = amatr(i,i) + ima*gders(i)
end do

rhs = fs
call linalg0_solve(kcheb,amatr,rhs)

valz  = rhs(kcheb)*exp(ima*gs(kcheb)) - rhs(1)*exp(ima*gs(1))
val   = real(valz)

end subroutine



subroutine levin_adapsin_diff1(vars,ier,eps,a,b,fun,par1,par2,par3,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)              :: vars
procedure(levin_fun4)           :: fun
double precision                :: val
!
!  Evaluate an integral of the form (2) given an external subroutine for
!  evaluating the functions f(x) and g(x).
!
!  Input parameters:
!    vars - the structure populated by the initialization code
!    eps - the desired precision for the calculations
!    fun - an external subroutine which provides the values of the functions g(x)
!      and g(x) 
!    par? - parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0  indicates successful execution
!      ier = 4  means that the maximum number of intervals was exceeded while
!        attempting to compute the integral
!
!    val - the value of the integral as calculated by this routine
!

double precision, allocatable    :: stack(:,:)
double precision, allocatable    :: ts0(:)
double precision, allocatable    :: gs0(:), fs0(:), gders0(:)
double complex                   :: ima

kcheb    = vars%kcheb
maxstack = vars%maxstack


allocate(stack(3,maxstack) )
allocate(ts0(kcheb), gs0(kcheb), fs0(kcheb), gders0(kcheb) )

ier        = 0
val        = 0
ind        = 0
ima        = (0.0d0,1.0d0)

ts0 = (b-a)/2 * vars%xscheb + (b+a)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_sindiff(vars,a,b,kcheb,ts0,fs0,gs0,gders0,val0)

nstack     = 1
stack(1,1) = a
stack(2,1) = b
stack(3,1) = val0

do while (nstack > 0)

a0     = stack(1,nstack)
b0     = stack(2,nstack)
val0   = stack(3,nstack)
nstack = nstack-1


c0     = (a0+b0)/2
ts0    = (c0-a0)/2 * vars%xscheb + (c0+a0)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_sindiff(vars,a0,c0,kcheb,ts0,fs0,gs0,gders0,valr)

ts0    = (b0-c0)/2 * vars%xscheb + (b0+c0)/2
call fun(kcheb,ts0,gs0,gders0,fs0,par1,par2,par3)
call levin_oneint_sindiff(vars,c0,b0,kcheb,ts0,fs0,gs0,gders0,vall)


dnorm = abs(val0-valr-vall)


if (dnorm .lt. eps) then
val = val + val0
else

if (nstack+2 .gt. maxstack) then
ier = 4
return
endif

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = c0
stack(3,nstack) = valr

nstack          = nstack+1
stack(1,nstack) = c0
stack(2,nstack) = b0
stack(3,nstack) = vall

endif

end do

end subroutine



subroutine levin_oneint_sindiff(vars,a,b,kcheb,ts,fs,gs,gders,val)
implicit double precision (a-h,o-z)
type(levin_vars_t)        :: vars
double precision          :: gs(kcheb), fs(kcheb), gders(kcheb)
double complex            :: ps(kcheb), amatr(kcheb,kcheb), rhs(kcheb)
double precision          :: ts(kcheb)
double complex            :: ima, valz

ima = (0.0d0,1.0d0)

amatr = 2.0d0/(b-a)*vars%adiff
!gders = ima*2.0d0/(b-a)*matmul(vars%adiff,gs)

do i=1,kcheb
amatr(i,i) = amatr(i,i) + ima*gders(i)
end do

rhs = fs
call linalg0_solve(kcheb,amatr,rhs)

valz  = rhs(kcheb)*exp(ima*gs(kcheb)) - rhs(1)*exp(ima*gs(1))
val   = imag(valz)

end subroutine

end module 
