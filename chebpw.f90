!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing and manipulating piecewise Chebyshev
!  expansions of univariate functions.  A piecewise Chebyshev expansion is a sum
!  of the form
!
! 
!         n-1  k-1                            (     2            a_i + a_{i+1}  )
!  f(t)=  sum  sum b    chi  (t)           T  ( ----------- x +  -------------  )         (1)
!         i=0  j=0  ij       [a_i,a_{i+1})  j ( a_{i+1}-a_i      a_i - a_{i+1}  )
!
!  where
!
!    a = a_0 < a_1 < a_2 < ... < a_n = b
!
!  is a partition of the interval [a,b], chi_I denotes the indicator function
!  on the interval I and T_j is the Chebyshev polynomial of degree j.  We
!  refer to the n intervals 
!
!    (a_0,a_1), ..., (a_{n-1},a_n)                                                        (2)
!
!  as the set of discretization intervals for the expansion (1), and the integer
!  k as the order of the expansion (1).  Together k and the collection of intervals
!  (2) is known as a "piecewise Chebyshev discretization scheme" for the expansion (1).
!
!  This code represents expansions of the form (1) via two different mechanism: using
!  the vector
!
!    (   b0{0}     )
!    (   b0{1}     )
!    (   ...       )
!    (   b0{k-1}   )                                                                      (3)
!    (   b1{0}     )
!    (   ...       )
!    (   b1{k-1}   )
!    (    ...      )
!    ( b{n-1}{k-1} )
!
!  of expansion coefficients and with the vector 
!
!    ( f( x_00 )          )
!    ( f( x_01 )          )
!    (    ...             )
!    ( f( x_0k-1}  )      )
!    ( f( x_10     )      )                                                               (4)
!    (    ...             )
!    ( f( x_1{k-1} )      )
!    (    ...             )
!    ( f( x_{n-1}0 )      )
!         ...
!    ( f( x_{n-1}n{k-1} ) )
!
!  of values of the expansion (1) at the "practical" Chebyshev nodes on the
!  discretization intervals.  That is, in (4),
!
!    x_i0, ..., xi(k-1)
!
!  denote the practical Chebyshev nodes on the interval [a_i,a_{i+1}].  
!
!  Note that since the left and right practical Chebyshev nodes coincide with the
!  endpoints of the intervals, the vector (4) contains redudant values of f(x).
!  The author has generally found that eliminating this redudancy is not worth
!  the ensuing effort.
!
!  The discretization intervals are arranged in a "flat" lists of intervals.  In
!  certain circumstances, this can make the evaluation of expansions expensive.
!  See chebpw_tree.f90 for a code which arranges the discretization intervals in a
!  tree in order to accelerate the evaluation of expansions.
!
!  The following subroutines can be used to construct a piecewise discretization
!  scheme:
!
!    chebpw_uniform - construct a chebscheme structure describing a piecewise
!      Chebyshev discretization scheme with equispaced intervals
! 
!    chebpw_specified - construct a chebscheme structure describing a discretization
!      scheme with a specified set of discretization intervals
!
!    chebpw_adap - adaptively construct a discretization scheme which represents
!      one or more real-valued user-supplied functions specified via an external 
!      subroutine and return vectors representing the input functions
!
!    chebpw_refine - refine an existing discretization until it represents a
!      specified function or collection of functions to a specified relative
!      accuracy
!
!    chebpw_merge - merge two existing discretization schemes given over the
!      same interval into a single discretization scheme
!
!  The following routines return information about a piecewise discretization
!  scheme:
!
!    chebpw_interval - return the extents of the interval on which the
!     discretization scheme is given
!
!    chebpw_ab - return the list of discretization intervals associated with
!      a piecewise discretization scheme
!
!    chebpw_quad - return the Clenshaw-Curtis quadrature rule associated with
!       a piecewise discretization scheme
!
!  The following routines manipulate expansions:
!
!    (*) chebpw_plot - produce a python script which generates a PDF file
!      displaying a plot of one or more expansions of the form (1)
!
!    chebpw_coefs - given the vector of values (4) of an expansion of the form (1),
!      compute the vector of expansion coefficients (3)
!
!    chebpw_roots - use the "colleague" matrix method to find the roots of a
!      real-valued expansion on the interval over which it is defined
!
!    chebpw_diff - given the vector (4) of values of an expansion of the form (1),
!      compute the vector of values (4) of the derivative of that expansion
!
!    chebpw_int - given the vector (4) of values of an expansion f(t) of the form (1)
!      and a point c in the interval [a,b],  compute the values of the function
!
!            t
!        \int   f(x) dx
!            c 
!
!     at the discretization nodes
!
!    chebpw_intl - given the vector (4) of values of an expansion of the form (1),
!      return the vector of values of an expansion of the function
!
!              t
!          \int  f(s) ds
!              a
!
!  The following routines allow for the evaluation of expansions:
!
!    chebpw_eval - evaluate one or more expansions of the form (1) given
!      the corresponding vectors of coefficients (3)
!
!    chebpw_evalder - evaluate one or more expansions of the form (1) and their
!      derivatives at one or more points given the corresponding vectors of 
!      coefficients (3)
!
!    chebpw_interp - evaluate one or more expansions of the form (1) given
!      the corresponding vectors of values (4)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Observations / To do:
!
!    - Write "adapdisc" routines which take an initial list of intervals (not a 
!      chebscheme structure) as input
!
!    - Write a chebpw_combine routine which combines a discretization on the interval
!      [a,c] with one on the interval [c,b] to form one on the interval [a,b]
!
!    - Further develop the plotting routines ... in particular, plot over a smaller
!      interval and sample the points more judiciously
!
!    - Write a routine which copies the data (nodes, weights, spectral matrices)
!      from an existing structure but then has a new set of intervals ... this will
!      speed up certain codes
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module chebpw

use utils
use chebyshev
use iso_c_binding

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The chebpw_scheme structure describes a piecewise Chebyshev discretization scheme
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type            chebpw_scheme

double precision                    :: a,b              ! endpoints of the global interval

! data for the local Chebyshev expansions
integer                             :: k                ! order of the local expansion on each interval

double precision, allocatable       :: xscheb(:)        ! the k Chebyshev nodes
double precision, allocatable       :: whtscheb(:)      ! the k Chebyshev weights
double precision, allocatable       :: acoefs(:,:)      ! (k,k) vals-to-coefs matrix

! The following are initialized as needed
double precision, allocatable       :: adiff(:,:)       ! spectral differnetiation matrix
double precision, allocatable       :: aintl(:,:)       ! left spectral integration matrix
double precision, allocatable       :: aintr(:,:)       ! right spectral integration matrix

! the list of discretization intervals
integer                             :: nints   
double precision, allocatable       :: ab(:,:) 

end type        chebpw_scheme


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The chebpw_params structure is used to pass optional parameters to the rouintes for
!  constructing discretization schemes; see chebpw_adap for more information
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type            chebpw_params
integer                             :: k
integer                             :: ntest
integer                             :: maxints
double precision                    :: epsint
end type        chebpw_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The chebw_plotopts structure is used to specifying various optional parameters to
!  the plotting routines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type            chebpw_plotopts

end type        chebpw_plotopts

!
!  Various interfaces for routines which take external subroutines
!  as input
!

interface       

subroutine chebpw_fun1(n,ts,vals,par1,par2,par3)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                             :: n
double precision                    :: ts(:)
double precision                    :: vals(:)
end subroutine

subroutine chebpw_fun2(n,ts,vals1,vals2,par1,par2,par3)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                             :: n
double precision                    :: ts(:)
double precision                    :: vals1(:)
double precision                    :: vals2(:)
end subroutine


subroutine chebpw_fun3(n,ts,nfuns,vals,par1,par2,par3)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                             :: n
double precision                    :: ts(:)
double precision                    :: vals(:,:)
end subroutine


subroutine chebpw_fun4(n,ts,vals,par1,par2,par3)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                             :: n
double precision                    :: ts(:)
double complex                      :: vals(:)
end subroutine

subroutine chebpw_fun5(n,ts,vals1,vals2,par1,par2,par3)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                             :: n
double precision                    :: ts(:)
double complex                      :: vals1(:), vals2(:)
end subroutine

subroutine chebpw_fun6(n,ts,nfuns,vals,par1,par2,par3)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                             :: n
double precision                    :: ts(:)
double complex                      :: vals(:,:)
end subroutine


subroutine chebpw_fun7(n,ts,vals,userptr)
use iso_c_binding
implicit double precision (a-h,o-z)
type(c_ptr)                         :: userptr
integer                             :: n
double precision                    :: ts(:)
double precision                    :: vals(:)
end subroutine


subroutine chebpw_fun8(n,ts,vals1,vals2,vals3,vals4,par1,par2,par3)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                             :: n
double precision                    :: ts(:)
double precision                    :: vals1(:)
double precision                    :: vals2(:)
double precision                    :: vals3(:)
double precision                    :: vals4(:)
end subroutine

subroutine chebpw_fun9(n,ts,vals1,vals2,vals3,par1,par2,par3)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                             :: n
double precision                    :: ts(:)
double complex                      :: vals1(:)
double complex                      :: vals2(:)
double complex                      :: vals3(:)
end subroutine

subroutine chebpw_fun10(n,ts,vals1,vals2,vals3,par1,par2,par3)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                             :: n
double precision                    :: ts(:)
double precision                    :: vals1(:)
double precision                    :: vals2(:)
double precision                    :: vals3(:)
end subroutine

subroutine chebpw_opfun1(n,ts,valsin,valsout,par1,par2,par3)
implicit double precision (a-h,o-z)
double precision                   :: ts(:), valsin(:), valsout(:)
end subroutine


end interface   

interface         chebpw_adap
module procedure  chebpw_adap1
module procedure  chebpw_adap2
module procedure  chebpw_adap3
module procedure  chebpw_adap4
module procedure  chebpw_adap5
module procedure  chebpw_adap6
module procedure  chebpw_adap7
module procedure  chebpw_adap8
module procedure  chebpw_adap9
module procedure  chebpw_adap10
end interface     chebpw_adap

interface         chebpw_refine
module procedure  chebpw_refine1
module procedure  chebpw_refine2
module procedure  chebpw_refine4
module procedure  chebpw_refine5
end interface     chebpw_refine

interface         chebpw_fit
module procedure  chebpw_fit1
module procedure  chebpw_fit2
end interface     chebpw_fit

interface         chebpw_plot
module procedure  chebpw_plot1
module procedure  chebpw_plot2
module procedure  chebpw_plot3
module procedure  chebpw_plot4
module procedure  chebpw_plot5
end interface     chebpw_plot


interface         chebpw_coefs
module procedure  chebpw_coefs1
module procedure  chebpw_coefs2
module procedure  chebpw_coefs3
module procedure  chebpw_coefs4
end interface     chebpw_coefs

interface         chebpw_diff
module procedure  chebpw_diff1
module procedure  chebpw_diff2
module procedure  chebpw_diff3
module procedure  chebpw_diff4
end interface     chebpw_diff

interface         chebpw_intl
module procedure  chebpw_intl1
module procedure  chebpw_intl2
module procedure  chebpw_intl3
module procedure  chebpw_intl4
end interface     chebpw_intl

interface         chebpw_int
module procedure  chebpw_int1
module procedure  chebpw_int2
end interface     chebpw_int

interface        chebpw_eval
module procedure chebpw_eval1
module procedure chebpw_eval2
module procedure chebpw_eval3
module procedure chebpw_eval4
module procedure chebpw_eval5
module procedure chebpw_eval6
module procedure chebpw_eval101
module procedure chebpw_eval102
module procedure chebpw_eval103
module procedure chebpw_eval104
module procedure chebpw_eval105
module procedure chebpw_eval106
end interface    

interface        chebpw_evalder
module procedure chebpw_evalder1
module procedure chebpw_evalder2
module procedure chebpw_evalder3
module procedure chebpw_evalder4
module procedure chebpw_evalder5
module procedure chebpw_evalder6
module procedure chebpw_evalder101
module procedure chebpw_evalder102
module procedure chebpw_evalder103
module procedure chebpw_evalder104
module procedure chebpw_evalder105
module procedure chebpw_evalder106
end interface    

interface        chebpw_interp
module procedure chebpw_interp1
module procedure chebpw_interp2
module procedure chebpw_interp3
module procedure chebpw_interp4
module procedure chebpw_interp5
module procedure chebpw_interp6
module procedure chebpw_interp7
module procedure chebpw_interp8
module procedure chebpw_interp101
module procedure chebpw_interp102
module procedure chebpw_interp103
module procedure chebpw_interp104
module procedure chebpw_interp105
module procedure chebpw_interp106
module procedure chebpw_interp107
module procedure chebpw_interp108
end interface    

contains


subroutine chebpw_uniform(chebscheme,a,b,nints,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
type(chebpw_params), optional              :: opts

!
!  Construct a chebscheme structure describing an equispaced piecewise
!  Chebyshev discretization scheme.
!
!  Input parameters;
!    (a,b) - the interval over which the expansions will be given
!    k - the number terms in the local Chebyshev expansions
!    nints - the number of equispaced intervals into which the interval should be
!      decomposed
!
!  Output parameters:
!    chebscheme - the structure describing the desired scheme
!

!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif


call chebpw_chebdata(chebscheme,k)

chebscheme%a        = a
chebscheme%b        = b
chebscheme%k        = k
chebscheme%nints    = nints

dlen = (b-a)/(nints+0.0d0)

allocate(chebscheme%ab(2,nints))
do int=1,nints
chebscheme%ab(1,int) = a + dlen*(int-1)
chebscheme%ab(2,int) = a + dlen*int
end do

end subroutine


subroutine chebpw_specified(chebscheme,k,nints,ab)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
double precision                           :: ab(:,:)
!
!  Construct a chebscheme structure describing a scheme with a specified
!  list of intervals.
!
!  Input parameters;
!    (a,b) - the interval over which the expansions will be given
!    k - the number terms in the local Chebyshev expansions
!    nints - the number of intervals
!    ab - the list of intervals
!
!  Output parameters:
!    chebscheme - the structure describing the desired scheme
!

call chebpw_chebdata(chebscheme,k)


chebscheme%a        = ab(1,1)
chebscheme%b        = ab(2,nints)
chebscheme%k        = k
chebscheme%nints    = nints

dlen = (b-a)/(nints+0.0d0)

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

end subroutine


subroutine chebpw_merge(chebin1,chebin2,chebout)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                        :: chebin1
type(chebpw_scheme)                        :: chebin2
type(chebpw_scheme), intent(out)           :: chebout
!
!  Merge two input discretization schemes given over the same interval to procedure
!  one which suffices to represent functions well-represented by either of the two
!  input schemes.
!
!        IMPORTANT WARNING: this code does not check that the two 
!        schemes are given over the same interval; if they are not, 
!        the resulting scheme might be invalid.
!
!  Input parameters:
!    chebin1, chebin2 - the input schemes
!
!  Output parameters:
!    chebout - the merged output scheme
!
!

double precision, allocatable :: xs1(:),xs2(:), ab(:,:)


eps0    = epsilon(0.0d0)

nints1  = chebin1%nints
nints2  = chebin2%nints

k1      = chebin1%k
k2      = chebin2%k

k       = max(k1,k2)
maxints = chebin1%nints + chebin2%nints 


!
!  Form a list with the endpoints from both collections of intervals, 
!  sort it and remove duplicates
!
allocate(xs1(maxints+2))
allocate(xs2(maxints+2))
xs1  = 0
xs2  = 0
nxs1 = 0

do int=1,nints1
nxs1      = nxs1+1
xs1(nxs1) = chebin1%ab(1,int)
end do


nxs1      = nxs1+1
xs1(nxs1) = chebin1%ab(2,nints1)


do int=1,nints2
nxs1      = nxs1+1
xs1(nxs1) = chebin2%ab(1,int)
end do

nxs1      = nxs1+1
xs1(nxs1) = chebin2%ab(2,nints2)

call quicksort(nxs1,xs1)


!
!  Remove duplicates
!

nxs2 = 0
xs2  = 0

do i=1,nxs1-1

if( xs1(i+1)-xs1(i) .lt. eps0*10) cycle

nxs2      = nxs2+1
xs2(nxs2) = xs1(i)
end do

nxs2      = nxs2+1
xs2(nxs2) = xs1(nxs1)


nints = nxs2-1
allocate(ab(2,nints))
do int=1,nints
ab(1,int) = xs2(int)
ab(2,int) = xs2(int+1)
end do

call chebpw_specified(chebout,k,nints,ab)

end subroutine



subroutine chebpw_default_opts(k,ntest,maxints,epsint)
implicit double precision (a-h,o-z)
k        = 16
ntest    = 4
maxints  = 2**16
epsint   = 1.0d-7
end subroutine



subroutine chebpw_chebdata(chebscheme,k)
implicit double precision (a-h,o-z)
type(chebpw_scheme)         :: chebscheme

call chebyshev_quad(k,chebscheme%xscheb,chebscheme%whtscheb)
call chebyshev_coefsmatrix(k,chebscheme%acoefs)

if (allocated(chebscheme%adiff)) deallocate(chebscheme%adiff)
if (allocated(chebscheme%aintl)) deallocate(chebscheme%aintl)
if (allocated(chebscheme%aintr)) deallocate(chebscheme%aintr)

call chebyshev_diffmatrix(k,chebscheme%adiff)
call chebyshev_intlmatrix(k,chebscheme%aintl)
call chebyshev_intrmatrix(k,chebscheme%aintr)

end subroutine


subroutine chebpw_coefs0(chebscheme,k,ntest,acoefs0)
implicit double precision (a-h,o-z)
type(chebpw_scheme)           :: chebscheme
double precision, allocatable :: acoefs0(:,:)


allocate(acoefs0(ntest+1,k))
acoefs0(1,:)         = chebscheme%acoefs(1,:)
acoefs0(2:ntest+1,:) = chebscheme%acoefs(k-ntest+1:k,:)

end subroutine



subroutine chebpw_fit1(k,acoefs,ntest,acoefs0,vals,dcoefs)
implicit double precision (a-h,o-z)
double precision      :: acoefs(:,:), acoefs0(:,:)
double precision      :: vals(:)
double precision      :: coefs(k)



coefs  = matmul(acoefs,vals)
dd1    = norm2(abs(coefs))
dd2    = norm2(abs(coefs(k-ntest+1:k)))

! this is included to handle the case of a function which is equal to 0
if (dd1 .eq. 0) dd1 = dd1+1

dcoefs = dd2/dd1


end subroutine


subroutine chebpw_fit2(k,acoefs,ntest,acoefs0,vals,dcoefs)
implicit double precision (a-h,o-z)
double precision      :: acoefs(:,:), acoefs0(:,:)
double complex        :: vals(:), coefs(k)

coefs  = matmul(acoefs,vals)
dd1    = norm2(abs(coefs(1:k)))
dd2    = norm2(abs(coefs(k-ntest+1:k)))

! this is included to handle the case of a function which is equal to 0
if (dd1 .eq. 0) dd1 = dd1+1

dcoefs = dd2/dd1

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Adaptive discretization routines 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chebpw_adap1(ier,chebscheme,eps,a,b,fun,par1,par2,par3,valsout,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun1)                     :: fun
type(chebpw_params), optional              :: opts
double precision, allocatable, intent(out) :: valsout(:)
!
!  Attempt to construct a Chebyshev discretization scheme which represents a single
!  real-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout - the 
!
double precision, allocatable              :: ab(:,:), ab0(:,:)
double precision, allocatable              :: xs(:), vals(:,:), coefs(:), acoefs0(:,:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals(k,maxints+1),coefs(ntest+1))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif


! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,vals(:,nints+1),par1,par2,par3)


ifsplit = 0
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals(:,nints+1),dcoefs)
if (dcoefs .gt. eps) ifsplit = 1

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

allocate(valsout(nints*k))
valsout = reshape(vals,[k*nints])

end subroutine



subroutine chebpw_adap2(ier,chebscheme,eps,a,b,fun,par1,par2,par3,valsout1,valsout2,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun2)                    :: fun
type(chebpw_params), optional              :: opts
double precision, allocatable, intent(out) :: valsout1(:), valsout2(:)
!
!  Attempt to construct a Chebyshev discretization scheme which represents two
!  real-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout - the 
!
double precision, allocatable              :: ab(:,:), ab0(:,:)
double precision, allocatable              :: xs(:), vals1(:,:), vals2(:,:)
double precision, allocatable              :: coefs(:), acoefs0(:,:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0


!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals1(k,maxints+1),vals2(k,maxints+1),coefs(ntest+1))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif


! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,vals1(:,nints+1),vals2(:,nints+1),par1,par2,par3)

ifsplit = 0
dcoefs1 = 1d300
dcoefs2 = 1d300


call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals1(:,nints+1),dcoefs1)
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals2(:,nints+1),dcoefs2)

if (dcoefs1 .gt. eps) ifsplit = 1
if (dcoefs2 .gt. eps) ifsplit = 1


if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

allocate(valsout1(nints*k))
allocate(valsout2(nints*k))

valsout1 = reshape(vals1,[k*nints])
valsout2 = reshape(vals2,[k*nints])


end subroutine



subroutine chebpw_adap3(ier,chebscheme,eps,a,b,nfuns,fun,par1,par2,par3,valsout,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun3)                     :: fun
type(chebpw_params), optional              :: opts
double precision, allocatable, intent(out) :: valsout(:,:)
!
!  Attempt to construct a Chebyshev discretization scheme which represents multiple
!  real-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout - each column is the vector of values of one of the discretized functions
!
double precision, allocatable              :: ab(:,:), ab0(:,:)
double precision, allocatable              :: xs(:), vals(:,:), coefs(:,:), acoefs0(:,:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals(k,nfuns),coefs(ntest+1,nfuns))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif

! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,nfuns,vals,par1,par2,par3)

ifsplit = 0
do j=1,nfuns
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals(:,j),dcoefs)
if (dcoefs .gt. eps) then
ifsplit = 1
exit
endif
end do

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

allocate(valsout(nints*k,nfuns))
do int=1,nints
i1 = (int-1)*k+1
i2 = int*k

a0 = ab(1,int)
b0 = ab(2,int)
xs = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,nfuns,valsout(i1:i2,:),par1,par2,par3)

end do

end subroutine


subroutine chebpw_adap4(ier,chebscheme,eps,a,b,fun,par1,par2,par3,valsout,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun4)                     :: fun
type(chebpw_params), optional              :: opts
double complex, allocatable, intent(out)   :: valsout(:)
!
!  Attempt to construct a Chebyshev discretization scheme which represents a single
!  complex-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout - the 
!
double precision, allocatable              :: ab(:,:), ab0(:,:), xs(:), acoefs0(:,:)
double complex, allocatable                :: vals(:,:), coefs(:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals(k,maxints+1),coefs(ntest+1))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif


! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,vals(:,nints+1),par1,par2,par3)

ifsplit = 0

call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals(:,nints+1),dcoefs)
if (dcoefs .gt. eps) ifsplit = 1

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

allocate(valsout(nints*k))
valsout = reshape(vals,[k*nints])

end subroutine



subroutine chebpw_adap5(ier,chebscheme,eps,a,b,fun,par1,par2,par3,valsout1,valsout2,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun5)                     :: fun
type(chebpw_params), optional              :: opts
double complex, allocatable, intent(out)   :: valsout1(:), valsout2(:)
!
!  Attempt to construct a Chebyshev discretization scheme which represents two
!  real-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout - the 
!
double precision, allocatable              :: ab(:,:), ab0(:,:), xs(:), acoefs0(:,:)
double complex, allocatable                :: vals1(:,:), vals2(:,:), coefs(:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

!
!  Read the optional parameters, if present
!


if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals1(k,maxints+1),vals2(k,maxints+1),coefs(ntest+1))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif


! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,vals1(:,nints+1),vals2(:,nints+1),par1,par2,par3)

ifsplit = 0

call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals1(:,nints+1),dcoefs1)
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals2(:,nints+1),dcoefs2)

if (dcoefs1 .gt. eps) ifsplit = 1
if (dcoefs2 .gt. eps) ifsplit = 1

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

allocate(valsout1(nints*k))
allocate(valsout2(nints*k))

valsout1 = reshape(vals1,[k*nints])
valsout2 = reshape(vals2,[k*nints])

end subroutine



subroutine chebpw_adap6(ier,chebscheme,eps,a,b,nfuns,fun,par1,par2,par3,valsout,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun6)                     :: fun
type(chebpw_params), optional              :: opts
double complex, allocatable, intent(out)   :: valsout(:,:)
!
!  Attempt to construct a Chebyshev discretization scheme which represents multiple
!  complex-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout - each column is the vector of values of one of the discretized functions
!
double precision, allocatable              :: ab(:,:), ab0(:,:), xs(:), acoefs0(:,:)
double precision, allocatable              :: stack(:,:), stack0(:,:)
double complex, allocatable                :: vals(:,:), coefs(:,:)

ier = 0

!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals(k,nfuns),coefs(ntest+1,nfuns))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif

! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,nfuns,vals,par1,par2,par3)

! coefs = matmul(acoefs0,vals)

ifsplit = 0
do j=1,nfuns
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals(:,j),dcoefs)
if (dcoefs .gt. eps) then
ifsplit = 1
exit
endif
end do

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

allocate(valsout(nints*k,nfuns))
do int=1,nints
i1 = (int-1)*k+1
i2 = int*k

a0 = ab(1,int)
b0 = ab(2,int)
xs = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,nfuns,valsout(i1:i2,:),par1,par2,par3)

end do

end subroutine


subroutine chebpw_adap7(ier,chebscheme,eps,a,b,fun,userptr,valsout,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun7)                     :: fun
type(c_ptr)                                :: userptr
type(chebpw_params), optional              :: opts
double precision, allocatable, intent(out) :: valsout(:)
!
!  Attempt to construct a Chebyshev discretization scheme which represents a single
!  real-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout - the 
!
double precision, allocatable              :: ab(:,:), ab0(:,:)
double precision, allocatable              :: xs(:), vals(:,:), coefs(:), acoefs0(:,:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals(k,maxints+1),coefs(ntest+1))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif


! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,vals(:,nints+1),userptr)

ifsplit = 0

call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals(:,nints+1),dcoefs)
if (dcoefs .gt. eps) ifsplit = 1

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

allocate(valsout(nints*k))
valsout = reshape(vals,[k*nints])

end subroutine



subroutine chebpw_adap8(ier,chebscheme,eps,a,b,fun,par1,par2,par3,valsout1,valsout2,&
  valsout3,valsout4,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun8)                    :: fun
type(chebpw_params), optional              :: opts
double precision, allocatable, intent(out) :: valsout1(:), valsout2(:)
double precision, allocatable, intent(out) :: valsout3(:), valsout4(:)

!
!  Attempt to construct a Chebyshev discretization scheme which represents four
!  real-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout? - the vectors of values of the functions
!
double precision, allocatable              :: ab(:,:), ab0(:,:), xs(:)
double precision, allocatable              :: vals1(:,:), vals2(:,:)
double precision, allocatable              :: vals3(:,:), vals4(:,:)
double precision, allocatable              :: coefs(:), acoefs0(:,:), coefs0(:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals1(k,maxints+1),vals2(k,maxints+1),coefs(ntest+1))
allocate(vals3(k,maxints+1),vals4(k,maxints+1),coefs0(k))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif

ifsplit = 0

! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,vals1(:,nints+1),vals2(:,nints+1),vals3(:,nints+1),vals4(:,nints+1),&
  par1,par2,par3)


call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals1(:,nints+1),dcoefs1)
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals2(:,nints+1),dcoefs2)
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals3(:,nints+1),dcoefs3)
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals4(:,nints+1),dcoefs4)

if (dcoefs1 .gt. eps) ifsplit = 1
if (dcoefs2 .gt. eps) ifsplit = 1
if (dcoefs3 .gt. eps) ifsplit = 1
if (dcoefs4 .gt. eps) ifsplit = 1

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

allocate(valsout1(nints*k))
allocate(valsout2(nints*k))
allocate(valsout3(nints*k))
allocate(valsout4(nints*k))

valsout1 = reshape(vals1,[k*nints])
valsout2 = reshape(vals2,[k*nints])
valsout3 = reshape(vals3,[k*nints])
valsout4 = reshape(vals4,[k*nints])

end subroutine



subroutine chebpw_adap9(ier,chebscheme,eps,a,b,fun,par1,par2,par3,valsout1,valsout2,&
  valsout3,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun9)                     :: fun
type(chebpw_params), optional              :: opts
double complex, allocatable, intent(out)   :: valsout1(:), valsout2(:)
double complex, allocatable, intent(out)   :: valsout3(:)

!
!  Attempt to construct a Chebyshev discretization scheme which represents four
!  real-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout? - the vectors of values of the functions
!
double precision, allocatable              :: ab(:,:), ab0(:,:), xs(:), acoefs0(:,:)
double complex, allocatable                :: vals1(:,:), vals2(:,:)
double complex, allocatable                :: vals3(:,:)
double complex, allocatable                :: coefs(:), coefs0(:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals1(k,maxints+1),vals2(k,maxints+1),coefs(ntest+1))
allocate(vals3(k,maxints+1),coefs0(k))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif

ifsplit = 0

! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,vals1(:,nints+1),vals2(:,nints+1),vals3(:,nints+1),par1,par2,par3)


call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals1(:,nints+1),dcoefs1)
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals2(:,nints+1),dcoefs2)
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals3(:,nints+1),dcoefs3)

if (dcoefs1 .gt. eps) ifsplit = 1
if (dcoefs2 .gt. eps) ifsplit = 1
if (dcoefs3 .gt. eps) ifsplit = 1

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)


allocate(valsout1(nints*k))
allocate(valsout2(nints*k))
allocate(valsout3(nints*k))

valsout1 = reshape(vals1,[k*nints])
valsout2 = reshape(vals2,[k*nints])
valsout3 = reshape(vals3,[k*nints])

end subroutine



subroutine chebpw_adap10(ier,chebscheme,eps,a,b,fun,par1,par2,par3,valsout1,valsout2,&
  valsout3,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebscheme
procedure(chebpw_fun10)                    :: fun
type(chebpw_params), optional              :: opts
double precision, allocatable, intent(out) :: valsout1(:), valsout2(:)
double precision, allocatable, intent(out) :: valsout3(:)

!
!  Attempt to construct a Chebyshev discretization scheme which represents four
!  real-valued input functions with a desired relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    (a,b) - interval on which the input function(s) are given
!    fun - a user-specified external functions which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebscheme - a structure describing the Chebyshev discretization scheme
!    valsout? - the vectors of values of the functions
!
double precision, allocatable              :: ab(:,:), ab0(:,:), xs(:), acoefs0(:,:)
double precision, allocatable              :: vals1(:,:), vals2(:,:)
double precision, allocatable              :: vals3(:,:)
double precision, allocatable              :: coefs(:), coefs0(:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

!
!  Read the optional parameters, if present
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
call chebpw_chebdata(chebscheme,k)
call chebpw_coefs0(chebscheme,k,ntest,acoefs0)

!
!  Adaptively refine the discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals1(k,maxints+1),vals2(k,maxints+1),coefs(ntest+1))
allocate(vals3(k,maxints+1),coefs0(k))

nints              = 0
nstack             = 1
stack(1,1)         = a
stack(2,1)         = b


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif

ifsplit = 0

! Check to see if the interval should be split
xs     = (b0-a0)/2*chebscheme%xscheb + (b0+a0)/2
call fun(k,xs,vals1(:,nints+1),vals2(:,nints+1),vals3(:,nints+1),par1,par2,par3)


call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals1(:,nints+1),dcoefs1)
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals2(:,nints+1),dcoefs2)
call chebpw_fit(k,chebscheme%acoefs,ntest,acoefs0,vals3(:,nints+1),dcoefs3)

if (dcoefs1 .gt. eps) ifsplit = 1
if (dcoefs2 .gt. eps) ifsplit = 1
if (dcoefs3 .gt. eps) ifsplit = 1

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebscheme%a     = a
chebscheme%b     = b
chebscheme%k     = k
chebscheme%nints = nints

allocate(chebscheme%ab(2,nints))
chebscheme%ab = ab(:,1:nints)

allocate(valsout1(nints*k))
allocate(valsout2(nints*k))
allocate(valsout3(nints*k))

valsout1 = reshape(vals1,[k*nints])
valsout2 = reshape(vals2,[k*nints])
valsout3 = reshape(vals3,[k*nints])

end subroutine


subroutine chebpw_refine1(ier,eps,chebin,chebout,fun,par1,par2,par3,valsout,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebout
type(chebpw_scheme)                        :: chebin
procedure(chebpw_fun1)                     :: fun
type(chebpw_params), optional              :: opts
double precision, allocatable, intent(out) :: valsout(:)
!
!  Adaptively refine a existing discretization scheme so that it represents
!  a user-specified function to a specified relative accuracy.
!
!  NOTE: the parameter k is ignored
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    chebin - a structure describing an existing discretization scheme
!    fun - a user-specified external function which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebout - a structure describing the new discretization scheme
!    valsout - the values of the input function at the nodes of the new
!      discretization scheme
!
double precision, allocatable              :: ab(:,:), ab0(:,:)
double precision, allocatable              :: xs(:), vals(:,:), coefs(:), acoefs0(:,:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0


!
!  Read the optional parameters, if present .. but we override the k
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
chebout%k = chebin%k
k         = chebin%k

allocate( chebout%xscheb(k), chebout%whtscheb(k), chebout%acoefs(k,k) )
allocate( chebout%adiff(k,k), chebout%aintl(k,k), chebout%aintr(k,k) )

chebout%xscheb   = chebin%xscheb
chebout%whtscheb = chebin%whtscheb
chebout%acoefs   = chebin%acoefs
chebout%aintl    = chebin%aintl
chebout%aintr    = chebin%aintr
chebout%adiff    = chebin%adiff

call chebpw_coefs0(chebout,k,ntest,acoefs0)

!
!  Adaptively refine the existing discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals(k,maxints+1),coefs(ntest+1))

nints              = 0
nstack             = 0

do int=chebin%nints,1,-1
nstack             = nstack+1
stack(1,nstack)    = chebin%ab(1,int)
stack(2,nstack)    = chebin%ab(2,int)
end do


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif


! Check to see if the interval should be split
xs     = (b0-a0)/2*chebout%xscheb + (b0+a0)/2
call fun(k,xs,vals(:,nints+1),par1,par2,par3)


ifsplit = 0
call chebpw_fit(k,chebout%acoefs,ntest,acoefs0,vals(:,nints+1),dcoefs)
if (dcoefs .gt. eps) ifsplit = 1


if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebout%a     = a
chebout%b     = b
chebout%k     = k
chebout%nints = nints

allocate(chebout%ab(2,nints))
chebout%ab = ab(:,1:nints)

allocate(valsout(nints*k))
valsout = reshape(vals,[k*nints])

end subroutine


subroutine chebpw_refine2(ier,eps,chebin,chebout,fun,par1,par2,par3,valsout1,valsout2,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebout
type(chebpw_scheme)                        :: chebin
procedure(chebpw_fun2)                     :: fun
type(chebpw_params), optional              :: opts
double precision, allocatable, intent(out) :: valsout1(:), valsout2(:)
!
!  Adaptively refine a existing discretization scheme so that it represents
!  a user-specified function to a specified relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    chebin - a structure describing an existing discretization scheme
!    fun - a user-specified external function which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebout - a structure describing the new discretization scheme
!    valsout - the values of the input function at the nodes of the new
!      discretization scheme
!
double precision, allocatable              :: ab(:,:), ab0(:,:)
double precision, allocatable              :: xs(:), vals1(:,:), vals2(:,:), coefs(:), acoefs0(:,:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
chebout%k = chebin%k
k         = chebin%k

allocate( chebout%xscheb(k), chebout%whtscheb(k), chebout%acoefs(k,k) )
allocate( chebout%adiff(k,k), chebout%aintl(k,k), chebout%aintr(k,k) )

chebout%xscheb   = chebin%xscheb
chebout%whtscheb = chebin%whtscheb
chebout%acoefs   = chebin%acoefs
chebout%aintl    = chebin%aintl
chebout%aintr    = chebin%aintr
chebout%adiff    = chebin%adiff

call chebpw_coefs0(chebout,k,ntest,acoefs0)

!
!  Adaptively refine the existing discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals1(k,maxints+1),coefs(ntest+1))
allocate(vals2(k,maxints+1))

nints              = 0
nstack             = 0

do int=chebin%nints,1,-1
nstack             = nstack+1
stack(1,nstack)    = chebin%ab(1,int)
stack(2,nstack)    = chebin%ab(2,int)
end do


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1


if (b0-a0 .lt. epsint) then
ier = 8
return
endif

! Check to see if the interval should be split
xs     = (b0-a0)/2*chebout%xscheb + (b0+a0)/2
call fun(k,xs,vals1(:,nints+1),vals2(:,nints+1),par1,par2,par3)

ifsplit = 0
call chebpw_fit(k,chebout%acoefs,ntest,acoefs0,vals1(:,nints+1),dcoefs1)
call chebpw_fit(k,chebout%acoefs,ntest,acoefs0,vals2(:,nints+1),dcoefs2)

if (dcoefs1 .gt. eps) ifsplit = 1
if (dcoefs2 .gt. eps) ifsplit = 1


if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebout%a     = a
chebout%b     = b
chebout%k     = k
chebout%nints = nints

allocate(chebout%ab(2,nints))
chebout%ab = ab(:,1:nints)

allocate(valsout1(nints*k))
valsout1 = reshape(vals1,[k*nints])
valsout2 = reshape(vals2,[k*nints])

end subroutine





subroutine chebpw_refine4(ier,eps,chebin,chebout,fun,par1,par2,par3,valsout,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebout
type(chebpw_scheme)                        :: chebin
procedure(chebpw_fun4)                     :: fun
type(chebpw_params), optional              :: opts
double complex, allocatable, intent(out  ) :: valsout(:)
!
!  Adaptively refine a existing discretization scheme so that it represents
!  a user-specified function to a specified relative accuracy.
!
!  NOTE: the parameter k is ignored
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    chebin - a structure describing an existing discretization scheme
!    fun - a user-specified external function which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!    opts - an optional structure specifying several algorithm parameters; its
!      entries are
!
!      k - number of Chebyshev nodes per interval 
!      ntest - the number of coefficients to examine when determining if
!        a local expansion is sufficient
!      maxints - the maximum number of intervals  
!      epsint - the length of the smallest allowable interval
!
!      SEE THE CHEBPW_DEFAULT_OPTS SUBROUTINE FOR THE DEFAULT VALUES
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebout - a structure describing the new discretization scheme
!    valsout - the values of the input function at the nodes of the new
!      discretization scheme
!
double precision, allocatable              :: ab(:,:), ab0(:,:), acoefs0(:,:), xs(:)
double complex, allocatable                :: vals(:,:), coefs(:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0


!
!  Read the optional parameters, if present .. but we override the k
!

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
chebout%k = chebin%k
k         = chebin%k

allocate( chebout%xscheb(k), chebout%whtscheb(k), chebout%acoefs(k,k) )
allocate( chebout%adiff(k,k), chebout%aintl(k,k), chebout%aintr(k,k) )

chebout%xscheb   = chebin%xscheb
chebout%whtscheb = chebin%whtscheb
chebout%acoefs   = chebin%acoefs
chebout%aintl    = chebin%aintl
chebout%aintr    = chebin%aintr
chebout%adiff    = chebin%adiff

call chebpw_coefs0(chebout,k,ntest,acoefs0)

!
!  Adaptively refine the existing discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals(k,maxints+1),coefs(ntest+1))

nints              = 0
nstack             = 0

do int=chebin%nints,1,-1
nstack             = nstack+1
stack(1,nstack)    = chebin%ab(1,int)
stack(2,nstack)    = chebin%ab(2,int)
end do


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif


! Check to see if the interval should be split
xs     = (b0-a0)/2*chebout%xscheb + (b0+a0)/2
call fun(k,xs,vals(:,nints+1),par1,par2,par3)


ifsplit = 0
call chebpw_fit(k,chebout%acoefs,ntest,acoefs0,vals(:,nints+1),dcoefs)
if (dcoefs .gt. eps) ifsplit = 1


if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebout%a     = a
chebout%b     = b
chebout%k     = k
chebout%nints = nints

allocate(chebout%ab(2,nints))
chebout%ab = ab(:,1:nints)

allocate(valsout(nints*k))
valsout = reshape(vals,[k*nints])

end subroutine


subroutine chebpw_refine5(ier,eps,chebin,chebout,fun,par1,par2,par3,valsout1,valsout2,opts)
implicit double precision (a-h,o-z)
type(chebpw_scheme), intent(out)           :: chebout
type(chebpw_scheme)                        :: chebin
procedure(chebpw_fun5)                     :: fun
type(chebpw_params), optional              :: opts
double complex, allocatable, intent(out)   :: valsout1(:), valsout2(:)
!
!  Adaptively refine a existing discretization scheme so that it represents
!  a user-specified function to a specified relative accuracy.
!
!  Input parameters:
!    eps - desired relative precision for the discretization
!    chebin - a structure describing an existing discretization scheme
!    fun - a user-specified external function which supplies the values of the
!      input functions
!    par?- user-supplied parameters which are passed to the external subroutine
!
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   means that the maximum number of intervals was exceeded
!       ier = 8   means that an interval smaller than the minimum
!                 allowable length was required
!       ier = 16  means that the internal stack used by this subroutine
!                 overflowed
!
!    chebout - a structure describing the new discretization scheme
!    valsout - the values of the input function at the nodes of the new
!      discretization scheme
!
double precision, allocatable              :: ab(:,:), ab0(:,:), xs(:), acoefs0(:,:)
double complex, allocatable                :: vals1(:,:), vals2(:,:), coefs(:)
double precision, allocatable              :: stack(:,:), stack0(:,:)

ier = 0

if ( present(opts) ) then
k         = opts%k
ntest     = opts%ntest
maxints   = opts%maxints
epsint    = opts%epsint
else
call chebpw_default_opts(k,ntest,maxints,epsint)
endif

!
!  Construct the Curtis-Clenshaw quadrature and the Chebyshev matrices
!
chebout%k = chebin%k
k         = chebin%k

allocate( chebout%xscheb(k), chebout%whtscheb(k), chebout%acoefs(k,k) )
allocate( chebout%adiff(k,k), chebout%aintl(k,k), chebout%aintr(k,k) )

ntest            = 4
chebout%xscheb   = chebin%xscheb
chebout%whtscheb = chebin%whtscheb
chebout%acoefs   = chebin%acoefs
chebout%aintl    = chebin%aintl
chebout%aintr    = chebin%aintr
chebout%adiff    = chebin%adiff

call chebpw_coefs0(chebout,k,ntest,acoefs0)

!
!  Adaptively refine the existing discretization
!

allocate(stack(2,maxints))
allocate(ab(2,maxints))

allocate(xs(k),vals1(k,maxints+1),coefs(ntest+1))
allocate(vals2(k,maxints+1))

nints              = 0
nstack             = 0

do int=chebin%nints,1,-1
nstack             = nstack+1
stack(1,nstack)    = chebin%ab(1,int)
stack(2,nstack)    = chebin%ab(2,int)
end do


do while (nstack > 0)

a0 = stack(1,nstack)
b0 = stack(2,nstack)
nstack=nstack-1

if (b0-a0 .lt. epsint) then
ier = 8
return
endif


! Check to see if the interval should be split
xs     = (b0-a0)/2*chebout%xscheb + (b0+a0)/2
call fun(k,xs,vals1(:,nints+1),vals2(:,nints+1),par1,par2,par3)

ifsplit = 0

call chebpw_fit(k,chebout%acoefs,ntest,acoefs0,vals1(:,nints+1),dcoefs1)
call chebpw_fit(k,chebout%acoefs,ntest,acoefs0,vals2(:,nints+1),dcoefs2)

if (dcoefs1 .gt. eps) ifsplit = 1
if (dcoefs2 .gt. eps) ifsplit = 1


if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxints) then
ier = 16
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints+1 .ge. maxints) then
ier = 4
return
endif

nints        = nints+1
ab(1,nints)  = a0
ab(2,nints)  = b0

endif

end do

chebout%a     = a
chebout%b     = b
chebout%k     = k
chebout%nints = nints

allocate(chebout%ab(2,nints))
chebout%ab = ab(:,1:nints)

allocate(valsout1(nints*k))
valsout1 = reshape(vals1,[k*nints])
valsout2 = reshape(vals2,[k*nints])

end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Information routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chebpw_info(chebscheme,k,nints)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                        :: chebscheme
!
!  Return some information about a chebyshev discretization scheme.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!
!  Output parameters:
!    k - the number of 
!    nints - the number of intervals
!

k         = chebscheme%k
nints     = chebscheme%nints

end subroutine


subroutine chebpw_quad(chebscheme,n,xs,whts)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                        :: chebscheme
double precision, allocatable, intent(out) :: xs(:), whts(:)
!
!  Return the quadrature obtained by combining the Curtis-Clenshaw on each discretization
!  interval.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!
!  Output parameters:
!    n - the number of nodes in the quadrature rule
!    xs, whts - the nodes and weights of the rule
!

k     = chebscheme%k
nints = chebscheme%nints

n     = nints*k

allocate(xs(n),whts(n))

do i=1,nints

a  = chebscheme%ab(1,i)
b  = chebscheme%ab(2,i)

i1 = (i-1)*k+1
i2 = i*k

xs(i1:i2)   = (b-a)/2 * chebscheme%xscheb + (b+a)/2
whts(i1:i2) = (b-a)/2 * chebscheme%whtscheb 
end do

end subroutine


subroutine chebpw_ab(chebscheme,nints,ab)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                        :: chebscheme
double precision, allocatable, intent(out) :: ab(:,:)
!
!  Return the list of discretization intervals.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!
!  Output parameters:
!    nints - the number of discretization intervals
!    ab - a (2,nints) array giving the list of discretization intervals
!

nints = chebscheme%nints

allocate(ab(2,nints))
ab = chebscheme%ab

end subroutine

subroutine chebpw_interval(chebscheme,a,b)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                        :: chebscheme
!
!  Return the extent of the interval on which the discretization scheme is given.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!
!  Output parameters:
!    a - the left-hand endpoint of the interval
!    b - the right-hand endpoint of the interval

!

a = chebscheme%ab(1,1)
b = chebscheme%ab(2,chebscheme%nints)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Routines for manipulating expansions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chebpw_coefs1(chebscheme,vals,coefs)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double precision                                  :: vals(:)
double precision, allocatable, intent(out)        :: coefs(:)
!
!  Compute the vector (3) of coefficients given the vector (4) of values.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4)
!
!  Output parameters:
!    coefs - the vector of coefficient (4)
!
nints = chebscheme%nints
k     = chebscheme%k

allocate(coefs(k*nints))

i1    = 1
i2    = k

do int=1,nints
coefs(i1:i2) = matmul(chebscheme%acoefs,vals(i1:i2))
i1 = i1+k
i2 = i2+k
end do

end subroutine

subroutine chebpw_coefs2(chebscheme,vals,coefs)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double complex                                    :: vals(:)
double complex, allocatable, intent(out)          :: coefs(:)
!
!  Compute the vector (3) of coefficients given the vector (4) of values.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4)
!
!  Output parameters:
!    coefs - the vector of coefficient (4)
!

nints = chebscheme%nints
k     = chebscheme%k

allocate(coefs(k*nints))

i1    = 1
i2    = k

do int=1,nints
coefs(i1:i2) = matmul(chebscheme%acoefs,vals(i1:i2))
i1 = i1+k
i2 = i2+k
end do

end subroutine


subroutine chebpw_coefs3(chebscheme,nfuns,vals,coefs)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double precision                                  :: vals(:,:)
double precision, allocatable, intent(out)        :: coefs(:,:)
!
!  Compute the vector (3) of coefficients given the vector (4) of values.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4)
!
!  Output parameters:
!    coefs - the vector of coefficient (4)
!

nints = chebscheme%nints
k     = chebscheme%k

allocate(coefs(k*nints,nfuns))

i1    = 1
i2    = k

do int=1,nints
coefs(i1:i2,:) = matmul(chebscheme%acoefs,vals(i1:i2,:))
i1 = i1+k
i2 = i2+k
end do

end subroutine


subroutine chebpw_coefs4(chebscheme,nfuns,vals,coefs)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double complex                                    :: vals(:,:)
double complex, allocatable, intent(out)          :: coefs(:,:)
!
!  Compute the vector (3) of coefficients given the vector (4) of values.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4)
!
!  Output parameters:
!    coefs - the vector of coefficient (4)
!

nints = chebscheme%nints
k     = chebscheme%k

allocate(coefs(nints*k,nfuns))

i1    = 1
i2    = k

do int=1,nints
coefs(i1:i2,:) = matmul(chebscheme%acoefs,vals(i1:i2,:))
i1 = i1+k
i2 = i2+k
end do

end subroutine


subroutine chebpw_int1(chebscheme,vals,c,valc,valsint)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double precision                                  :: vals(:)
double precision, allocatable, intent(out)        :: valsint(:)
!
!  Compute the values of the function
!
!                     t
!    g(t) = valc + int    f(s) ds
!                     c
!
!  with c a point in [a,b].
!  
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    valsint - the vector of values (4) of the antiderivative of the expansion
!     which is 0 at the left-hand endpoint of [a,b]
!

double precision     ::    coefs0(chebscheme%k), vals0(chebscheme%k)
double precision     ::    polsc(chebscheme%k+1)


nints = chebscheme%nints
k     = chebscheme%k

allocate(valsint(k*nints))

if (.not. allocated(chebscheme%aintl)) then
call chebyshev_intlmatrix(chebscheme%k,chebscheme%aintl)
endif

if (.not. allocated(chebscheme%aintr)) then
call chebyshev_intrmatrix(chebscheme%k,chebscheme%aintr)
endif

!
!  Find the interval containing the point c 
!

do intc=1,nints-1
if (c .lt. chebscheme%ab(2,intc)) exit
end do

! ! Handle the special cases where intc==1 or intc==nints

! if (intc .eq. 1) then

! val0 = valc
! do int=1,nints
! a0 = chebscheme%ab(1,int)
! b0 = chebscheme%ab(2,int)
! i1 = (int-1)*k + 1
! i2 = int*k
! valsint(i1:i2) = val0 + (b0-a0)/2 * matmul(chebscheme%aintl,vals(i1:i2))
! val0           = valsint(i2)
! end do
! return
! endif


! if (intc .eq. nints) then
! val0 = valc
! do int=nints,1,-1
! a0 = chebscheme%ab(1,int)
! b0 = chebscheme%ab(2,int)
! i1 = (int-1)*k + 1
! i2 = int*k
! valsint(i1:i2) = val0 + (b0-a0)/2 * matmul(chebscheme%aintr,vals(i1:i2))
! val0           = valsint(i1)
! end do
! return
! endif


! Handle the general case

a0 = chebscheme%ab(1,intc)
b0 = chebscheme%ab(2,intc)

i1 = (intc-1)*k+1
i2 = intc*k

cc   = (2*c - (b0+a0) ) /(b0-a0)


!
!  Evaluate the integral on the interval containing c
!

vals0  = matmul(chebscheme%aintl,vals(i1:i2))
coefs0 = matmul(chebscheme%acoefs,vals(i1:i2))

call chebs(k+1,cc,polsc)

vals0 = vals0 + coefs0(1)*(-cc-1)
vals0 = vals0 + coefs0(2)*(0.5d0-cc**2/2)

dsign  = -1.0d0
do i=2,k-1
dd1    = (1.0d0 /(2.0d0*(i+1)) - 1.0d0/(2.0d0*(i-1)) ) * dsign
dd2    = 1.0d0 /(2.0d0*(i+1))*polsc(i+2) - 1.0d0/(2.0d0*(i-1))*polsc(i) 
vals0  = vals0 + coefs0(i+1)*(dd1-dd2)
dsign  = -dsign
end do

valsint(i1:i2) = valc + (b0-a0)/2*vals0
vall           = valsint(i1)
valr           = valsint(i2)


!
!  Integrate forward from the interval containing c
!

val0 = valr
do int=intc+1,nints
a0 = chebscheme%ab(1,int)
b0 = chebscheme%ab(2,int)
i1 = (int-1)*k + 1
i2 = int*k
valsint(i1:i2) = val0 + (b0-a0)/2 * matmul(chebscheme%aintl,vals(i1:i2))
val0           = valsint(i2)
end do


!
!  Now integrate backward from the interval containing c
!

val0 = vall
do int=intc-1,1,-1
a0 = chebscheme%ab(1,int)
b0 = chebscheme%ab(2,int)
i1 = (int-1)*k + 1
i2 = int*k
valsint(i1:i2) = val0 + (b0-a0)/2 * matmul(chebscheme%aintr,vals(i1:i2))
val0           = valsint(i1)
end do

end subroutine


subroutine chebpw_int2(chebscheme,vals,c,valc,valsint)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double complex                                    :: vals(:)
double complex, allocatable, intent(out)          :: valsint(:)
double complex                                    :: valc
!
!  Compute the values of the function
!
!                     t
!    g(t) = valc + int    f(s) ds
!                     c
!
!  with c a point in [a,b].
!  
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    valsint - the vector of values (4) of the antiderivative of the expansion
!     which is 0 at the left-hand endpoint of [a,b]
!

double precision     ::    polsc(chebscheme%k+1)
double complex       ::    vals0(chebscheme%k), coefs0(chebscheme%k), val0, vall, valr


nints = chebscheme%nints
k     = chebscheme%k

allocate(valsint(k*nints))

if (.not. allocated(chebscheme%aintl)) then
call chebyshev_intlmatrix(chebscheme%k,chebscheme%aintl)
endif

if (.not. allocated(chebscheme%aintr)) then
call chebyshev_intrmatrix(chebscheme%k,chebscheme%aintr)
endif

!
!  Find the interval containing the point c 
!

do intc=1,nints-1
if (c .lt. chebscheme%ab(2,intc)) exit
end do

! Handle the special cases where intc==1 or intc==nints

! if (intc .eq. 1) then

! val0 = valc
! do int=1,nints
! a0 = chebscheme%ab(1,int)
! b0 = chebscheme%ab(2,int)
! i1 = (int-1)*k + 1
! i2 = int*k
! valsint(i1:i2) = val0 + (b0-a0)/2 * matmul(chebscheme%aintl,vals(i1:i2))
! val0           = valsint(i2)
! end do
! return
! endif


! if (intc .eq. nints) then
! val0 = valc
! do int=nints,1,-1
! a0 = chebscheme%ab(1,int)
! b0 = chebscheme%ab(2,int)
! i1 = (int-1)*k + 1
! i2 = int*k
! valsint(i1:i2) = val0 + (b0-a0)/2 * matmul(chebscheme%aintr,vals(i1:i2))
! val0           = valsint(i1)
! end do
! return
! endif


! Handle the general case

a0 = chebscheme%ab(1,intc)
b0 = chebscheme%ab(2,intc)

i1 = (intc-1)*k+1
i2 = intc*k

cc   = (2*c - (b0+a0) ) /(b0-a0)


!
!  Evaluate the integral on the interval containing c
!

vals0  = matmul(chebscheme%aintl,vals(i1:i2))
coefs0 = matmul(chebscheme%acoefs,vals(i1:i2))

call chebs(k+1,cc,polsc)

vals0 = vals0 + coefs0(1)*(-cc-1)
vals0 = vals0 + coefs0(2)*(0.5d0-cc**2/2)

dsign  = -1.0d0
do i=2,k-1
dd1    = (1.0d0 /(2.0d0*(i+1)) - 1.0d0/(2.0d0*(i-1)) ) * dsign
dd2    = 1.0d0 /(2.0d0*(i+1))*polsc(i+2) - 1.0d0/(2.0d0*(i-1))*polsc(i) 
vals0  = vals0 + coefs0(i+1)*(dd1-dd2)
dsign  = -dsign
end do

valsint(i1:i2) = valc + (b0-a0)/2*vals0
vall           = valsint(i1)
valr           = valsint(i2)


!
!  Integrate forward from the interval containing c
!

val0 = valr
do int=intc+1,nints
a0 = chebscheme%ab(1,int)
b0 = chebscheme%ab(2,int)
i1 = (int-1)*k + 1
i2 = int*k
valsint(i1:i2) = val0 + (b0-a0)/2 * matmul(chebscheme%aintl,vals(i1:i2))
val0           = valsint(i2)
end do


!
!  Now integrate backward from the interval containing c
!

val0 = vall
do int=intc-1,1,-1
a0 = chebscheme%ab(1,int)
b0 = chebscheme%ab(2,int)
i1 = (int-1)*k + 1
i2 = int*k
valsint(i1:i2) = val0 + (b0-a0)/2 * matmul(chebscheme%aintr,vals(i1:i2))
val0           = valsint(i1)
end do

end subroutine



subroutine chebpw_diff1(chebscheme,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double precision                                  :: vals(:)
double precision, allocatable, intent(out)        :: ders(:)
!
!  Perform spectral differentiation on an expansion; that is, given the vector (4)
!  of values of the an expansion of the form (1), return the vector of values
!  of the derivative of the expansion.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    ders - the vector of values (4) of the derivative of the expansion
!

nints = chebscheme%nints
k     = chebscheme%k

allocate(ders(k*nints))

if (.not. allocated(chebscheme%adiff)) then
call chebyshev_diffmatrix(chebscheme%k,chebscheme%adiff)
endif

i1    = 1
i2    = k

do int=1,nints
a            = chebscheme%ab(1,int)
b            = chebscheme%ab(2,int)
ders(i1:i2)  = 2/(b-a)*matmul(chebscheme%adiff,vals(i1:i2))
i1           = i1+k
i2           = i2+k
end do

end subroutine

subroutine chebpw_diff2(chebscheme,nfuns,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double precision                                  :: vals(:,:)
double precision, allocatable, intent(out)        :: ders(:,:)

!
!  Perform spectral differentiation on an expansion; that is, given the vector (4)
!  of values of the an expansion of the form (1), return the vector of values
!  of the derivative of the expansion.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    ders - the vector of values (4) of the derivative of the expansion
!

nints = chebscheme%nints
k     = chebscheme%k

allocate(ders(k*nints,nfuns))

if (.not. allocated(chebscheme%adiff)) then
call chebyshev_diffmatrix(chebscheme%k,chebscheme%adiff)
endif

i1    = 1
i2    = k

do int=1,nints
a              = chebscheme%ab(1,int)
b              = chebscheme%ab(2,int)
ders(i1:i2,:)  = 2/(b-a)*matmul(chebscheme%adiff,vals(i1:i2,:))
i1             = i1+k
i2             = i2+k
end do

end subroutine



subroutine chebpw_diff3(chebscheme,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double complex                                    :: vals(:)
double complex, allocatable, intent(out)          :: ders(:)
!
!  Perform spectral differentiation on an expansion; that is, given the vector (4)
!  of values of the an expansion of the form (1), return the vector of values
!  of the derivative of the expansion.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    ders - the vector of values (4) of the derivative of the expansion
!

nints = chebscheme%nints
k     = chebscheme%k

allocate(ders(k*nints))

if (.not. allocated(chebscheme%adiff)) then
call chebyshev_diffmatrix(chebscheme%k,chebscheme%adiff)
endif

i1    = 1
i2    = k

do int=1,nints
a            = chebscheme%ab(1,int)
b            = chebscheme%ab(2,int)
ders(i1:i2)  = 2/(b-a)*matmul(chebscheme%adiff,vals(i1:i2))
i1           = i1+k
i2           = i2+k
end do

end subroutine

subroutine chebpw_diff4(chebscheme,nfuns,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double complex                                    :: vals(:,:)
double complex, allocatable, intent(out)          :: ders(:,:)
!
!  Perform spectral differentiation on an expansion; that is, given the vector (4)
!  of values of the an expansion of the form (1), return the vector of values
!  of the derivative of the expansion.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    ders - the vector of values (4) of the derivative of the expansion
!

nints = chebscheme%nints
k     = chebscheme%k

allocate(ders(k*nints,nfuns))

if (.not. allocated(chebscheme%adiff)) then
call chebyshev_diffmatrix(chebscheme%k,chebscheme%adiff)
endif

i1    = 1
i2    = k

do int=1,nints
a              = chebscheme%ab(1,int)
b              = chebscheme%ab(2,int)
ders(i1:i2,:)  = 2/(b-a)*matmul(chebscheme%adiff,vals(i1:i2,:))
i1             = i1+k
i2             = i2+k
end do

end subroutine



subroutine chebpw_intl1(chebscheme,vals,valsint)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double precision                                  :: vals(:)
double precision, allocatable, intent(out)        :: valsint(:)
!
!  Perform "left" spectral integration on a single real-valued expansion.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    valsint - the vector of values (4) of the antiderivative of the expansion
!     which is 0 at the left-hand endpoint of [a,b]
!

nints = chebscheme%nints
k     = chebscheme%k

allocate(valsint(k*nints))

if (.not. allocated(chebscheme%aintl)) then
call chebyshev_intlmatrix(chebscheme%k,chebscheme%aintl)
endif

i1    = 1
i2    = k
val0  = 0

do int=1,nints
a              = chebscheme%ab(1,int)
b              = chebscheme%ab(2,int)
valsint(i1:i2) = val0 + (b-a)/2*matmul(chebscheme%aintl,vals(i1:i2))
val0           = valsint(i2)
i1             = i1+k
i2             = i2+k
end do

end subroutine

subroutine chebpw_intl2(chebscheme,nfuns,vals,valsint)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double precision                                  :: vals(:,:)
double precision, allocatable, intent(out)        :: valsint(:,:)
!
!  Perform "left" spectral integration on a collection of real-valued expansions.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    valsint - the vector of values (4) of the antiderivative of the expansion
!     which is 0 at the left-hand endpoint of [a,b]
!

double precision :: vals0(nfuns)

nints = chebscheme%nints
k     = chebscheme%k

allocate(valsint(k*nints,nfuns))

if (.not. allocated(chebscheme%aintl)) then
call chebyshev_intlmatrix(chebscheme%k,chebscheme%aintl)
endif

i1     = 1
i2     = k
vals0  = 0

do int=1,nints
a                = chebscheme%ab(1,int)
b                = chebscheme%ab(2,int)
valsint(i1:i2,:) = (b-a)/2*matmul(chebscheme%aintl,vals(i1:i2,:))

do i=1,nfuns
valsint(:,i) = valsint(:,i)+vals0(i)
vals0(i)     = valsint(i2,i)
end do

i1               = i1+k
i2               = i2+k
end do

end subroutine


subroutine chebpw_intl3(chebscheme,vals,valsint)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double complex                                    :: vals(:)
double complex, allocatable, intent(out)          :: valsint(:)
!
!  Perform "left" spectral integration on a single complex-valued expansion.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    valsint - the vector of values (4) of the antiderivative of the expansion
!     which is 0 at the left-hand endpoint of [a,b]
!
double complex                                    :: val0

nints = chebscheme%nints
k     = chebscheme%k

allocate(valsint(k*nints))

if (.not. allocated(chebscheme%aintl)) then
call chebyshev_intlmatrix(chebscheme%k,chebscheme%aintl)
endif

i1    = 1
i2    = k
val0  = 0

do int=1,nints
a              = chebscheme%ab(1,int)
b              = chebscheme%ab(2,int)
valsint(i1:i2) = val0 + (b-a)/2*matmul(chebscheme%aintl,vals(i1:i2))
val0           = valsint(i2)
i1             = i1+k
i2             = i2+k
end do

end subroutine

subroutine chebpw_intl4(chebscheme,nfuns,vals,valsint)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                               :: chebscheme
double complex                                    :: vals(:,:)
double complex, allocatable, intent(out)          :: valsint(:,:)
!
!  Perform "left" spectral integration on a collection of complex-valued expansions.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4) of the expansion
!
!  Output parameters:
!    valsint - the vector of values (4) of the antiderivative of the expansion
!     which is 0 at the left-hand endpoint of [a,b]
!

double complex :: vals0(nfuns)

nints = chebscheme%nints
k     = chebscheme%k

allocate(valsint(k*nints,nfuns))

if (.not. allocated(chebscheme%aintl)) then
call chebyshev_intlmatrix(chebscheme%k,chebscheme%aintl)
endif

i1     = 1
i2     = k
vals0  = 0

do int=1,nints
a                = chebscheme%ab(1,int)
b                = chebscheme%ab(2,int)
valsint(i1:i2,:) = (b-a)/2*matmul(chebscheme%aintl,vals(i1:i2,:))

do i=1,nfuns
valsint(:,i) = valsint(:,i)+vals0(i)
vals0(i)     = valsint(i2,i)
end do

i1               = i1+k
i2               = i2+k
end do

end subroutine


subroutine chebpw_roots(chebscheme,coefs,nroots,roots)
implicit double precision (a-h,o-z)
type(chebpw_scheme)                        :: chebscheme
double precision                           :: coefs(:)
double precision, intent(out), allocatable :: roots(:)
!
!  Use a colleague matrix approach to find all of the roots of a real-valued
!  expansion of the form (1) on the interval [a,b] given its vector of expansion 
!  coefficients (3).
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - the vector of coefficients (3) of an expansion of the form (1)
!
!  Output parameters:
!    nroots - the number of roots of the expansion on [a,b]
!    roots - an array containing the roots
!

double precision, allocatable :: roots0(:), roots00(:)

nints  = chebscheme%nints
k      = chebscheme%k
nroots = 0
x0     = 1d300
eps0   = sqrt(epsilon(0.0))

allocate(roots0(k))
allocate(roots00(nints*k))

do int=1,nints
i1 = (int-1)*k+1
i2 = int*k
a  = chebscheme%ab(1,int)
b  = chebscheme%ab(2,int)

call chebyshev_roots(k,coefs(i1:i2),nroots0,roots0)

roots0 = roots0*(b-a)/2 + (b+a)/2

do i=1,nroots0
x               = roots0(i)

if ( abs(x-x0) .lt. eps0) cycle

nroots          = nroots+1
roots00(nroots) = x
x0              = x
end do

! roots00(nroots+1:nroots+nroots0) = roots0(1:nroots0)
! nroots                           = nroots+nroots0

end do


allocate(roots(nroots))
roots = roots00(1:nroots)

end subroutine




subroutine chebpw_plot1(chebscheme,vals,a,b,scriptname,filename)
implicit double precision (a-h,o-z)
type(chebpw_scheme)            :: chebscheme
double precision               :: vals(:)
character(len=*)               :: scriptname,filename
!
!  Produce a simple plot of a real-valued function represented via a piecewise 
!  Chebyshev    discretization scheme.  More explicitly, generate a python script which
!  produces a PDF file with the desired plot, and exceute the python script
!  in the background.
!
!  Input parameters:
!    chebscheme - structure describing the piecewise Chebyshev scheme
!    vals - the vector of values of the expansion
!    (a,b) - the interval on which to plot the function
!    scriptname - the name for the python script
!    filename - the name for the PDF file
!
!  Output parameters:
!    N/A
!

double precision, allocatable :: xs(:), ys(:)

nints = chebscheme%nints
k     = chebscheme%k

n     = 1000
allocate(xs(n),ys(n))

do i=1,n
x     = a + (b-a)*(i-1.0d0)/(n-1.0d0)
xs(i) = x
call chebpw_interp(chebscheme,vals,x,ys(i))
end do


x1 = 0.0d0
x2 = -1.0d0
y1 = 0.0d0
y2  =-1.0d0

iflogx = 0
iflogy = 0
nfuns  = 1

call plot_functions(scriptname,filename,nfuns,"","",iflogx,iflogy,x1,x2,y1,y2, &
   n,xs,ys,"","","")

end subroutine


subroutine chebpw_plot2(chebscheme,vals,a,b,scriptname,filename)
implicit double precision (a-h,o-z)
type(chebpw_scheme)            :: chebscheme
double complex                 :: vals(:)
character(len=*)               :: scriptname,filename
!
!  Produce a simple plot of a real-valued function represented via a piecewise 
!  Chebyshev    discretization scheme.  More explicitly, generate a python script which
!  produces a PDF file with the desired plot, and exceute the python script
!  in the background.
!
!  Input parameters:
!    chebscheme - structure describing the piecewise Chebyshev scheme
!    vals - the vector of values of the expansion
!    (a,b) - the interval on which to plot the function
!    scriptname - the name for the python script
!    filename - the name for the PDF file
!
!  Output parameters:
!    N/A
!

double precision, allocatable :: xs(:), ys(:,:)
double complex                :: val

nints = chebscheme%nints
k     = chebscheme%k

n     = 1000
allocate(xs(n),ys(2,n))

do i=1,n
x     = a + (b-a)*(i+0.0d0)/(n+1.0d0)
xs(i) = x
call chebpw_interp(chebscheme,vals,x,val)

ys(1,i) = real(val)
ys(2,i) = imag(val)

end do



x1 = 0.0d0
x2 = -1.0d0
y1 = 0.0d0
y2  =-1.0d0
iflogx = 0
iflogy = 0
nfuns  = 2

call plot_functions(scriptname,filename,nfuns,"","",iflogx,iflogy,x1,x2,y1,y2, &
   n,xs,ys,"","","")


end subroutine


subroutine chebpw_plot3(chebscheme,vals,scriptname,filename)
implicit double precision (a-h,o-z)
type(chebpw_scheme)            :: chebscheme
double precision               :: vals(:)
character(len=*)               :: scriptname,filename
!
!  Produce a simple plot of a real-valued function represented via a piecewise 
!  Chebyshev    discretization scheme.  More explicitly, generate a python script which
!  produces a PDF file with the desired plot, and exceute the python script
!  in the background.
!
!  Input parameters:
!    chebscheme - structure describing the piecewise Chebyshev scheme
!    vals - the vector of values of the expansion
!    scriptname - the name for the python script
!    filename - the name for the PDF file
!
!  Output parameters:
!    N/A
!

double precision, allocatable :: xs(:), ys(:)

nints = chebscheme%nints
k     = chebscheme%k
a     = chebscheme%a
b     = chebscheme%b


n     = 1000
allocate(xs(n),ys(n))

do i=1,n
x     = a + (b-a)*(i-1.0d0)/(n-1.0d0)
xs(i) = x
call chebpw_interp(chebscheme,vals,x,ys(i))
end do


x1 = 0.0d0
x2 = -1.0d0
y1 = 0.0d0
y2  =-1.0d0
iflogx = 0
iflogy = 0
nfuns  = 1

call plot_functions(scriptname,filename,nfuns,"","",iflogx,iflogy,x1,x2,y1,y2, &
   n,xs,ys,"","","")


end subroutine


subroutine chebpw_plot4(chebscheme,vals,scriptname,filename)
implicit double precision (a-h,o-z)
type(chebpw_scheme)            :: chebscheme
double complex                 :: vals(:)
character(len=*)               :: scriptname,filename
!
!  Produce a simple plot of a real-valued function represented via a piecewise 
!  Chebyshev    discretization scheme.  More explicitly, generate a python script which
!  produces a PDF file with the desired plot, and exceute the python script
!  in the background.
!
!  Input parameters:
!    chebscheme - structure describing the piecewise Chebyshev scheme
!    vals - the vector of values of the expansion
!    scriptname - the name for the python script
!    filename - the name for the PDF file
!
!  Output parameters:
!    N/A
!

double precision, allocatable :: xs(:), ys(:,:)
double complex                :: val

nints = chebscheme%nints
k     = chebscheme%k
a     = chebscheme%a
b     = chebscheme%b

n     = 1000
allocate(xs(n),ys(2,n))

do i=1,n
x     = a + (b-a)*(i+0.0d0)/(n+1.0d0)
xs(i) = x
call chebpw_interp(chebscheme,vals,x,val)

ys(1,i) = real(val)
ys(2,i) = imag(val)

end do



x1 = 0.0d0
x2 = -1.0d0
y1 = 0.0d0
y2  =-1.0d0
iflogx = 0
iflogy = 0
nfuns  = 2

call plot_functions(scriptname,filename,nfuns,"","",iflogx,iflogy,x1,x2,y1,y2, &
   n,xs,ys,"","","")


end subroutine


subroutine chebpw_plot5(chebscheme,vals1,vals2,a,b,scriptname,filename)
implicit double precision (a-h,o-z)
type(chebpw_scheme)            :: chebscheme
double precision               :: vals1(:), vals2(:)
character(len=*)               :: scriptname,filename
!
!  Produce a simple plot of two real-valued functions represented via a piecewise 
!  Chebyshev discretization scheme.  More explicitly, generate a python script which
!  produces a PDF file with the desired plot, and exceute the python script
!  in the background.
!
!  Input parameters:
!    chebscheme - structure describing the piecewise Chebyshev scheme
!    vals - the vector of values of the expansion
!    (a,b) - the interval on which to plot the function
!    scriptname - the name for the python script
!    filename - the name for the PDF file
!
!  Output parameters:
!    N/A
!

double precision, allocatable :: xs(:), ys(:,:)

nints = chebscheme%nints
k     = chebscheme%k

n     = 5000
allocate(xs(n),ys(2,n))

do i=1,n
x     = a + (b-a)*(i-1.0d0)/(n-1.0d0)
xs(i) = x
call chebpw_interp(chebscheme,vals1,vals2,x,ys(1,i),ys(2,i))
end do


x1 = 0.0d0
x2 = -1.0d0
y1 = 0.0d0
y2  =-1.0d0
iflogx = 0
iflogy = 0
nfuns  = 2

call plot_functions(scriptname,filename,nfuns,"","",iflogx,iflogy,x1,x2,y1,y2, &
   n,xs,ys,"","","")

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Routines for evaluating one or more functions at a single point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine chebpw_eval1(chebscheme,coefs,x,val)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs(:), val
!
!  Evaluate a real-valued expansion of the form (1) given the vector (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - the vector of coefficients (3)
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the function at the point x
! 


nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k

idx = i2
xxx = 2*xx
a0  = coefs(idx)
a1  = coefs(idx-1)
idx = idx-2

do i=2,k-1
y   = a1
a1  = coefs(idx)-a0
a0  = y + a0*xxx
idx = idx-1
end do

val = a1 + a0 * xx;

end subroutine


subroutine chebpw_eval2(chebscheme,nfuns,coefs,x,vals)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs(:,:), vals(:)
!
!  Evaluate several real-valued expansions of the form (1) given their vectors (3) of
!  coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - an array whose jth vector is the vector of coefficients (3)
!     for the jth expansion
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals - the values of the expansions at the point x
! 

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k

xxx = 2*xx

do j=1,nfuns
idx = i2
a0  = coefs(idx,j)
a1  = coefs(idx-1,j)
idx = idx-2

do i=2,k-1
y   = a1
a1  = coefs(idx,j)-a0
a0  = y + a0*xxx
idx = idx-1
end do

vals(j) = a1 + a0 * xx;
end do

end subroutine

subroutine chebpw_eval3(chebscheme,coefs1,coefs2,x,val1,val2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs1(:), coefs2(:)
!
!  Evaluate two real-valued expansions of the form (1) given their vectors (3) of 
!  coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs1 - the coefficients vector for the first expansion
!    coefs2 - the coefficients vector for the second expansion
!    x - the point in [a,b] at which to evaluate the functions
!
!  Output parameters:
!    val1 - the value of the first expansion
!    val2 - the value of the second expansion


nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k

xxx = 2*xx

idx = i2
a0  = coefs1(idx)
a1  = coefs1(idx-1)

b0  = coefs2(idx)
b1  = coefs2(idx-1)

idx = idx-2

do i=2,k-1
y   = a1
z   = b1

a1  = coefs1(idx)-a0
a0  = y + a0*xxx

b1  = coefs2(idx)-b0
b0  = z + b0*xxx

idx = idx-1
end do

val1 = a1 + a0 * xx;
val2 = b1 + b0 * xx;

end subroutine


subroutine chebpw_eval4(chebscheme,coefs,x,val)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: coefs(:), val
!
!  Evaluate a complex-valued expansion of the form (1) given the vector (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - the vector of coefficients (3)
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the function at the point x
! 

double complex :: a0, a1, y

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k

idx = i2
xxx = 2*xx
a0  = coefs(idx)
a1  = coefs(idx-1)
idx = idx-2

do i=2,k-1
y   = a1
a1  = coefs(idx)-a0
a0  = y + a0*xxx
idx = idx-1
end do

val = a1 + a0 * xx;

end subroutine


subroutine chebpw_eval5(chebscheme,nfuns,coefs,x,vals)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: coefs(:,:), vals(:)
!
!  Evaluate several complex-valued expansions of the form (1) given their vectors 
!  (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - an array whose jth vector is the vector of coefficients (3)
!     for the jth expansion
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals - the values of the expansions at the point x
! 

double complex :: a0, a1, y

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k

xxx = 2*xx

do j=1,nfuns
idx = i2
a0  = coefs(idx,j)
a1  = coefs(idx-1,j)
idx = idx-2

do i=2,k-1
y   = a1
a1  = coefs(idx,j)-a0
a0  = y + a0*xxx
idx = idx-1
end do

vals(j) = a1 + a0 * xx;
end do

end subroutine

subroutine chebpw_eval6(chebscheme,coefs1,coefs2,x,val1,val2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: coefs1(:), coefs2(:), val1, val2
!
!  Evaluate two complex-valued expansions of the form (1) given their vectors (3) of 
!  coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs1 - the coefficients vector for the first expansion
!    coefs2 - the coefficients vector for the second expansion
!    x - the point in [a,b] at which to evaluate the functions
!
!  Output parameters:
!    val1 - the value of the first expansion
!    val2 - the value of the second expansion

double complex :: a0, a1, y
double complex :: b0, b1, z

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k

xxx = 2*xx

idx = i2
a0  = coefs1(idx)
a1  = coefs1(idx-1)

b0  = coefs2(idx)
b1  = coefs2(idx-1)

idx = idx-2

do i=2,k-1
y   = a1
z   = b1

a1  = coefs1(idx)-a0
a0  = y + a0*xxx

b1  = coefs2(idx)-b0
b0  = z + b0*xxx

idx = idx-1
end do

val1 = a1 + a0 * xx;
val2 = b1 + b0 * xx;

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Routines for evaluating the functions and their derivatives at a single point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine chebpw_evalder1(chebscheme,coefs,x,val,der)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs(:), val
!
!  Evaluate a real-valued expansion of the form (1) and its derivative at a single point
!  given the vector (3) of coefficients.
!
!  Input parameters:
!    chebcheme - structure describing a piecewise discretization scheme
!    coefs - the vector of coefficients (3)
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the expansion at the point x
!    der - the derivative of the expansion at the point 
! 
double precision :: pols(chebscheme%k), ders(chebscheme%k)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k


!  Chebnomial evaluation 
call chebders(k,xx,pols,ders)
val = dot_product(pols,coefs(i1:i2))
der = 2/(b-a)*dot_product(ders,coefs(i1:i2))

if (itype .eq. 3) then
dd  = sqrt(2/(b-a))
val = val * dd
der = der * dd
endif

end subroutine


subroutine chebpw_evalder2(chebscheme,nfuns,coefs,x,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs(:,:), vals(:), ders(:)
!
!  Evaluate several real-valued expansions of the form (1) given their vectors (3) of
!  coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a piecewise discretization scheme
!    coefs - an array whose jth vector is the vector of coefficients (3)
!     for the jth expansion
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals - the values of the expansions at the point x
! 

double precision :: pols(chebscheme%k), polders(chebscheme%k)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k


! Evaluation via chebnomials
call chebders(k,xx,pols,polders)
vals = matmul(pols,coefs(i1:i2,:))
ders = 2/(b-a)*matmul(polders,coefs(i1:i2,:))

if (itype .eq. 3) then
dd   = sqrt(2/(b-a))
vals = vals * dd
ders = ders * dd
endif

end subroutine



subroutine chebpw_evalder3(chebscheme,coefs1,coefs2,x,val1,val2,der1,der2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs1(:), coefs2(:), val1, val2, der1, der2
!
!  Evaluate two real-valued expansions of the form (1) given their vectors (3) of 
!  coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a piecewise discretization scheme
!    coefs1 - the coefficients vector for the first expansion
!    coefs2 - the coefficients vector for the second expansion
!    x - the point in [a,b] at which to evaluate the functions
!
!  Output parameters:
!    val1 - the value of the first expansion
!    val2 - the value of the second expansion


double precision :: pols(chebscheme%k), polders(chebscheme%k)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k

! Chebnomial Evaluation

call chebders(k,xx,pols,polders)
val1 = dot_product(pols,coefs1(i1:i2))
val2 = dot_product(pols,coefs2(i1:i2))
der1 = 2/(b-a)*dot_product(polders,coefs1(i1:i2))
der2 = 2/(b-a)*dot_product(polders,coefs2(i1:i2))

if (itype .eq. 3) then
dd   = sqrt(2/(b-a))
val1 = val1 * dd
val2 = val2 * dd
der1 = der1 * dd
der2 = der2 * dd
endif

end subroutine


subroutine chebpw_evalder4(chebscheme,coefs,x,val,der)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: coefs(:), val, der
!
!  Evaluate a real-valued expansion of the form (1) and its derivative at a single point
!  given the vector (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a piecewise discretization scheme
!    coefs - the vector of coefficients (3)
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the expansion at the point x
!    der - the derivative of the expansion at the point 
! 
double precision :: pols(chebscheme%k), ders(chebscheme%k)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k


!  Chebnomial evaluation 
call chebders(k,xx,pols,ders)
val = dot_product(pols,coefs(i1:i2))
der = 2/(b-a)*dot_product(ders,coefs(i1:i2))

if (itype .eq. 3) then
dd   = sqrt(2/(b-a))
val = val * dd
der = der * dd
endif

end subroutine


subroutine chebpw_evalder5(chebscheme,nfuns,coefs,x,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: coefs(:,:), vals(:), ders(:)
!
!  Evaluate several real-valued expansions of the form (1) given their vectors (3) of
!  coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a piecewise discretization scheme
!    coefs - an array whose jth vector is the vector of coefficients (3)
!     for the jth expansion
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals - the values of the expansions at the point x
! 

double precision :: pols(chebscheme%k), polders(chebscheme%k)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k



! Evaluation chebnomials
call chebders(k,xx,pols,polders)
vals = matmul(pols,coefs(i1:i2,:))
ders = 2/(b-a)*matmul(polders,coefs(i1:i2,:))

if (itype .eq. 3) then
dd   = sqrt(2/(b-a))
vals = vals * dd
ders = ders * dd
endif

end subroutine



subroutine chebpw_evalder6(chebscheme,coefs1,coefs2,x,val1,val2,der1,der2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: coefs1(:), coefs2(:), val1, val2, der1, der2
!
!  Evaluate two real-valued expansions of the form (1) given their vectors (3) of 
!  coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a piecewise discretization scheme
!    coefs1 - the coefficients vector for the first expansion
!    coefs2 - the coefficients vector for the second expansion
!    x - the point in [a,b] at which to evaluate the functions
!
!  Output parameters:
!    val1 - the value of the first expansion
!    val2 - the value of the second expansion


double precision :: pols(chebscheme%k), polders(chebscheme%k)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a = chebscheme%ab(1,int)
b = chebscheme%ab(2,int)

xx   = (2*x - (b+a) ) /(b-a)

i1   = 1 + (int-1)*k
i2   = int*k

! Chebnomial Evaluation
call chebders(k,xx,pols,polders)
val1 = dot_product(pols,coefs1(i1:i2))
val2 = dot_product(pols,coefs2(i1:i2))
der1 = 2/(b-a)*dot_product(polders,coefs1(i1:i2))
der2 = 2/(b-a)*dot_product(polders,coefs2(i1:i2))

if (itype .eq. 3) then
dd   = sqrt(2/(b-a))
val1 = val1 * dd
val2 = val2 * dd
der1 = der1 * dd
der2 = der2 * dd
endif

end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Routines for evaluating one or more expansions at a single point using barycentric
!  interpolation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chebpw_interp1(chebscheme,vals,x,val)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: val, vals(:)
!
!  Evaluate a real-valued expansion represented via the vector of values (4) at a 
!  specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4)
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the expansion at the point x
! 

double precision       :: whts(chebscheme%k)
data pi / 3.14159265358979323846264338327950288d0 /

eps0 = epsilon(0.0d0)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a  = chebscheme%ab(1,int)
b  = chebscheme%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)
i1 = 1 + (int-1)*k
i2 = int*k

!
!  Perform barcyentric interpolation
!

dsign = 1.0d0
do i=1,k
!xxx     =  cos(pi*(k-i)/(k-1.0d0))
dd      = xx - chebscheme%xscheb(i)
if ( abs(dd) .le. eps0) then
val = vals(i1+i-1)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2

dsum2 = sum(whts)
val    = dot_product(whts,vals(i1:i2)) / dsum2

return

end subroutine


subroutine chebpw_interp2(chebscheme,nfuns,vals,x,valsout)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: valsout(:), vals(:,:)
!
!  Evaluate a collection of real-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals - an array whose jth column is the vector of values (4)
!      for the jth expansion
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals - the values of the expansions at the point x
! 

double precision       :: whts(chebscheme%k)
data pi / 3.14159265358979323846264338327950288d0 /


eps0 = epsilon(0.0d0)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a  = chebscheme%ab(1,int)
b  = chebscheme%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)
i1 = 1 + (int-1)*k
i2 = int*k

!
!  Perform barcyentric interpolation
!

dsign = 1.0d0
do i=1,k
!xxx     =  cos(pi*(k-i)/(k-1.0d0))
dd      = xx - chebscheme%xscheb(i)
if ( abs(dd) .le. eps0) then
valsout = vals(i1+i-1,:)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2

dsum2    = sum(whts)
valsout  = matmul(whts,vals(i1:i2,:)) / dsum2

return

end subroutine


subroutine chebpw_interp3(chebscheme,vals1,vals2,x,val1,val2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: vals1(:), vals2(:), val1, val2
!
!  Evaluate two real-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals1 - the vector (4) of values of the first expansion
!    vals2 - the vector (4) of values of the second expansion
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val1 - the value of the first expansion at x
!    val2 - the value of the second expansion at x

! 

double precision       :: whts(chebscheme%k)
data pi / 3.14159265358979323846264338327950288d0 /


eps0 = epsilon(0.0d0)
nints = chebscheme%nints
k     = chebscheme%k


!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a  = chebscheme%ab(1,int)
b  = chebscheme%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)
i1 = 1 + (int-1)*k
i2 = int*k


!
!  Perform barcyentric interpolation
!

dsign = 1.0d0
do i=1,k
!xxx     =  cos(pi*(k-i)/(k-1.0d0))
dd      = xx - chebscheme%xscheb(i)
if ( abs(dd) .le. eps0) then
! valsout = vals(i1+i-1,:)
val1 = vals1(i1+i-1)
val2 = vals2(i1+i-1)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2

dsum2    = sum(whts)
val1     = dot_product(whts,vals1(i1:i2)) / dsum2
val2     = dot_product(whts,vals2(i1:i2)) / dsum2

return
end subroutine


subroutine chebpw_interp4(chebscheme,vals,x,val)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: val, vals(:)
!
!  Evaluate a complex-valued expansion represented via the vector of values (4) at a 
!  specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4)
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the expansion at the point x
! 

double precision       :: whts(chebscheme%k)
data pi / 3.14159265358979323846264338327950288d0 /

eps0 = epsilon(0.0d0)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a  = chebscheme%ab(1,int)
b  = chebscheme%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)
i1 = 1 + (int-1)*k
i2 = int*k

!
!  Perform barcyentric interpolation
!

dsign = 1.0d0
do i=1,k
!xxx     =  cos(pi*(k-i)/(k-1.0d0))
dd      = xx - chebscheme%xscheb(i)
if ( abs(dd) .le. eps0) then
val = vals(i1+i-1)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2

dsum2 = sum(whts)
val    = dot_product(whts,vals(i1:i2)) / dsum2

return

end subroutine


subroutine chebpw_interp5(chebscheme,nfuns,vals,x,valsout)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex    :: valsout(:), vals(:,:)
!
!  Evaluate a collection of complex-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals - an array whose jth column is the vector of values (4)
!      for the jth expansion
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals - the values of the expansions at the point x
! 

double precision       :: whts(chebscheme%k)
data pi / 3.14159265358979323846264338327950288d0 /


eps0 = epsilon(0.0d0)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a  = chebscheme%ab(1,int)
b  = chebscheme%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)
i1 = 1 + (int-1)*k
i2 = int*k

!
!  Perform barcyentric interpolation
!

dsign = 1.0d0
do i=1,k
!xxx     =  cos(pi*(k-i)/(k-1.0d0))
dd      = xx - chebscheme%xscheb(i)
if ( abs(dd) .le. eps0) then
valsout = vals(i1+i-1,:)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2

dsum2    = sum(whts)
valsout  = matmul(whts,vals(i1:i2,:)) / dsum2

return

end subroutine


subroutine chebpw_interp6(chebscheme,vals1,vals2,x,val1,val2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: vals1(:), vals2(:), val1, val2
!
!  Evaluate two complex-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals1 - the vector (4) of values of the first expansion
!    vals2 - the vector (4) of values of the second expansion
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val1 - the value of the first expansion at x
!    val2 - the value of the second expansion at x

! 

double precision       :: whts(chebscheme%k)
data pi / 3.14159265358979323846264338327950288d0 /


eps0 = epsilon(0.0d0)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a  = chebscheme%ab(1,int)
b  = chebscheme%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)
i1 = 1 + (int-1)*k
i2 = int*k

!
!  Perform barcyentric interpolation
!

dsign = 1.0d0
do i=1,k
dd      = xx - chebscheme%xscheb(i)
if ( abs(dd) .le. eps0) then
val1 = vals1(i1+i-1)
val2 = vals2(i1+i-1)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2

dsum2    = sum(whts)
val1     = dot_product(whts,vals1(i1:i2)) / dsum2
val2     = dot_product(whts,vals2(i1:i2)) / dsum2

return
end subroutine


subroutine chebpw_interp7(chebscheme,vals1,vals2,vals3,vals4,x,val1,val2,val3,val4)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: vals1(:), vals2(:), vals3(:), vals4(:),val1, val2, val3, val4
!
!  Evaluate four real-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals? - the vector (4) of values of the expansions
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val? - the values of the expansions at x

! 

double precision       :: whts(chebscheme%k)
data pi / 3.14159265358979323846264338327950288d0 /


eps0 = epsilon(0.0d0)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a  = chebscheme%ab(1,int)
b  = chebscheme%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)
i1 = 1 + (int-1)*k
i2 = int*k

!
!  Perform barcyentric interpolation
!

dsign = 1.0d0
do i=1,k
dd      = xx - chebscheme%xscheb(i)
if ( abs(dd) .le. eps0) then
val1 = vals1(i1+i-1)
val2 = vals2(i1+i-1)
val3 = vals3(i1+i-1)
val4 = vals4(i1+i-1)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2

dsum2    = sum(whts)
val1     = dot_product(whts,vals1(i1:i2)) / dsum2
val2     = dot_product(whts,vals2(i1:i2)) / dsum2
val3     = dot_product(whts,vals3(i1:i2)) / dsum2
val4     = dot_product(whts,vals4(i1:i2)) / dsum2

return
end subroutine



subroutine chebpw_interp8(chebscheme,vals1,vals2,vals3,x,val1,val2,val3)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: vals1(:), vals2(:), vals3(:), val1, val2, val3
!
!  Evaluate three real-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals? - the vector (4) of values of the expansions
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val? - the values of the expansions at x

! 

double precision       :: whts(chebscheme%k)
data pi / 3.14159265358979323846264338327950288d0 /


eps0 = epsilon(0.0d0)

nints = chebscheme%nints
k     = chebscheme%k

!
!  Find the interval via a brute force search
!

do int = 1,nints-1
b = chebscheme%ab(2,int)
if (x .le. b) exit
end do

a  = chebscheme%ab(1,int)
b  = chebscheme%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)
i1 = 1 + (int-1)*k
i2 = int*k

!
!  Perform barcyentric interpolation
!

dsign = 1.0d0
do i=1,k
dd      = xx - chebscheme%xscheb(i)
if ( abs(dd) .le. eps0) then
val1 = vals1(i1+i-1)
val2 = vals2(i1+i-1)
val3 = vals3(i1+i-1)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2

dsum2    = sum(whts)
val1     = dot_product(whts,vals1(i1:i2)) / dsum2
val2     = dot_product(whts,vals2(i1:i2)) / dsum2
val3     = dot_product(whts,vals3(i1:i2)) / dsum2

return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Routines for interpolating evaluating one or more expansions at multiple points;
!  these are just convenience wrappers around the standard evaluation routines which
!  are not accelerated in any way.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine chebpw_interp101(chebscheme,vals,n,xs,valsout)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: valsout(:), vals(:)
double precision    :: xs(:)
!
!  Evaluate a real-valued expansion represented via the vector of values (4) at a 
!  collection of specified points.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4)
!    n - the number of evaluation points
!    x - the list of evaluation points
!
!  Output parameters:
!    vals - the values of the expansion 
! 

do i=1,n
call chebpw_interp1(chebscheme,vals,xs(i),valsout(i))
end do


end subroutine


subroutine chebpw_interp102(chebscheme,nfuns,vals,n,xs,valsout)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: valsout(:,:), vals(:,:)
double precision    :: xs(:)
!
!  Evaluate a collection of real-valued expansions represented via their vectors 
!  of values (4) at a collection of points.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals - an array whose jth column is the vector of values (4)
!      for the jth expansion
!    n - the number of evaluation points
!    xs - the list of points in [a,b] at which to evaluate the function
!
!  Output parameters:
!    valsout - the ith row of this array will contain the values of the
!     expansions at the ith input point
! 

do i=1,n
call chebpw_interp2(chebscheme,nfuns,vals,xs(i),valsout(i,:))
end do

end subroutine


subroutine chebpw_interp103(chebscheme,vals1,vals2,n,xs,valsout1,valsout2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: vals1(:), vals2(:), valsout1(:), valsout2(:)
double precision    :: xs(:)
!
!  Evaluate two real-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals1 - the vector (4) of values of the first expansion
!    vals2 - the vector (4) of values of the second expansion
!    n - the number of evaluation points
!    xs - the list of points in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals1 - the values of the first expansion at x
!    vals2 - the values of the second expansion at x
! 

do i=1,n
call chebpw_interp3(chebscheme,vals1,vals2,xs(i),valsout1(i),valsout2(i))
end do

end subroutine



subroutine chebpw_interp104(chebscheme,vals,n,xs,valsout)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: valsout(:), vals(:)
double precision    :: xs(:)
!
!  Evaluate a complex-valued expansion represented via the vector of values (4) at a 
!  collection of specified points.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    vals - the vector of values (4)
!    n - the number of evaluation points
!    x - the list of evaluation points
!
!  Output parameters:
!    vals - the values of the expansion 
! 

do i=1,n
call chebpw_interp4(chebscheme,vals,xs(i),valsout(i))
end do

end subroutine


subroutine chebpw_interp105(chebscheme,nfuns,vals,n,xs,valsout)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: valsout(:,:), vals(:,:)
double precision    :: xs(:)
!
!  Evaluate a collection of complex-valued expansions represented via their vectors 
!  of values (4) at a collection of points.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals - an array whose jth column is the vector of values (4)
!      for the jth expansion
!    n - the number of evaluation points
!    xs - the list of points in [a,b] at which to evaluate the function
!
!  Output parameters:
!    valsout - the ith row of this array will contain the values of the
!     expansions at the ith input point
! 

do i=1,n
call chebpw_interp5(chebscheme,nfuns,vals,xs(i),valsout(i,:))
end do

end subroutine


subroutine chebpw_interp106(chebscheme,vals1,vals2,n,xs,valsout1,valsout2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: vals1(:), vals2(:), valsout1(:), valsout2(:)
double precision    :: xs(:)
!
!  Evaluate two complex-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals1 - the vector (4) of values of the first expansion
!    vals2 - the vector (4) of values of the second expansion
!    n - the number of evaluation points
!    xs - the list of points in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals1 - the values of the first expansion at x
!    vals2 - the values of the second expansion at x
! 

do i=1,n
call chebpw_interp6(chebscheme,vals1,vals2,xs(i),valsout1(i),valsout2(i))
end do

end subroutine

subroutine chebpw_interp107(chebscheme,vals1,vals2,vals3,vals4,n,xs,valsout1, &
  valsout2,valsout3,valsout4)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: vals1(:), vals2(:), valsout1(:), valsout2(:)
double precision    :: vals3(:), vals4(:), valsout3(:), valsout4(:)
double precision    :: xs(:)
!
!  Evaluate two real-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals1 - the vector (4) of values of the first expansion
!    vals2 - the vector (4) of values of the second expansion
!    n - the number of evaluation points
!    xs - the list of points in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals1 - the values of the first expansion at x
!    vals2 - the values of the second expansion at x
! 

do i=1,n
call chebpw_interp7(chebscheme,vals1,vals2,vals3,vals4,xs(i),&
  valsout1(i),valsout2(i),valsout3(i),valsout4(i))
end do

end subroutine


subroutine chebpw_interp108(chebscheme,vals1,vals2,vals3,n,xs,valsout1, &
  valsout2,valsout3)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: vals1(:), vals2(:), valsout1(:), valsout2(:)
double precision    :: vals3(:), valsout3(:)
double precision    :: xs(:)
!
!  Evaluate three real-valued expansions represented via their vectors 
!  of values (4) at a specified point.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    nfuns - the number of expansions to evaluate
!    vals? - the vector (4) of values of the expansions
!    n - the number of evaluation points
!    xs - the list of points in [a,b] at which to evaluate the function
!
!  Output parameters:
!    valsout? - the values of the expansions
! 

do i=1,n
call chebpw_interp8(chebscheme,vals1,vals2,vals3,xs(i),valsout1(i),valsout2(i),&
  valsout3(i))
end do

end subroutine




subroutine chebpw_eval101(chebscheme,coefs,n,xs,vals)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: xs(:)
double precision    :: coefs(:), vals(:)
!
!  Evaluate a real-valued expansion of the form (1) at a collection of points
!  given the vector (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - the vector of coefficients (3)
!    n - the number of points at which to evaluate the expansion
!    xs - the list of evaluation points
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals - the values of the function at the specified points
! 

do i=1,n
call  chebpw_eval1(chebscheme,coefs,xs(i),vals(i))
end do

end subroutine


subroutine chebpw_eval102(chebscheme,nfuns,coefs,n,xs,vals)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs(:,:), vals(:,:)
double precision    :: xs(:)
!
!  Evaluate several real-valued expansions of the form (1) at a collection of
!  points given their vectors (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - an array whose jth vector is the vector of coefficients (3)
!     for the jth expansion
!    n - the number of points at which to evaluate the expansions
!    xs - the list of evaluation points
!
!  Output parameters:
!    vals - the ith row of this (n,nfuns) array gives the values of the
!      expansions at the ith input point
! 

do i=1,n
call chebpw_eval2(chebscheme,nfuns,coefs,xs(i),vals(i,:))
end do


end subroutine

subroutine chebpw_eval103(chebscheme,coefs1,coefs2,n,xs,vals1,vals2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs1(:), coefs2(:),vals1(:),vals2(:)
double precision    :: xs(:)
!
!  Evaluate two real-valued expansions of the form (1) at a collection of points
!  given their vectors (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs1 - the coefficients vector for the first expansion
!    coefs2 - the coefficients vector for the second expansion
!    n - the number of input points
!    xs - the points in [a,b] at which to evaluate the functions
!
!  Output parameters:
!    vals1 - the values of the first expansion
!    vals2 - the values of the second expansion

do i=1,n
call chebpw_eval3(chebscheme,coefs1,coefs2,xs(i),vals1(i),vals2(i))
end do

end subroutine



subroutine chebpw_eval104(chebscheme,coefs,n,xs,vals)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: xs(:)
double complex      :: coefs(:), vals(:)
!
!  Evaluate a complex-valued expansion of the form (1) at a collection of points
!  given the vector (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - the vector of coefficients (3)
!    n - the number of points at which to evaluate the expansion
!    xs - the list of evaluation points
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    vals - the values of the function at the specified points
! 

do i=1,n
call  chebpw_eval4(chebscheme,coefs,xs(i),vals(i))
end do

end subroutine


subroutine chebpw_eval105(chebscheme,nfuns,coefs,n,xs,vals)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: coefs(:,:), vals(:,:)
double precision    :: xs(:)
!
!  Evaluate several complex-valued expansions of the form (1) at a collection of
!  points given their vectors (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - an array whose jth vector is the vector of coefficients (3)
!     for the jth expansion
!    n - the number of points at which to evaluate the expansions
!    xs - the list of evaluation points
!
!  Output parameters:
!    vals - the ith row of this (n,nfuns) array gives the values of the
!      expansions at the ith input point
! 

do i=1,n
call chebpw_eval5(chebscheme,nfuns,coefs,xs(i),vals(i,:))
end do


end subroutine

subroutine chebpw_eval106(chebscheme,coefs1,coefs2,n,xs,vals1,vals2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex    :: coefs1(:), coefs2(:),vals1(:),vals2(:)
double precision    :: xs(:)
!
!  Evaluate two complex-valued expansions of the form (1) at a collection of points
!  given their vectors (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs1 - the coefficients vector for the first expansion
!    coefs2 - the coefficients vector for the second expansion
!    n - the number of input points
!    xs - the points in [a,b] at which to evaluate the functions
!
!  Output parameters:
!    vals1 - the values of the first expansion
!    vals2 - the values of the second expansion

do i=1,n
call chebpw_eval6(chebscheme,coefs1,coefs2,xs(i),vals1(i),vals2(i))
end do

end subroutine



subroutine chebpw_evalder101(chebscheme,coefs,n,xs,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: xs(:)
double precision    :: coefs(:), vals(:), ders(:)
!
!  Evalderuate a real-valued expansion of the form (1) at a collection of points
!  given the vector (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - the vector of coefficients (3)
!    n - the number of points at which to evalderuate the expansion
!    xs - the list of evalderuation points
!    x - the point in [a,b] at which to evalderuate the function
!
!  Output parameters:
!    vals - the values of the function at the specified points
! 

do i=1,n
call  chebpw_evalder1(chebscheme,coefs,xs(i),vals(i),ders(i))
end do

end subroutine


subroutine chebpw_evalder102(chebscheme,nfuns,coefs,n,xs,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs(:,:), vals(:,:), ders(:,:)
double precision    :: xs(:)
!
!  Evalderuate several real-valued expansions of the form (1) at a collection of
!  points given their vectors (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - an array whose jth vector is the vector of coefficients (3)
!     for the jth expansion
!    n - the number of points at which to evalderuate the expansions
!    xs - the list of evalderuation points
!
!  Output parameters:
!    vals - the ith row of this (n,nfuns) array gives the values of the
!      expansions at the ith input point
! 

do i=1,n
call chebpw_evalder2(chebscheme,nfuns,coefs,xs(i),vals(i,:),ders(i,:))
end do


end subroutine

subroutine chebpw_evalder103(chebscheme,coefs1,coefs2,n,xs,vals1,vals2,ders1,ders2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: coefs1(:), coefs2(:),vals1(:),vals2(:),ders1(:),ders2(:)
double precision    :: xs(:)
!
!  Evalderuate two real-valued expansions of the form (1) at a collection of points
!  given their vectors (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs1 - the coefficients vector for the first expansion
!    coefs2 - the coefficients vector for the second expansion
!    n - the number of input points
!    xs - the points in [a,b] at which to evalderuate the functions
!
!  Output parameters:
!    vals1 - the values of the first expansion
!    vals2 - the values of the second expansion

do i=1,n
call chebpw_evalder3(chebscheme,coefs1,coefs2,xs(i),vals1(i),vals2(i),ders1(i),ders2(i))
end do

end subroutine



subroutine chebpw_evalder104(chebscheme,coefs,n,xs,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double precision    :: xs(:)
double complex      :: coefs(:), vals(:), ders(:)
!
!  Evalderuate a complex-valued expansion of the form (1) at a collection of points
!  given the vector (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - the vector of coefficients (3)
!    n - the number of points at which to evalderuate the expansion
!    xs - the list of evalderuation points
!    x - the point in [a,b] at which to evalderuate the function
!
!  Output parameters:
!    vals - the values of the function at the specified points
! 

do i=1,n
call  chebpw_evalder4(chebscheme,coefs,xs(i),vals(i),ders(i))
end do

end subroutine


subroutine chebpw_evalder105(chebscheme,nfuns,coefs,n,xs,vals,ders)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex      :: coefs(:,:), vals(:,:), ders(:,:)
double precision    :: xs(:)
!
!  Evalderuate several complex-valued expansions of the form (1) at a collection of
!  points given their vectors (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs - an array whose jth vector is the vector of coefficients (3)
!     for the jth expansion
!    n - the number of points at which to evalderuate the expansions
!    xs - the list of evalderuation points
!
!  Output parameters:
!    vals - the ith row of this (n,nfuns) array gives the values of the
!      expansions at the ith input point
! 

do i=1,n
call chebpw_evalder5(chebscheme,nfuns,coefs,xs(i),vals(i,:),ders(i,:))
end do


end subroutine

subroutine chebpw_evalder106(chebscheme,coefs1,coefs2,n,xs,vals1,vals2,ders1,ders2)
implicit double precision (a-h,o-z)
type(chebpw_scheme) :: chebscheme
double complex    :: coefs1(:), coefs2(:),vals1(:),vals2(:),ders1(:),ders2(:)
double precision    :: xs(:)
!
!  Evalderuate two complex-valued expansions of the form (1) at a collection of points
!  given their vectors (3) of coefficients.
!
!  Input parameters:
!    chebscheme - structure describing a Chebyshev discretization scheme
!    coefs1 - the coefficients vector for the first expansion
!    coefs2 - the coefficients vector for the second expansion
!    n - the number of input points
!    xs - the points in [a,b] at which to evalderuate the functions
!
!  Output parameters:
!    vals1 - the values of the first expansion
!    vals2 - the values of the second expansion

do i=1,n
call chebpw_evalder6(chebscheme,coefs1,coefs2,xs(i),vals1(i),vals2(i),ders1(i),ders2(i))
end do

end subroutine

end module
