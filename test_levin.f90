module        test_levin_functions

use utils
use adapgauss
use levin

contains

subroutine testfun(n,ts,gs,fs,dlambda,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double complex              :: gs(:), fs(:)
double complex              :: ima
data ima / (0.0d0,1.0d0) /
data pi  / 3.14159265358979323846264338327950288d0 /

gs = ima*dlambda*ts**2 - log(ts)
fs = 1

end subroutine

subroutine testfun2(n,ts,vals,dlambda,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double complex              :: vals(:)
double complex              :: ima
double complex              :: gs(n), fs(n)
data ima / (0.0d0,1.0d0) /
data pi  / 3.14159265358979323846264338327950288d0 /
call testfun(n,ts,gs,fs,dlambda,par2,par3)
vals  = fs * exp(gs)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine testfun_diff(n,ts,gs,gders,fs,dlambda,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double complex              :: gs(:), gders(:), fs(:)
double complex              :: ima
data ima / (0.0d0,1.0d0) /
data pi  / 3.14159265358979323846264338327950288d0 /


gs    = ima*dlambda*ts**2
gders = ima*dlambda*2*ts
fs    = ts**2

end subroutine

subroutine testfun2_diff(n,ts,vals,dlambda,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double complex              :: vals(:)
double complex              :: ima
double complex              :: gs(n), fs(n), gders(n)
data ima / (0.0d0,1.0d0) /
data pi  / 3.14159265358979323846264338327950288d0 /
call testfun_diff(n,ts,gs,gders,fs,dlambda,par2,par3)
vals  = fs * exp(gs)
end subroutine



subroutine testfun_cos(n,ts,gs,fs,dlambda,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double precision            :: gs(:), fs(:)
data pi  / 3.14159265358979323846264338327950288d0 /

gs = dlambda*ts
fs = 1

end subroutine

subroutine testfun2_cos(n,ts,vals,dlambda,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double precision            :: vals(:)
double precision            :: gs(n), fs(n)
data pi  / 3.14159265358979323846264338327950288d0 /

call testfun_cos(n,ts,gs,fs,dlambda,par2,par3)
vals  = fs * cos(gs)
end subroutine


subroutine testfun_sin(n,ts,gs,fs,dlambda,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double precision            :: gs(:), fs(:)
data pi  / 3.14159265358979323846264338327950288d0 /

gs = dlambda*ts
fs = 1.0d0/ts

end subroutine

subroutine testfun2_sin(n,ts,vals,dlambda,par2,par3)
implicit double precision (a-h,o-z)
double precision            :: ts(:)
double precision            :: vals(:)
double precision            :: gs(n), fs(n)
data pi  / 3.14159265358979323846264338327950288d0 /

call testfun_sin(n,ts,gs,fs,dlambda,par2,par3)
vals  = fs * sin(gs)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Test the adaptive integration code by comparison with adapgauss
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine levin_test_adapexp()
implicit double precision (a-h,o-z)

type(levin_vars_t)           :: levin_vars
type(adapgauss_vars_t)       :: adapgauss_vars
double complex               :: ima, val, val0, val1, val2, val01, val02

!
!  Initialize the Levin code
!

ima      = (0.0d0,1.0d0)
pi       = acos(-1.0d0)
eps0     = epsilon(0.0d0)
k        = 12
maxstack = 1000
call levin_init(levin_vars,k,maxstack)


!
!  Initialize adapgauss code
!
maxstack = 1000000
nlege    = 30
call  adapgauss_init(adapgauss_vars,maxstack,nlege)

a       =  1.0d-12
b       =  1.0d0
dlambda =  200.0d0
eps     =  1.0d-13

call elapsed(t1)
call levin_adapexp(levin_vars,ier,eps,a,b,testfun,dlambda,par2,par3,val1)
call elapsed(t2)
tlevin = t2-t1


if (ier .ne. 0) then
print *,"after levin_adap, ier = ",ier
stop
endif

call elapsed(t1)
call adapgauss_int(adapgauss_vars,ier,eps,a,b,testfun2,dlambda,par2,par3,val01)
call elapsed(t2)
tgauss = t2-t1

a       = -1.0d0
b       = -1.0d-12


call levin_adapexp(levin_vars,ier,eps,a,b,testfun,dlambda,par2,par3,val2)
call adapgauss_int(adapgauss_vars,ier,eps,a,b,testfun2,dlambda,par2,par3,val02)

val = val1+val2
val0 = val01+val02

derr = abs(val-val0)


write (*,"(A,1(D24.16,1X))") "a             = ",a
write (*,"(A,1(D24.16,1X))") "b             = ",b
write (*,"(A,1(D24.16,1X))") "dlambda       = ",dlambda
write (*,"(A,D24.16,' , ',D24.16)") "levin val     = ",val
write (*,"(A,D24.16,' , ',D24.16)") "adapgauss val = ",val0
write (*,"(A,1(D24.16,1X))") "abs error     = ",derr
write (*,"(A,1(D24.16,1X))") "time levin    = ",tlevin
write (*,"(A,1(D24.16,1X))") "time gauss    = ",tgauss
write (*,"(A,1(D24.16,1X))") "ratio         = ",tgauss/tlevin
print *,""
print *,""


end subroutine



subroutine levin_test_adapexp_diff()
implicit double precision (a-h,o-z)

type(levin_vars_t)           :: levin_vars
type(adapgauss_vars_t)       :: adapgauss_vars

double complex               :: ima, val, val0

!
!  Initialize the Levin code
!

ima      = (0.0d0,1.0d0)
pi       = acos(-1.0d0)
eps0     = epsilon(0.0d0)
k        = 12
maxstack = 1000
call levin_init(levin_vars,k,maxstack)

!
!  Initialize adapgauss code
!
maxstack = 1000000
nlege    = 30
call  adapgauss_init(adapgauss_vars,maxstack,nlege)

a       = -1.0d0
b       =  1.0d0
dlambda =  100.0d0
eps     =  1.0d-13

call elapsed(t1)
call levin_adapexp_diff1(levin_vars,ier,eps,a,b,testfun_diff,dlambda,par2,par3,val)
call elapsed(t2)
tlevin = t2-t1

if (ier .ne. 0) then
print *,"after levin_adap, ier = ",ier
stop
endif

call elapsed(t1)
call adapgauss_int(adapgauss_vars,ier,eps,a,b,testfun2_diff,dlambda,par2,par3,val0)
call elapsed(t2)
tgauss = t2-t1

derr   = abs(val-val0)

write (*,"(A,1(D24.16,1X))") "a             = ",a
write (*,"(A,1(D24.16,1X))") "b             = ",b
write (*,"(A,1(D24.16,1X))") "dlambda       = ",dlambda
write (*,"(A,D24.16,' , ',D24.16)") "levin val     = ",val
write (*,"(A,D24.16,' , ',D24.16)") "adapgauss val = ",val0
write (*,"(A,1(D24.16,1X))") "abs error     = ",derr
write (*,"(A,1(D24.16,1X))") "time levin    = ",tlevin
write (*,"(A,1(D24.16,1X))") "time gauss    = ",tgauss
write (*,"(A,1(D24.16,1X))") "ratio         = ",tgauss/tlevin
print *,""
print *,""

end subroutine


subroutine levin_test_adapcos()
implicit double precision (a-h,o-z)

type(levin_vars_t)           :: levin_vars
type(adapgauss_vars_t)       :: adapgauss_vars

double precision             :: ima, val, val0

!
!  Initialize the Levin code
!

ima      = (0.0d0,1.0d0)
pi       = acos(-1.0d0)
eps0     = epsilon(0.0d0)
k        = 12
maxstack = 1000
call levin_init(levin_vars,k,maxstack)

!
!  Initialize adapgauss code
!
maxstack = 1000000
nlege    = 30
call  adapgauss_init(adapgauss_vars,maxstack,nlege)
!


a       = -1.0d0
b       =  1.0d0
dlambda =  1000.0d0
eps     =  1.0d-14

call elapsed(t1)
call levin_adapcos(levin_vars,ier,eps,a,b,testfun_cos,dlambda,par2,par3,val)
call elapsed(t2)
tlevin = t2-t1

if (ier .ne. 0) then
print *,"after levin_adap, ier = ",ier
stop
endif

call elapsed(t1)
call adapgauss_int(adapgauss_vars,ier,eps,a,b,testfun2_cos,dlambda,par2,par3,val0)
call elapsed(t2)
tgauss = t2-t1

derr = abs(val-val0)


write (*,"(A,1(D24.16,1X))") "a             = ",a
write (*,"(A,1(D24.16,1X))") "b             = ",b
write (*,"(A,1(D24.16,1X))") "dlambda       = ",dlambda
write (*,"(A,D24.16)") "levin val     = ",val
write (*,"(A,D24.16)") "adapgauss val = ",val0
write (*,"(A,1(D24.16,1X))") "abs error     = ",derr
write (*,"(A,1(D24.16,1X))") "time levin    = ",tlevin
write (*,"(A,1(D24.16,1X))") "time gauss    = ",tgauss
write (*,"(A,1(D24.16,1X))") "ratio         = ",tgauss/tlevin
print *,""
print *,""


end subroutine


subroutine levin_test_adapsin()
implicit double precision (a-h,o-z)

type(levin_vars_t)           :: levin_vars
type(adapgauss_vars_t)       :: adapgauss_vars

double precision             :: ima, val, val0

!
!  Initialize the Levin code
!

ima      = (0.0d0,1.0d0)
pi       = acos(-1.0d0)
eps0     = epsilon(0.0d0)
k        = 12
maxstack = 1000
call levin_init(levin_vars,k,maxstack)

!
!  Initialize adapgauss code
!
maxstack = 1000000
nlege    = 30
call  adapgauss_init(adapgauss_vars,maxstack,nlege)
!


a       =  1.0d-15
b       =  1.0d0
dlambda =  1000000.0d0
eps     =  1.0d-14

call elapsed(t1)
call levin_adapsin(levin_vars,ier,eps,a,b,testfun_sin,dlambda,par2,par3,val)
call elapsed(t2)
tlevin = t2-t1

if (ier .ne. 0) then
print *,"after levin_adap, ier = ",ier
stop
endif

call elapsed(t1)
call adapgauss_int(adapgauss_vars,ier,eps,a,b,testfun2_sin,dlambda,par2,par3,val0)
call elapsed(t2)
tgauss = t2-t1

derr = abs(val-val0)


write (*,"(A,1(D24.16,1X))") "a             = ",a
write (*,"(A,1(D24.16,1X))") "b             = ",b
write (*,"(A,1(D24.16,1X))") "dlambda       = ",dlambda
write (*,"(A,D24.16)") "levin val     = ",val
write (*,"(A,D24.16)") "adapgauss val = ",val0
write (*,"(A,1(D24.16,1X))") "abs error     = ",derr
write (*,"(A,1(D24.16,1X))") "time levin    = ",tlevin
write (*,"(A,1(D24.16,1X))") "time gauss    = ",tgauss
write (*,"(A,1(D24.16,1X))") "ratio         = ",tgauss/tlevin
print *,""
print *,""


end subroutine

end module

program       test_levin

use utils
use chebyshev
use chebpw
use levin
use test_levin_functions

implicit double precision (a-h,o-z)

call levin_test_adapexp()
! call levin_test_adapcos()
!call levin_test_adapsin()

!call levin_test_adapexp_diff()


end program
