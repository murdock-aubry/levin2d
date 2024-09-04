module levin2d_lege_funs

contains

subroutine testfun(n,xs,ys,f,g,dg)
    implicit double precision (a-h,o-z)
    double precision              :: xs(:),ys(:)
    double complex                :: f(n), g(2, n), dg(2, 2, n)
    double complex                :: ima
    data  ima / (0.0d0,1.0d0) /

    f           = xs ** 2 + ys ** 2

    g(1, :)     = (xs ** 2 + ys ** 2)
    g(2, :)     = (xs ** 2 - ys ** 2)

    dg(1, 1, :) = 2 * xs
    dg(1, 2, :) = 2 * ys
    dg(2, 1, :) = 2 * xs
    dg(2, 2, :) = (-2) * ys


    ! g(1, :)  = ima * xs
    ! g(2, :)  = ima * ys

    ! dg(1, 1, :) = ima * 1.0d0
    ! dg(1, 2, :) = 0.0d0
    ! dg(2, 1, :) = 0.0d0
    ! dg(2, 2, :) = ima * 1.0d0


    ! f  = 1.0d0/(1.0d0+xs**2+ys**2)
end subroutine

subroutine testfun_scalar(n,xs,ys,f,g,gx,gy,par)
    implicit double precision (a-h,o-z)
    double precision              :: xs(:),ys(:)
    double complex                :: f(:), g(:), gx(:), gy(:)
    double complex                :: ima
    data  ima / (0.0d0,1.0d0) /


    ! f  = (xs)**2+(ys)**2-xs
    ! g  = (xs**2 + ys**2 + 4*ys + 4*xs)
    ! gx = (2*xs+4)
    ! gy = (2*ys+4)

    f  = xs ** 2 + ys ** 2
    g  = (xs ** 2 + ys ** 2)
    gx = 2 * xs
    gy = 2 * ys

    ! f = 2.0d0 + xs + ys
    ! g = xs + ys
    ! gx = 1.0d0
    ! gy = 1.0d0


    ! f  = 1.0d0/(1.0d0+xs**2+ys**2)
end subroutine

end module



program test_levin2d_lege_final

    use utils
    use linalg0
    use chebyshev
    use legendre
    use levin 
    use levin2d_lege
    use levin2d_final
    use levin2d_lege_funs
    use levin_lege


    implicit double precision (a-h,o-z)


    type(levin2d_lege_vars_t)            :: vars
    double complex                  :: ima, val, val0
    double precision                :: omega(2)

    double complex, allocatable     :: fvals(:), fders(:), diffx(:, :), evalmat(:, :)
    double precision, allocatable   :: xs(:), ys(:)


    norder = 11
    ncheb  = norder
    ima    = (0.0d0,1.0d0)
    eps    = 1.0d-12

    call init_levin2d_lege(vars,norder,ncheb)

    x1 = -1.0d0
    x2 =  1.0d0
    y1 = -1.0d0
    y2 =  1.0d0

    ! omega(1) = 100.0d0
    ! omega(2) = 100.0d0

    omega(1) = 1.0d0
    omega(2) = 1.0d0


    ! xs = vars%xs
    ! ys = vars%ys 
    ! diffx = vars%diffx
    ! evalmat = vars%evalmat


    ! ! Test FULL
    ! call levin2d_lege_onesquare(vars,testfun,omega,x1,x2,y1,y2,val)
    ! ! call levin2d_lege_adap(vars,ier,eps,testfun,omega,x1,x2,y1,y2,val)
    ! val0 = 1.3006147020258621279d0 + 1.7488261589112582301d0 * ima

    ! print *, "value = ", val 
    ! print *, "true = ", val0 

    ! print *, "value error = ", abs(val - val0)

    ! stop


    ! Test SCALAR

    omega_scalar = 1.0d0

    call levin2d_lege_adap_scalar(vars,ier,eps,testfun_scalar,omega_scalar,x1,x2,y1,y2,val,dtime,par)
    ! call levin2d_lege_onesquare_scalar(vars,testfun_scalar,omega_scalar,x1,x2,y1,y2,val,par)

    ! call levin2d_onesquare_scalar(vars,testfun_scalar,omega_scalar,x1,x2,y1,y2,val,par)

    ! val0 = 0.00008933547229637002d0 - 0.00002963898624439625d0*ima
    val0 = 1.4699167271423303991d0 + 1.9770516793978609568d0*ima
    ! val0 = 0.0005149977484043343598d0 - 0.0032473780162911775495d0 * ima
    print *, "value = ", val
    print *, "true = ", val0
    print *,"value error = ", abs(val-val0)
    print *, ""


end program
