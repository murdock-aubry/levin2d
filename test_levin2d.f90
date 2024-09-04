module levin2d_funs

contains

subroutine testfun(n,xs,ys,f,g,gx,gy,par1,par2,par3)
    implicit double precision (a-h,o-z)
    double precision              :: xs(:),ys(:)
    double complex                :: f(:), g(:), gx(:), gy(:)
    double complex                :: ima
    data  ima / (0.0d0,1.0d0) /


    ! f  = (xs)**2+(ys)**2-xs
    ! g  = (xs**2 + ys**2 + 4*ys + 4*xs)
    ! gx = (2*xs+4)
    ! gy = (2*ys+4)

    ! f  = xs ** 2 + ys ** 2
    ! g  = (xs ** 2 + ys ** 2)
    ! gx = 2 * xs
    ! gy = 2 * ys

    f  = 1.0d0
    g  = (1.0d0 + xs) * (1.0d0 + ys ** 2)
    gx = 1.0d0 + ys ** 2
    gy = 2 * ys * (1.0d0 + xs)


end subroutine

subroutine testfun_gauss(n,ts,vals,par1,par2,par3)
    implicit double precision (a-h,o-z)
    double precision         :: ts(:,:)
    double complex           :: vals(:), ima
    double precision         :: xs(n), ys(n)

    ima = (0.0d0, 1.0d0)

    xs = ts(:, 1)
    ys = ts(:, 2)

    vals = exp(ima * par1 * (1.0d0 + xs) * (1.0d0 + ys ** 2))

end subroutine

end module



program test_levin2d

    use utils
    use linalg0
    use chebyshev
    use chebpw
    use levin2d
    use adapgauss2d
    use levin2d_funs
    use levin


    implicit double precision (a-h,o-z)
    type(levin2d_vars_t)            :: vars
    type(adapgauss2d_vars_t)        :: vars_gauss
    double complex                  :: ima, val, val_gauss


    norder = 11 ! resulting in (norder + 1)(norder + 2) / 2 basis functions
    ncheb  = norder - 2 ! resulting in ncheb**2 quadrature nodes
    ima    = (0.0d0,1.0d0)
    eps    = 1.0d-10
    nlege    = 12
    maxstack = 10000

    call init_levin2d(vars,norder,ncheb)
    call adapgauss2d_init(vars_gauss,maxstack,nlege)


    omega = 1000.0d0

    x1 = -1.0d0
    x2 =  1.0d0
    y1 = -1.0d0
    y2 =  1.0d0
    
    call elapsed(t1)
    call adapgauss2d_int(vars_gauss, ier, eps, x1, x2, y1, y2, testfun_gauss, omega, par2, par3, val_gauss)
    call elapsed(t2)
    dtime_gauss = t2 - t1

    call levin2d_adap_exp(vars,ier,eps,testfun,omega,x1,x2,y1,y2,dtime,nboxes,par1,par2,par3,val)


    print *, "value = ", val, val_gauss
    print *, "diff = ", abs(val - val_gauss)
    print *, "time gauss/levin = ", dtime_gauss / dtime 
    print *, ""


end program
