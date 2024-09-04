module test_adapgauss2d_functions

    contains

    subroutine testfun1(n,ts,vals,omega,alpha,par3)
        implicit double precision (a-h,o-z)
        double precision         :: ts(:,:)
        double complex           :: vals(:), ima
        double precision         :: xs(n), ys(n)

        ima = (0.0d0, 1.0d0)

        xs = ts(:, 1)
        ys = ts(:, 2)

        vals = xs ** 2 * cos(xs * ys) + 2.0d0
 
    end subroutine

    subroutine testfun2(n,ts,vals,omega,alpha,par3)
        implicit double precision (a-h,o-z)
        double precision         :: ts(:,:)
        double complex           :: vals(:), ima
        double precision         :: xs(n), ys(n)

        ima = (0.0d0, 1.0d0)

        xs = ts(:, 1)
        ys = ts(:, 2)

        vals = (xs ** 2 + ys ** 2) * exp(xs) * exp(ima * omega * (xs + ys))

    end subroutine
    
end module
    
program test_adapgauss
    
    use utils
    use legendre
    use adapgauss2d
    use test_adapgauss2d_functions
    
    implicit double precision (a-h,o-z)
    type(adapgauss2d_vars_t)        :: adapgauss2d_vars
    double complex                  :: ima, zval, val0
    
    ima = (0.0d0,1.0d0)
    pi  = acos(-1.0d0)
    eps = epsilon(0.0d0)*100

    nlege    = 12
    maxstack = 10000
    call adapgauss2d_init(adapgauss2d_vars,maxstack,nlege)

    a = -1.0d0
    b = 1.0d0
    c = -1.0d0
    d = 1.0d0

    ! 
    ! TEST 1
    ! 
    
    call elapsed(t1)
    call adapgauss2d_int(adapgauss2d_vars,ier,eps,a,b,c,d,testfun1,omega,alpha,par3,zval)
    call elapsed(t2)

    val0 = 9.2046747157590271570d0 + 0.0d0 * ima

    print *, "time1 = ", t2 - t1 
    print *, "error1 = ", abs(zval - val0)


    ! 
    ! TEST 2
    ! 

    alpha   = 4
    omega   = 100.0d0 


    call elapsed(t1)
    call adapgauss2d_int(adapgauss2d_vars,ier,eps,a,b,c,d,testfun2,omega,alpha,par3,zval)
    call elapsed(t2)
    val0 = 0.00030156953765510744948d0 + 0.00040881849452999488919d0 * ima

    print *, "time2 = ", t2 - t1 
    print *, "error2 = ", abs(zval - val0)
    
end program