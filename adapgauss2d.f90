module adapgauss2d

use utils
use legendre
use iso_c_binding

type       adapgauss2d_vars_t
    integer                                    :: maxstack
    integer                                    :: nlege
    double precision, allocatable              :: xslege(:), yslege(:), whtslege(:)
    double precision, allocatable              :: xs(:,:), whts(:)
end type   adapgauss2d_vars_t

interface        

subroutine adapgauss2d_fun(n,ts,vals,par1,par2,par3)
    use iso_c_binding
    implicit double precision (a-h,o-z)
    double precision :: ts(:, :)
    double complex   :: vals(:)
    double complex   :: alpha
    type(c_ptr)      :: userptr
end subroutine

end interface

contains 

subroutine adapgauss2d_init(vars,maxstack0,nlege0)
    implicit double precision (a-h,o-z)
    type(adapgauss2d_vars_t)             :: vars
    integer, optional                    :: maxstack0
    integer, optional                    :: nlege0
    !
    !  Iniitalize the structure containing any data needed by the other
    !  procedures in this module. 
    !
    !  If this routine is called with no arguments other than vars,
    !  reasonable defaults are chosen.
    !
    !  Input parameters:
    !   nlege0 - the order of the Gauss-Legendre quadrature to use
    !   maxstack0 - the size of the stack used to conduct the adaptive
    !     integration procedure
    !
    !  Output parameters:
    !    vars - the structure which holds all of the necessary data
    !      for the other routines in this module
    !
    
    if (.not. present(maxstack0)) then
        maxstack = 10000
        nlege    = 12
    else
        maxstack = maxstack0
        nlege    = nlege0
    endif
    
    vars%nlege     = nlege
    vars%maxstack  = maxstack
    
    call legendre_quad(nlege,vars%xslege,vars%whtslege)

    allocate(vars%xs(nlege**2, 2), vars%whts(nlege**2))

    idx = 1
    do i = 1, nlege
        do j = 1, nlege
            vars%xs(idx, 1) = vars%xslege(i)
            vars%xs(idx, 2) = vars%xslege(j)
            vars%whts(idx)  = vars%whtslege(i) * vars%whtslege(j)

            idx = idx + 1
        end do 
    end do
    
end subroutine


subroutine adapgauss2d_int(vars,ier,eps,a,b,c,d,fun,par1,par2,par3,val)
    implicit double precision (a-h,o-z)
    type(adapgauss2d_vars_t)               :: vars
    procedure(adapgauss2d_fun)             :: fun
    double complex                         :: ima, alpha
    !
    !  Evaluate an integral of the form (1) with a complex-valued integrand
    !  using a piecewise Gauss-Legendre method.
    !
    !  Input parameters:
    !    vars - the structure populated by the initialization routine
    !    eps - the desired precision for the calculation
    !    (a,b) - the interval over which the integral is to be computed
    !    fun - an external subroutine conforming to the adapgauss_fun_2d interface
    !    par? - user-supplied parameters which are passed to fun
    !
    !  Ouptut parameters:
    !    ier - an error return code;
    !       ier = 0  indicates successful execution
    !       ier = 8  means too many intervals were required to evaluate the integral
    !              
    !
    !    val - the value of the integral as calculated by this routine
    !

    double precision, allocatable       :: stack(:,:), xs(:, :), whts(:)
    double complex, allocatable         :: vals(:)
    double complex                      :: val, val0, val00, val1, val2, val3, val4

    ier = 0
    nlege    = vars%nlege
    maxstack = vars%maxstack
    ima = (0.0d0, 1.0d0)

    allocate(stack(6,maxstack),xs(nlege ** 2, 2), whts(nlege ** 2), vals(nlege ** 2))

    xs(:, 1) = min(b,max(a,(b-a)/2 * vars%xs(:, 1) + (b+a)/2))
    xs(:, 2) = min(d,max(c,(d-c)/2 * vars%xs(:, 2) + (d+c)/2))
    whts = vars%whts * (b-a)*(d-c)/4

    call fun(nlege**2,xs,vals,par1,par2,par3)
    val0 = sum(vals*whts)

    

    nstack     = 1
    stack(1,1) = a
    stack(2,1) = b
    stack(3,1) = c
    stack(4,1) = d
    stack(5,1) = real(val0)
    stack(6,1) = imag(val0)
    val        = 0

    idx = 0

    do while (nstack > 0) 
        idx   = idx + 1
        a0    = stack(1,nstack)
        b0    = stack(2,nstack)
        c0    = stack(3,nstack)
        d0    = stack(4,nstack)
        val0r = stack(5,nstack)
        val0i = stack(6,nstack)
        val0  = val0r + ima*val0i

        

        nstack = nstack-1

        x0 = (a0+b0)/2
        y0 = (c0+d0)/2

        
        

        ! Bottom left
        xs(:, 1) = min(x0,max(a0,(x0-a0)/2 * vars%xs(:, 1) + (x0+a0)/2))
        xs(:, 2) = min(y0,max(c0,(y0-c0)/2 * vars%xs(:, 2) + (y0+c0)/2))
        whts     = vars%whts * (x0-a0) * (y0-c0) / 4
        call fun(nlege**2,xs,vals,par1,par2,par3)
        val1 = sum(vals*whts)
        
        

        ! Bottom right
        xs(:, 1) = min(b0,max(x0,(b0-x0)/2 * vars%xs(:, 1) + (b0+x0)/2))
        xs(:, 2) = min(y0,max(c0,(y0-c0)/2 * vars%xs(:, 2) + (y0+c0)/2))
        whts     = vars%whts * (b0-x0) * (y0-c0) / 4
        call fun(nlege**2,xs,vals,par1,par2,par3)
        val2 = sum(vals*whts)


        ! Top left
        xs(:, 1) = min(x0,max(a0,(x0-a0)/2 * vars%xs(:, 1) + (x0+a0)/2))
        xs(:, 2) = min(d0,max(y0,(d0-y0)/2 * vars%xs(:, 2) + (d0+y0)/2))
        whts     = vars%whts * (x0-a0) * (d0-y0) / 4
        call fun(nlege**2,xs,vals,par1,par2,par3)
        val3 = sum(vals*whts)


        ! Top right
        xs(:, 1) = min(b0,max(x0,(b0-x0)/2 * vars%xs(:, 1) + (b0+x0)/2))
        xs(:, 2) = min(d0,max(y0,(d0-y0)/2 * vars%xs(:, 2) + (d0+y0)/2))
        whts     = vars%whts * (b0-x0) * (d0-y0) / 4
        call fun(nlege**2,xs,vals,par1,par2,par3)
        val4 = sum(vals*whts)


        val00 = val1 + val2 + val3 + val4 
    
        
        diff = abs(val0-val00)



        
        ifsplit = 0


        if (diff .le. eps) ifplit = 1
        print *, idx, ifplit, a0, b0, c0, d0, diff

        if (diff .le. eps) then
            val = val  + val0
            ifsplit = 1
        else
            if (nstack+4 .gt. maxstack) then
                ier = 8
                return
            endif


            ! bottom left
            nstack          = nstack+1
            stack(1,nstack) = a0
            stack(2,nstack) = x0
            stack(3,nstack) = c0
            stack(4,nstack) = y0
            stack(5,nstack) = real(val1)
            stack(6,nstack) = imag(val1)


            ! bottom right
            nstack          = nstack+1
            stack(1,nstack) = x0
            stack(2,nstack) = b0
            stack(3,nstack) = c0
            stack(4,nstack) = y0
            stack(5,nstack) = real(val2)
            stack(6,nstack) = imag(val2)

            ! top left
            nstack          = nstack+1
            stack(1,nstack) = a0
            stack(2,nstack) = x0
            stack(3,nstack) = y0
            stack(4,nstack) = d0
            stack(5,nstack) = real(val3)
            stack(6,nstack) = imag(val3)

            ! top right
            nstack          = nstack+1
            stack(1,nstack) = x0
            stack(2,nstack) = b0
            stack(3,nstack) = y0
            stack(4,nstack) = d0
            stack(5,nstack) = real(val4)
            stack(6,nstack) = imag(val4)


        endif

    end do



end subroutine

end module