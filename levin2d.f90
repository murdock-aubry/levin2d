!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!  This module contains code for efficiently evaluating integrals of the form
! 
!        b     d
!    \int  \int  f(x, y) exp (i \omega g(x, y) ) dx dy                                      (1)
!        a     c
! 
! where g: \R^2 \to \R and \omega \in \R, via the "Two dimensional adaptive Levin method". The 
! functions f(x) and g(x) are specified by the user via an external subroutine or piecewise 
! Chebyshev expansions. 
!   
!   init_levin2d - initialize the code -- this entails populating a structure which
!      is then passed to the other codes in this file
! 
!   levin2d_adap_exp - adaptively compute integrals of the form (1)
! 
!   levin2d_onesquare_exp - Estimate the value of (1) over a single subrectangle 
! 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module levin2d

    use utils
    use linalg0
    use chebyshev
    use levin
    use chebpw
    
    interface
    
    subroutine levin2d_fun(n,xs,ys,f,g,gx,gy,par1,par2,par3)
      double precision              :: xs(:),ys(:)
      double precision              :: par1, par2, par3
      double complex                :: f(:), g(:), gx(:), gy(:)
    end subroutine
    
    end interface
    
    type        levin2d_vars_t
    type(levin_vars_t)            :: levin_vars
    integer                       :: norder
    integer                       :: ncoefs
    integer                       :: iquad
    integer                       :: ibasis
    
    integer                       :: nquad
    double precision, allocatable :: xs(:)
    double precision, allocatable :: ys(:)
    double precision, allocatable :: whts(:)
    
    integer                       :: ncheb
    double precision, allocatable :: xscheb(:)
    double precision, allocatable :: whtscheb(:)
    double precision, allocatable :: adiffcheb(:,:)
    double precision, allocatable :: acoefscheb(:,:)
    double complex, allocatable   :: diffx(:, :), diffy(:, :)
    double complex, allocatable   :: evalmat(:, :)
    double complex, allocatable   :: xleft(:, :), xright(:, :), yleft(:, :), yright(:, :)
    double complex, allocatable   :: xleft_interp(:, :), xright_interp(:, :), yleft_interp(:, :), yright_interp(:, :)

    
    end type    levin2d_vars_t
    
    contains
    
    subroutine init_levin2d(vars,norder,ncheb)
      implicit double precision (a-h,o-z)
      type(levin2d_vars_t), intent(out)             :: vars
      double precision, allocatable                 :: whts(:), xs(:)
    
      vars%norder = norder

      !
      ! Get the Chebyshev data
      !
      
      vars%ncheb  = ncheb
      call chebyshev_quad(vars%ncheb,vars%xscheb,vars%whtscheb)
      call chebyshev_coefsmatrix(vars%ncheb,vars%acoefscheb)
      call chebyshev_diffmatrix(vars%ncheb,vars%adiffcheb)
    
      call chebyshev_quad(norder + 1,xs,whts)
      nquad  = (norder + 1) ** 2
    
      allocate(vars%xs(nquad), vars%ys(nquad), vars%whts(nquad))

      vars%nquad = nquad
      idx        = 0
      do i=1,norder + 1
          do j=1,norder + 1
              idx            = idx + 1
              vars%xs(idx)   = xs(i)
              vars%ys(idx)   = xs(j)
          end do
      end do
    
      
      ncoefs = (norder + 1) * (norder + 2) / 2
    
      vars%ncoefs = ncoefs

      call levin2d_diffmat(norder,ncoefs,vars%xs,vars%ys, vars%evalmat, vars%diffx, vars%diffy)
      call levin2d_bdyval(norder,ncoefs,ncheb,vars%xscheb,vars%xleft, vars%xright, vars%yleft, vars%yright)
      
      call levin2d_bdy_interp(norder,ncoefs,ncheb,vars%xs,vars%ys,vars%xscheb,&
                vars%xleft_interp, vars%xright_interp, vars%yleft_interp, vars%yright_interp)
    
      call levin_init(vars%levin_vars,vars%ncheb, maxstack)
    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    


    subroutine levin2d_adap_exp(vars,ier,eps,fun,omega,x1,x2,y1,y2,dtime,nboxes,par1,par2,par3,val)
      implicit double precision (a-h,o-z)
      type(levin2d_vars_t)               :: vars
      procedure(levin2d_fun)             :: fun
      double complex                     :: val
      double precision                   :: omega
      !
      ! Adaptively evaluate integrals of the form (1) using the two-dimensional adaptive Levin method.
      !
      ! Input parameters:
      !   (x1, x2, y1, y2) - the rectangle over which to integrate
      !   vars             - type variable initialized by levin2d_init
      !   eps              - tollerance parameter >0
      !   fun              - external subroutine used to evaluate the functions f, g, g_x, g_y.
      !   omega            - scalar value used in the exponent of the integrand of (1)
      !   
      ! Output parameters:
      !   val              - the estimated value of the integral (1)
      !
      double precision, allocatable      :: stack(:,:)
      double complex                     :: val0, val00, val1, val2, val3, val4, ima
    
      ier      = 0
      ima      = (0.0d0,1.0d0)
      maxstack = 10000
      val      = 0
      nboxes   = 0
      ind      = 0
      dtime    = 0

      call elapsed(t1)
      call levin2d_onesquare_exp(vars,fun,omega,x1,x2,y1,y2,val0,par1,par2,par3)
      call elapsed(t2)

      dtime = dtime + (t2 - t1)
    
      allocate(stack(6,maxstack))
      nstack     = 1
      stack(1,1) = x1
      stack(2,1) = x2
      stack(3,1) = y1
      stack(4,1) = y2
      stack(5,1) = real(val0)
      stack(6,1) = imag(val0)
    
      do while (nstack > 0 )
    
        ifsplit = 0
        diff    = 0
        a0      = stack(1,nstack)
        b0      = stack(2,nstack)
        c0      = stack(3,nstack)
        d0      = stack(4,nstack)
        valr    = stack(5,nstack)
        vali    = stack(6,nstack)
        val00   = valr + ima*vali
        nstack  = nstack-1
      
      
        x0     = (a0+b0)/2
        y0     = (c0+d0)/2
        
        call elapsed(t1)
        call levin2d_onesquare_exp(vars,fun,omega,a0,x0,c0,y0,val1,par1,par2,par3)
        call levin2d_onesquare_exp(vars,fun,omega,x0,b0,c0,y0,val2,par1,par2,par3)
        call levin2d_onesquare_exp(vars,fun,omega,x0,b0,y0,d0,val3,par1,par2,par3)
        call levin2d_onesquare_exp(vars,fun,omega,a0,x0,y0,d0,val4,par1,par2,par3)
        call elapsed(t2)

        dtime = dtime + (t2 - t1)

        val0  = val1+val2+val3+val4
        diff  = abs(val00-val0)
      
        if (diff .gt. eps) ifsplit = 1

        ind    = ind+1
      
        write (*,"(I0,1X,I0,1X,5(ES25.15),2X,I1.1)")  ind,nboxes,a0,b0,c0,d0,diff,ifsplit
        write (1111,"(I0,1X,I0,1X,5(ES25.15),2X,I1.1)")  ind,nboxes,a0,b0,c0,d0,diff,ifsplit
        
        
        if (ifsplit .eq. 0) then
          nboxes = nboxes+1
          val    = val + val00
        else
          if (nstack+4 .ge. maxstack) then
            ier = 4
            return
          endif
        
          nstack          = nstack+1
          stack(1,nstack) = a0
          stack(2,nstack) = x0
          stack(3,nstack) = c0
          stack(4,nstack) = y0
          stack(5,nstack) = real(val1)
          stack(6,nstack) = imag(val1)
        
          nstack          = nstack+1
          stack(1,nstack) = x0
          stack(2,nstack) = b0
          stack(3,nstack) = c0
          stack(4,nstack) = y0
          stack(5,nstack) = real(val2)
          stack(6,nstack) = imag(val2)
        
          nstack          = nstack+1
          stack(1,nstack) = x0
          stack(2,nstack) = b0
          stack(3,nstack) = y0
          stack(4,nstack) = d0
          stack(5,nstack) = real(val3)
          stack(6,nstack) = imag(val3)
        
          nstack          = nstack+1
          stack(1,nstack) = a0
          stack(2,nstack) = x0
          stack(3,nstack) = y0
          stack(4,nstack) = d0
          stack(5,nstack) = real(val4)
          stack(6,nstack) = imag(val4)

        endif
      end do
    
    end subroutine
    
    
    subroutine levin2d_onesquare_exp(vars,fun,omega,x1,x2,y1,y2,val,par1,par2,par3)
      implicit double precision (a-h,o-z)
      type(levin2d_vars_t)           :: vars
      procedure(levin2d_fun)         :: fun
      double complex                 :: val
      double precision               :: omega
      !
      ! Estimate the value of integrals of the form (1) over a single subrectangle
      !
      ! Input parameters:
      !   (x1, x2, y1, y2) - the rectangle over which to integrate
      !   vars             - type variable initialized by levin2d_init
      !   eps              - tollerance parameter >0
      !   fun              - external subroutine used to evaluate the functions f, g, g_x, g_y.
      !   omega            - scalar value used in the exponent of the integrand of (1)
      !   
      ! Output parameters:
      !   val              - the estimated value of the integral (1) over the single subrectangle.
      !
      double complex                 :: ima, val1, val2
      double precision, allocatable  :: xs(:), ys(:)
      double precision, allocatable  :: xsbdy(:), ysbdy(:)
      double complex, allocatable    :: f(:), g(:), gx(:), gy(:)
      double complex, allocatable    :: amatr(:,:),p(:)
      double complex, allocatable    :: gbdy(:)
      double complex, allocatable    :: pbdy(:)
    
      ima  = (0.0d0,1.0d0)
      eps  = 1.0d-13
    
      norder = vars%norder
      ncoefs = vars%ncoefs
      nquad  = vars%nquad
      ncheb  = vars%ncheb
      
      allocate(xs(nquad),ys(nquad))
      allocate(amatr(nquad,ncoefs),p(ncoefs))
      allocate(f(nquad),g(nquad),gx(nquad),gy(nquad))
      allocate(xsbdy(ncheb),ysbdy(ncheb),pbdy(ncheb))
      allocate(gbdy(ncheb))
    
      xs = (x2-x1)/2 * vars%xs + (x2+x1)/2
      ys = (y2-y1)/2 * vars%ys + (y2+y1)/2
    
      call fun(nquad,xs,ys,f,g,gx,gy,par1,par2,par3)

      ! Decide direction of vector-field solution based on magntidue of \grad(g)
      idir = 0
      if (minval(abs(gx)) < minval(abs(gy))) then
        idir = 1
      end if

      ! write (*,"(I0,1X,I0,1X,5(ES25.15),2X,I1.1)")  ind,nboxes,a0,b0,c0,d0,idir
      write (2222,"(4(ES25.15),2X,I1.1, I0)")  x1,x2,y1,y2,idir



      !!!!!!!!!! Construct amatr !!!!!!!!!!
      amatr = 0

      if (idir == 0) then
        amatr = vars%diffx * 2 / (x2 - x1)
        do i=1,nquad
          amatr(i,:) = amatr(i,:) + ima * omega * gx(i) * vars%evalmat(i,:)
        end do

      else
        amatr = vars%diffy * 2 / (y2 - y1)
        do i=1,nquad
          amatr(i,:) = amatr(i,:) + ima * omega * gy(i) * vars%evalmat(i,:)
        end do
      end if

      call linalg0_qrsolve(eps,nquad,ncoefs,amatr,p,f)

      xsbdy = 0
      ysbdy = 0
      gbdy = 0
      pbdy = 0

      if (idir == 0) then 
        ! Evaluate the integral on the left boundary
        xsbdy   = x1
        ysbdy   = (y2-y1)/2 * vars%xscheb + (y2+y1)/2
        gbdy = ima * omega * matmul(vars%yleft_interp, g)
        pbdy = - matmul(vars%yleft, p)
        call levin_adapexp4(vars%levin_vars, ier, eps, y1, y2, gbdy, pbdy, val1)

        ! Evaluate the integral on the right boundary
        xsbdy   = x2
        ysbdy   = (y2-y1)/2 * vars%xscheb + (y2+y1)/2
        gbdy = ima * omega * matmul(vars%yright_interp, g)
        pbdy = matmul(vars%yright, p)
        call levin_adapexp4(vars%levin_vars, ier, eps, y1, y2, gbdy, pbdy, val2)

      else
        ! Evaluate the integral on the top boundary
        xsbdy   = (x2-x1)/2 * vars%xscheb + (x2+x1)/2
        ysbdy   = y2   
        gbdy = ima * omega * matmul(vars%xright_interp, g)
        pbdy = matmul(vars%xright, p)
        call levin_adapexp4(vars%levin_vars, ier, eps, x1, x2, gbdy, pbdy, val1)

        ! Evaluate the integral on the bottom boundary
        xsbdy   = (x2-x1)/2 * vars%xscheb + (x2+x1)/2
        ysbdy   = y1    
        gbdy = ima * omega * matmul(vars%xleft_interp, g)
        pbdy = -matmul(vars%xleft, p)
        call levin_adapexp4(vars%levin_vars, ier, eps, x1, x2, gbdy, pbdy, val2)

      end if
      
      val = val1 + val2

    end subroutine


    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,vals)
      implicit double precision (a-h,o-z)
      double precision                  :: vals(:)
      double precision                  :: polsx(0:norder),    polsy(0:norder)
    
      xx   = (2*x - (x2+x1) ) /(x2-x1)
      yy   = (2*y - (y2+y1) ) /(y2-y1)
  
      call chebs(norder+1,xx,polsx)
      call chebs(norder+1,yy,polsy)
  
      idx = 0
      do nn=0,norder
        do i=0,nn
          j=nn-i
          idx = idx+1
          vals(idx) = polsx(i)*polsy(j)
        end do
      end do  
    end subroutine
    
    
    subroutine levin2d_evalbasis2(norder,x1,x2,y1,y2,x,y,vals,derxs,derys)
      implicit double precision (a-h,o-z)
      double precision                  :: vals(:), derxs(:), derys(:)
      double precision                  :: polsx(0:norder),    polsy(0:norder)
      double precision                  :: poldersx(0:norder), poldersy(0:norder)
    
      xx   = (2*x - (x2+x1) ) /(x2-x1)
      yy   = (2*y - (y2+y1) ) /(y2-y1)
    
      call chebders(norder+1,xx,polsx,poldersx)
      call chebders(norder+1,yy,polsy,poldersy)
      poldersx=poldersx*2/(x2-x1)
      poldersy=poldersy*2/(y2-y1)
      
      idx = 0
      do nn=0,norder
        do i=0,nn
          j=nn-i
          idx = idx+1
      
          vals(idx)         =  polsx(i)*polsy(j)
          derxs(idx)        =  poldersx(i)*polsy(j)   
          derys(idx)        =  polsx(i)*poldersy(j)
        end do
      end do
    
    end subroutine
        
    
    subroutine levin2d_diffmat(norder, ncoefs, xs, ys, evalmat, diffx, diffy)
      ! A subroutine which returns the spectral differentiation matrices for bivariate Chebyshev expansions of the form 
      implicit double precision (a-h, o-z)
      double precision                      :: xs(:), ys(:)
      integer                               :: norder, ncoefs
      double complex, allocatable           :: evalmat(:, :), diffx(:, :), diffy(:, :)
      double precision, allocatable         :: v1(:), v2(:), v3(:)
    
      nquad = size(xs)

      x1 = -1.0d0
      x2 = 1.0d0
      y1 = -1.0d0
      y2 = 1.0d0
    
      allocate(v1(ncoefs),v2(ncoefs), v3(ncoefs))
      allocate(evalmat(nquad, ncoefs), diffx(nquad, ncoefs), diffy(nquad, ncoefs))
      
      diffx = 0
      diffy = 0
    
      do l=1,nquad
        x = xs(l)
        y = ys(l)
        idx = 0 
    
        call levin2d_evalbasis2(norder,x1,x2,y1,y2,x,y,v1,v2,v3)
    
        evalmat(l, :) = v1
        diffx(l,:)    = v2
        diffy(l,:)    = v3
      end do
    end subroutine


    subroutine levin2d_evalmat(norder, ncoefs, xs, ys, evalmat)
      implicit double precision (a-h, o-z)
    
      double precision                      :: xs(:), ys(:)
      integer                               :: norder, ncoefs
      double complex, allocatable           :: evalmat(:, :), diffx(:, :), diffy(:, :)
      double precision, allocatable         :: v1(:), v2(:), v3(:)
    
      nquad = size(xs)

      x1 = -1.0d0
      x2 = 1.0d0
      y1 = -1.0d0
      y2 = 1.0d0
    
      allocate(v1(ncoefs),v2(ncoefs), v3(ncoefs))
      allocate(evalmat(nquad, ncoefs), diffx(nquad, ncoefs), diffy(nquad, ncoefs))
    
      
      diffx = 0
      diffy = 0
    
      do l=1,nquad
        x = xs(l)
        y = ys(l)
        idx = 0 
    
    
        call levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,v1)
    
        evalmat(l, :) = v1
    
      end do
    end subroutine
    
    
    subroutine levin2d_bdyval(norder, ncoefs, ncheb, xsbdy, xleft, xright, yleft, yright)
      implicit double precision (a-h, o-z)
      double precision, allocatable               :: v1(:)
      double complex, allocatable                 :: xleft(:, :), xright(:, :), yleft(:, :), yright(:, :)
      double precision                            :: xsbdy(:)
    
      x1 = -1.0d0
      x2 = 1.0d0
      y1 = -1.0d0
      y2 = 1.0d0
    
      x0 = (x1+x2)/2
      y0 = (y1+y2)/2
    
      r  = (x2-x1)/2
      s  = (y2-y1)/2
    
      allocate(xleft(ncheb, ncoefs), xright(ncheb, ncoefs), yleft(ncheb, ncoefs), yright(ncheb, ncoefs))
      allocate(v1(ncoefs))
    
      xleft = 0
      xright = 0
      yleft = 0
      yright = 0
    
      do l = 1, ncheb
        ! xleft
        x       = xsbdy(l)
        y       = y1
        call levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,v1)
        xleft(l, :) = v1
    
    
        ! xright
        x       = xsbdy(l)
        y       = y2
        call levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,v1)
        xright(l, :) = v1
    
        ! yleft
        x       = x1
        y       = xsbdy(l)
        call levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,v1)
        yleft(l, :) = v1
    
    
        ! yright
        x       = x2
        y       = xsbdy(l)
        call levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,v1)
        yright(l, :) = v1
    
      end do
    end subroutine 


    subroutine levin2d_coefmat(norder, ncoefs, xs, ys, coefmat)
      implicit double precision (a-h, o-z)
      double precision                        :: xs(:), ys(:)
      double complex, allocatable             :: coefmat(:, :), evalmat(:, :)

      eps = epsilon(0.0d0) * 1000

      nquad = size(xs)

      allocate(coefmat(ncoefs, nquad))

      call levin2d_evalmat(norder, ncoefs, xs, ys, evalmat)

      call linalg0_qrtinvert_c(eps, nquad, ncoefs, evalmat, coefmat)

    end subroutine
  

    subroutine levin2d_bdy_interp(norder, ncoefs, ncheb, xs, ys, xsbdy, xleft, xright, yleft, yright)
      implicit double precision (a-h, o-z)
      double precision                            :: xs(:), ys(:)

      double precision, allocatable               :: v1(:)
      double complex, allocatable                 :: xleft(:, :), xright(:, :), yleft(:, :), yright(:, :)
      double complex, allocatable                 :: aleft(:, :), aright(:, :), bleft(:, :), bright(:, :)
      double precision                            :: xsbdy(:)
      double complex, allocatable                 :: coefmat(:, :)

      x1 = -1.0d0
      x2 = 1.0d0
      y1 = -1.0d0
      y2 = 1.0d0
    
      x0 = (x1+x2)/2
      y0 = (y1+y2)/2
    
      r  = (x2-x1)/2
      s  = (y2-y1)/2
    
      allocate(aleft(ncheb, ncoefs), aright(ncheb, ncoefs), bleft(ncheb, ncoefs), bright(ncheb, ncoefs))
      allocate(xleft(ncheb, nquad), xright(ncheb, nquad), yleft(ncheb, nquad), yright(ncheb, nquad))

      allocate(v1(ncoefs))

      call levin2d_coefmat(norder, ncoefs, xs, ys, coefmat)
    
      xleft = 0
      xright = 0
      yleft = 0
      yright = 0

      aleft = 0
      aright = 0
      bleft = 0
      bright = 0
    
      do l = 1, ncheb
        ! xleft
        x       = xsbdy(l)
        y       = y1
        call levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,v1)
        aleft(l, :) = v1
    
    
        ! xright
        x       = xsbdy(l)
        y       = y2
        call levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,v1)
        aright(l, :) = v1
    
        ! yleft
        x       = x1
        y       = xsbdy(l)
        call levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,v1)
        bleft(l, :) = v1
    
    
        ! yright
        x       = x2
        y       = xsbdy(l)
        call levin2d_evalbasis1(norder,x1,x2,y1,y2,x,y,v1)
        bright(l, :) = v1
    
      end do


      xleft  = matmul(aleft, coefmat)
      xright = matmul(aright, coefmat)
      yleft  = matmul(bleft, coefmat)
      yright = matmul(bright, coefmat)
    
    end subroutine 
    

    end module
    