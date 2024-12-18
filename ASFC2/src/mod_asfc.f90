
!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: may, 03 2019
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **A fast and accurate way of determining the area-scale fractal analysis complexity parameter Asfc**
!<  </span>

module asfc
use data_arch,     only : I4, R8, UN, EPS_R8
use miscellaneous, only : get_unit
use surfile,       only : read_surf, write_surf, SCALE_SURF
use bspline,       only : db2ink, db2val
use minpack,       only : lmder1
use least_squares, only : moindres_carres_lineaire
use stat_mom,      only : moment_stat, calc_moments
!$ use omp_lib
implicit none

!> {!ASFC2/src/inc_doc/Asfc.md!}

private

integer(kind=I4), parameter :: npp   = 128   !! *number of points for the asfc determination*
integer(kind=I4), parameter :: long0 =   8   !! *roughest grid*
integer(kind=I4), parameter :: larg0 =   8   !! *roughest grid*

integer(kind=I4), parameter :: kx = 3        !! *x bspline order*
integer(kind=I4), parameter :: ky = 3        !! *y bspline order*

integer(kind=I4), parameter :: nb_beta = 4   !! *number of parameters to be determined when optimizing*

integer(kind=I4), parameter :: lin_all     = 0        !! *linear interpolation between grids*
integer(kind=I4), parameter :: spl_all     = 2        !! *spline interpolation between grids*
integer(kind=I4), parameter :: hermite     = 4        !! *Hermite interpolation between grids*

integer(kind=I4), parameter :: method_asfc = hermite  !! *default*

integer(kind=I4) :: unit_out_lin  !! *output unit for the linear case*
integer(kind=I4) :: unit_out_spl  !! *output unit for the spline case*
integer(kind=I4) :: unit_out_her  !! *output unit for the Hermite case*

logical(kind=I4), parameter :: live             = .false.   !! *default for file outputs*
logical(kind=I4), parameter :: out_lin          = live      !! *default for file outputs, linear case*
logical(kind=I4), parameter :: out_spl          = live      !! *default for file outputs, spline case*
logical(kind=I4), parameter :: out_her          = live      !! *default for file outputs, Hermite case*
logical(kind=I4), parameter :: out_ter          = live      !! *default for terminal output*

public :: calcul_asfc_hermite, indice_fractal

contains

   subroutine calcul_aire(tab_in, long, larg, hx, hy, aire)
   !================================================================================================
   !! Return the area of a surface
   implicit none
   integer(kind=I4), intent(in )                            :: long     !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg     !! *2D array width*
   real   (kind=R8), intent(in )                            :: hx       !! *increment along x*
   real   (kind=R8), intent(in )                            :: hy       !! *increment along y*
   real   (kind=R8), intent(out)                            :: aire     !! *computed area*
   real   (kind=R8), intent(in ), dimension(1:long,1:larg)  :: tab_in   !! *surface array*

      integer(kind=I4) :: i, j
      real   (kind=R8) :: z1, z2, z3, z4, si

      si = 1!SCALE_IMG%si

      ! Raisonnement sur chaque carré du domaine
      aire = 0.
      do j = 1, larg -1
      do i = 1, long -1

         z1 = tab_in(i   , j   )*si
         z2 = tab_in(i   , j +1)*si
         z3 = tab_in(i +1, j +1)*si
         z4 = tab_in(i +1, j   )*si

         aire = aire +0.5_R8*( sqrt( UN +((z1-z2)/hx)**2 +((z1-z4)/hy)**2 ) + &
                               sqrt( UN +((z3-z2)/hy)**2 +((z3-z4)/hx)**2 ) )

      enddo
      enddo
      aire = aire/( (long -1)*(larg -1) )

   return
   endsubroutine calcul_aire


   subroutine calcul_asfc_lin_all(tab_in, scal, asfc_res)
   !================================================================================================
   !! Return the *asfc* of a surface. The different grids are obtained by linear interpolation
   implicit none
   type(SCALE_SURF), intent(in )                                      :: scal     !! *surface characteristics*
   real   (kind=R8), intent(in ), dimension(1:scal%xres, 1:scal%yres) :: tab_in   !! *input surface*
   real   (kind=R8), intent(out), dimension(1:2)                      :: asfc_res !! *result: asfc, adjustment factor*

      real(kind=R8), allocatable, dimension(:) :: x        ! x points in original grid
      real(kind=R8), allocatable, dimension(:) :: y        ! y points in original grid

      integer(kind=I4) :: long_new              ! number of points in x dimension for new grid
      integer(kind=I4) :: larg_new              ! number of points in y dimension for new grid

      real(kind=R8), allocatable, dimension(:)   :: x_new    ! new grid x points
      real(kind=R8), allocatable, dimension(:)   :: y_new    ! new grid y points
      real(kind=R8), allocatable, dimension(:,:) :: tab_ou   ! new grid function evaluations

      real   (kind=R8) :: rr
      integer(kind=I4) :: i, ii, j, jj, ip, long_tmp, larg_tmp

      integer(kind=I4) :: nb_pt

      real(kind=R8), dimension(1:nb_beta) :: beta
      real(kind=R8), dimension(1:npp)     :: vec_l

      real(kind=R8), dimension(1:npp) :: vec_x  ! points coordinates
      real(kind=R8), dimension(1:npp) :: vec_y  ! points coordinates

      real   (kind=R8) :: asfc1, asfc2, aire_lin, xm, xp, ym, yp, hx, hy, hhx, hhy, h1, h2, h3, h4, width, height
      logical(kind=I4) :: new_it
      integer(kind=I4) :: long, larg

      long = scal%xres
      larg = scal%yres

      if (out_lin) open(unit=unit_out_lin, file="out/asfc_lin_lin_all.txt")

      width  = scal%lx
      height = scal%ly

      hx =  width/(long -1)
      hy = height/(larg -1)

      ! définition d'abscisses pour l'interpolation par splines
      allocate( x(1:long), y(1:larg) )
      do i = 1, long
         x(i) = hx*real(i-1, kind=R8)
      enddo
      do j = 1, larg
         y(j) = hy*real(j-1, kind=R8)
      enddo

      rr =          (real(long0, kind=R8)/long)**(UN/npp)    ! facteur de réduction pour aller du maillage initial au maillage minimal avec npp points
      rr = max( rr, (real(larg0, kind=R8)/larg)**(UN/npp) )  ! facteur de réduction pour aller du maillage initial au maillage minimal avec npp points

      allocate( x_new(1:long),              &
                y_new(1:larg),              &  ! nouvelles abscisses
               tab_ou(1:long, 1:larg) )

      x_new(1:long)          = x(1:long)
      y_new(1:larg)          = y(1:larg)
      tab_ou(1:long, 1:larg) = tab_in(1:long, 1:larg)

      ! pour chaque réduction de maillage, calcul du maillage résultant et de l'aire relative associée
      !.............................................................
      long_new = long
      larg_new = larg
      hhx      = hx
      hhy      = hy
      new_it   = .true.
      nb_pt    = 1
      do ip = 1, npp +1

         if (new_it) then
            long_tmp = long_new
            larg_tmp = larg_new

            ! calcul de l'aire relative
            call calcul_aire(tab_ou(1:long_tmp, 1:larg_tmp), long_tmp, larg_tmp, hhx, hhy, aire_lin)

            if (nb_pt>1) then
               vec_x(nb_pt -1) = log( (hhx*1e6)*(hhy*1e6)/2 )
               vec_l(nb_pt -1) = log(aire_lin)
            endif

            if (out_lin .and. nb_pt>1) write(unit_out_lin, *) vec_x(nb_pt -1), vec_l(nb_pt -1)

            if (out_ter) write(*, *) aire_lin, long_tmp, larg_tmp
         endif

         if (ip == npp +1) then
            deallocate( x_new, y_new, x, y, tab_ou )
            exit
         endif

         long_new = nint(long*(rr**ip))   ! nb points en suite géométrique
         larg_new = nint(larg*(rr**ip))

         new_it = .true.
         if (long_new==long_tmp .or. larg_new==larg_tmp) then
            new_it = .false.
            cycle ! à découper trop fin, on peut tomber sur les mêmes entiers
         endif

         hhx =  width/(long_new -1)
         hhy = height/(larg_new -1)

         deallocate( x_new, y_new, tab_ou )
           allocate( x_new(1:long_new),              &
                     y_new(1:larg_new),              &  ! nouvelles abscisses
                    tab_ou(1:long_new, 1:larg_new) )

         nb_pt = nb_pt +1

         do i = 1, long_new
            x_new(i) = hhx*real(i-1, kind=R8)
         enddo

         do j = 1, larg_new
            y_new(j) = hhy*real(j-1, kind=R8)
         enddo

         do j = 1, larg_new
            jj = locate(n=larg, xx=y(1:larg), x=y_new(j))

            ym = y_new(j) -y(jj) ; yp = y_new(j) -y(jj +1)
            ym = ym/hy ; yp = yp/hy

            do i = 1, long_new
               ii = locate(n=long, xx=x(1:long), x=x_new(i))

               xm = x_new(i) -x(ii) ; xp = x_new(i) -x(ii +1)
               xm = xm/hx ; xp = xp/hx

               h1 = tab_in(ii   , jj   )
               h2 = tab_in(ii +1, jj   )
               h3 = tab_in(ii +1, jj +1)
               h4 = tab_in(ii   , jj +1)

               tab_ou(i, j) = h1*xp*yp -h2*xm*yp +h3*xm*ym -h4*xp*ym
            enddo

         enddo

      enddo
      if (out_lin) close(unit_out_lin)
      !call system('python pyt/filetoplot.py &')

      !.............................................................

      nb_pt = nb_pt -1

      vec_y(1:nb_pt) = vec_l(1:nb_pt)
      call interpolate_asfc( bt = beta(1:nb_beta),    &
                           n_bt = nb_beta,            &
                           n_pt = nb_pt,              &
                            asf = asfc1,              &
                             r2 = asfc2)

      asfc_res = [-1000*asfc1, asfc2]

   return
   contains

      subroutine interpolate_asfc(bt, n_bt, n_pt, asf, r2)
      !================================================================================================
      !< @note Function that fits the data points for the asfc determination.<br/>
      !        \(f\_boltz(x_i)=\beta_2 + \beta_3 \tanh \left( \dfrac{x_i -\beta_1}{\beta_4} \right) \)
      !
      !  @endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      integer(kind=I4), intent(in )                    :: n_bt !! *number of parameters*
      integer(kind=I4), intent(in )                    :: n_pt !! *data vector length*
      real   (kind=R8), intent(out), dimension(1:n_bt) :: bt   !! *vector \(\beta\) of parameters*
      real   (kind=R8), intent(out)                    :: asf  !! *Asfc number*
      real   (kind=R8), intent(out)                    :: r2   !! *correlation number to assess validity of the Asfc calculus*

         real   (kind=R8), dimension(1:n_pt, 1:n_bt) :: jacob
         real   (kind=R8), dimension(1:n_pt)         :: v_x, v_y, res_y, pentes
         integer(kind=I4)                            :: i0, i, j, info
         real   (kind=R8)                            :: delta1, delta2, y_mean

         v_x(1:n_pt) = vec_x(1:n_pt)

         ! smoothing
         v_y(1) = vec_y(1)
         do i = 1 + 1, n_pt - 1

            v_y(i) = 0.25_R8 * ( vec_y(i - 1) + 2 * vec_y(i) + vec_y(i + 1) )

         enddo
         v_y(n_pt) = vec_y(n_pt)

         call init_beta_boltz(   bt = bt(1:n_bt),     & !
                               n_bt = n_bt,           & !
                                v_x = v_x(1:n_pt),    & !
                                v_y = v_y(1:n_pt),    & !
                               n_pt = n_pt)             !

         res_y(1:n_pt)         = 0._R8
         jacob(1:n_pt, 1:n_bt) = 0._R8

         call lmder1(   fcn = lmder1_f,               & !
                          m = n_pt,                   & !
                          n = n_bt,                   & !
                          x = bt(1:n_bt),             & !
                       fvec = res_y(1:n_pt),          & !
                       fjac = jacob(1:n_pt, 1:n_bt),  & !
                     ldfjac = n_pt,                   & !
                        tol = 1.0e-8_R8,              & !
                       info = info)                     !

         delta1 = 0._R8
         delta2 = 0._R8
         y_mean = sum( vec_y(1:n_pt) )/n_pt
         do i = 1, n_pt
            delta1 = delta1 +( vec_y(i) -f_boltz(xi = vec_x(i),   & !
                                               beta = bt(1:n_bt), & !
                                             n_beta = n_bt)       & !
                             )**2
            delta2 = delta2 +( vec_y(i) -y_mean )**2
         enddo
         r2 = 1._R8 -delta1/delta2

   !~       if (r2<0) then
   !~       do i=1,n_pt
   !~       write(99,*) v_x(i), v_y(i)
   !~       enddo
   !~       stop 'error'
   !~       endif

         i0 = locate(n = n_pt,        & !
                    xx = v_x(1:n_pt), & !
                     x = bt(1))         !

         j = 0
         do i = i0 -5, i0 +5
            if (i<1 .or. i>n_pt) cycle
            j = j +1
            pentes(j) = +df_boltz(  xi = v_x(i),     & !
                                  beta = bt(1:n_bt), & !
                                n_beta = n_bt,       & !
                                  ivar = 0)            !
         enddo

         asf = sum( pentes(1:j)/j )

      return
      endsubroutine interpolate_asfc

      subroutine lmder1_f(m, n, x, fvec, fjac, ldfjac, iflag)
      !================================================================================================
      !< @note Function called by [[lmder1]] as part of **minpack**, modified by
      !        [Burkardt](https://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html) <br/>
      !        According *iflag* value it calculates the function [[f_boltz]] at the data points
      !        or the jacobian.
      !
      !  @endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      integer(kind=I4), intent(in   )                           :: m       !! *number of points*
      integer(kind=I4), intent(in   )                           :: n       !! *number of parameters*
      integer(kind=I4), intent(in   )                           :: ldfjac  !! *leading dimension of fjac, which must be at least n*
      integer(kind=I4), intent(in   )                           :: iflag   !! *which calculus to perform*
      real   (kind=R8), intent(out  ), dimension(1:m)           :: fvec    !! *vector of f_boltz(xi)*
      real   (kind=R8), intent(out  ), dimension(1:ldfjac, 1:n) :: fjac    !! *jacobian*
      real   (kind=R8), intent(inout), dimension(1:n)           :: x       !! *parameter values*

         integer(kind=I4) :: i, k

         select case (iflag)
            case(0)
               continue

            case(1)
               do i = 1, m
                  fvec(i) = f_boltz(     xi = vec_x(i),  &  !
                                       beta = x(1:n),    &  !
                                     n_beta = n)         &  !
                            -vec_y(i)
               enddo

            case(2)
               do i = 1, m
               do k = 1, n
                  fjac(i, k) = df_boltz(xi = vec_x(i), & !
                                      beta = x(1:n),   & !
                                      ivar = k,        & !
                                    n_beta = n)          !
               enddo
               enddo

            case default
               write ( *, '(a)' ) ' '
               write ( *, '(a)' ) 'LMDER1_F - Fatal error!'
               write ( *, '(a,i6)' ) '  Called with unexpected value of IFLAG = ', iflag
               stop

         endselect
      return
      endsubroutine lmder1_f

   endsubroutine calcul_asfc_lin_all


   subroutine calcul_asfc_spl_all(tab_in, scal, asfc_res)
   !================================================================================================
  !! Return the *asfc* of a surface. The different grids are obtained by spline of degree 3
   implicit none
   type(SCALE_SURF), intent(in )                                      :: scal     !! *surface characteristics*
   real   (kind=R8), intent(in ), dimension(1:scal%xres, 1:scal%yres) :: tab_in   !! *input surface*
   real   (kind=R8), intent(out), dimension(1:2)                      :: asfc_res !! *result: asfc, adjustment factor*

      integer(kind=I4), parameter :: idx = 0    ! [[db2val]] input
      integer(kind=I4), parameter :: idy = 0    ! [[db2val]] input

      real(kind=R8), allocatable, dimension(:) :: x        ! x points in original grid
      real(kind=R8), allocatable, dimension(:) :: y        ! y points in original grid

      integer(kind=I4) :: long_new              ! number of points in x dimension for new grid
      integer(kind=I4) :: larg_new              ! number of points in y dimension for new grid

      real(kind=R8), allocatable, dimension(:)   :: x_new    ! new grid x points
      real(kind=R8), allocatable, dimension(:)   :: y_new    ! new grid y points
      real(kind=R8), allocatable, dimension(:,:) :: tab_ou   ! new grid function evaluations

      real(kind=R8), allocatable, dimension(:,:) :: coeff, coeff_tmp
      real(kind=R8), allocatable, dimension(:)   :: tx, tx_tmp       ! x knots
      real(kind=R8), allocatable, dimension(:)   :: ty, ty_tmp       ! y knots

      real(kind=R8)    :: val, rr
      integer(kind=I4) :: i, j, k, ip, long_tmp, larg_tmp
      integer(kind=I4) :: iflag
      integer(kind=I4) :: inbvx, inbvy, iloy
      integer(kind=I4) :: nb_pt

      integer(kind=I4), dimension(1:8) :: d1, d2

      real(kind=R8), dimension(1:nb_beta) :: beta
      real(kind=R8), dimension(1:npp)     :: vec_l, vec_s

      real(kind=R8), dimension(1:npp) :: vec_x  ! points coordinates
      real(kind=R8), dimension(1:npp) :: vec_y  ! points coordinates

      real(kind=R8), dimension(1:8)       :: gx, gy, vf

      real(kind=R8) :: asfc1, asfc2, aire, aire_lin, aire_tmp, x1, x2, y1, y2, val_x, val_y, width, height, hx, hy, hhx, hhy

      logical(kind=I4) :: new_it

      integer(kind=I4) :: long, larg

      long = scal%xres
      larg = scal%yres

      if (out_lin) open(unit=unit_out_lin, file="out/asfc_lin_spl_all.txt")
      if (out_spl) open(unit=unit_out_spl, file="out/asfc_spl_spl_all.txt")

      width  = scal%lx
      height = scal%ly

      hx =  width/(long -1)
      hy = height/(larg -1)

      ! définition d'abscisses pour l'interpolation par splines
      allocate( x(1:long), y(1:larg) )
      do i = 1, long
         x(i) = hx*real(i-1, kind=R8)
      enddo
      do j = 1, larg
         y(j) = hy*real(j-1, kind=R8)
      enddo

      rr =          (real(long0, kind=R8)/long)**(UN/npp)    ! facteur de réduction pour aller du maillage initial au maillage minimal avec npp points
      rr = max( rr, (real(larg0, kind=R8)/larg)**(UN/npp) )  ! facteur de réduction pour aller du maillage initial au maillage minimal avec npp points

      allocate( x_new(1:long),              &
                y_new(1:larg),              &  ! nouvelles abscisses
               tab_ou(1:long, 1:larg) )

      x_new(1:long)          = x(1:long)
      y_new(1:larg)          = y(1:larg)
      tab_ou(1:long, 1:larg) = tab_in(1:long, 1:larg)

      ! pour chaque réduction de maillage, calcul du maillage résultant et de l'aire relative associée
      !.............................................................
      long_new = long
      larg_new = larg
      hhx      = hx
      hhy      = hy
      nb_pt    = 1
      new_it   = .true.
      do ip = 1, npp +1

         if (new_it) then

            long_tmp = long_new
            larg_tmp = larg_new

            ! calcul de l'aire relative
            call calcul_aire(tab_ou(1:long_tmp, 1:larg_tmp), long_tmp, larg_tmp, hhx, hhy, aire_lin)!, rx=real(long-1, kind=R8)/(long_tmp-1), &
                                                                                          ! ry=real(larg-1, kind=R8)/(larg_tmp-1))

            allocate( coeff_tmp(1:long_tmp, 1:larg_tmp) )

            allocate( tx_tmp(1:(long_tmp +kx)),  &
                      ty_tmp(1:(larg_tmp +ky)) )

            iflag = 0
            call db2ink(   x = x_new(1:long_tmp),                 &
                          nx = long_tmp,                          &
                           y = y_new(1:larg_tmp),                 &
                          ny = larg_tmp,                          &
                         fcn = tab_ou(1:long_tmp, 1:larg_tmp),    &
                          kx = kx,                                &
                          ky = ky,                                &
                          tx = tx_tmp(1:(long_tmp +kx)),          &
                          ty = ty_tmp(1:(larg_tmp +ky)),          &
                       bcoef = coeff_tmp(1:long_tmp, 1:larg_tmp), &
                       iflag = iflag)

            if (iflag/=1) then
               write(*,*) iflag, long_tmp, larg_tmp, kx, ky
               stop 'error calling db2ink'
            endif

            if (ip==1) then
               allocate( coeff(1:long, 1:larg) )
               allocate( tx(1:(long +kx)),  &
                         ty(1:(larg +ky)) )
               coeff(1:long, 1:larg) = coeff_tmp(1:long, 1:larg)
                  tx(1:(long +kx))   =    tx_tmp(1:(long +kx))
                  ty(1:(larg +ky))   =    ty_tmp(1:(larg +ky))
            endif

            inbvx = 1
            inbvy = 1
            iloy  = 1

            aire_tmp = 0._R8
!~             !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NB_PROCS) IF(MULTI_PROCS2)
!~             !$OMP DO SCHEDULE (STATIC,(larg_tmp-1)/NB_PROCS) PRIVATE(i, k, iflag, val, y1, y2, x1, x2, d1, d2, gx, gy, vf) FIRSTPRIVATE(inbvx, inbvy, iloy) REDUCTION(+:aire_tmp)
               do j = 1, larg_tmp -1
                  y1 = y_new(j) +(hhy/2)*(UN -UN/sqrt(3._R8))
                  y2 = y_new(j) +(hhy/2)*(UN +UN/sqrt(3._R8))

                  do i = 1, long_tmp -1
                     x1 = x_new(i) +(hhx/2)*(UN -UN/sqrt(3._R8))
                     x2 = x_new(i) +(hhx/2)*(UN +UN/sqrt(3._R8))

                     d1(1:8) = [ 1,  0,  1,  0,  1,  0,  1,  0]
                     d2(1:8) = [ 0,  1,  0,  1,  0,  1,  0,  1]
                     gx(1:8) = [x1, x1, x1, x1, x2, x2, x2, x2]
                     gy(1:8) = [y1, y1, y2, y2, y1, y1, y2, y2]

                     do k = 1, 8
                        call db2val(xval = gx(k),                             &
                                    yval = gy(k),                             &
                                     idx = d1(k),                             &
                                     idy = d2(k),                             &
                                      tx = tx_tmp(1:(long_tmp +kx)),          &
                                      ty = ty_tmp(1:(larg_tmp +ky)),          &
                                      nx = long_tmp,                          &
                                      ny = larg_tmp,                          &
                                      kx = kx,                                &
                                      ky = ky,                                &
                                   bcoef = coeff_tmp(1:long_tmp, 1:larg_tmp), &
                                       f = vf(k),                             &
                                   iflag = iflag,                             &
                                   inbvx = inbvx,                             &
                                   inbvy = inbvy,                             &
                                    iloy = iloy)
                        if (iflag/=0) then
                           write(*,*) iflag
                           stop 'error calling db2val'
                        endif
                     enddo
                     do k = 1, 4
                        val_x = vf(2*k -1)
                        val_y = vf(2*k   )
                        aire_tmp = aire_tmp +sqrt( UN +val_x**2 +val_y**2 )
                     enddo

                  enddo
               enddo
!~             !$OMP END DO
!~             !$OMP END PARALLEL

            aire = aire_tmp*( hhx/2 )*( hhy/2 )
            aire = aire/(width*height)

            vec_x(nb_pt) = log( (hhx*1e6)*(hhy*1e6)/2 )
            vec_l(nb_pt) = log(aire_lin)
            vec_s(nb_pt) = log(aire    )

            if (out_lin .and. nb_pt>1) write(unit_out_lin, *) vec_x(nb_pt -1), vec_l(nb_pt -1)
            if (out_spl .and. nb_pt>1) write(unit_out_spl, *) vec_x(nb_pt -1), vec_s(nb_pt -1)

            if (out_ter) write(*, *) hhx, hhy, aire, ip

         endif

         if (ip == npp +1) then
            deallocate( x_new, y_new, x, y, tab_ou, coeff, tx, ty, coeff_tmp, tx_tmp, ty_tmp )
            exit
         endif

         long_new = nint(long*(rr**ip))   ! nb points en suite géométrique
         larg_new = nint(larg*(rr**ip))

         new_it = .true.
         if (long_new==long_tmp .or. larg_new==larg_tmp) then
            new_it = .false.
            cycle ! à découper trop fin, on peut tomber sur les mêmes entiers
         endif

         hhx =  width/(long_new -1)
         hhy = height/(larg_new -1)

         nb_pt = nb_pt +1

         deallocate( x_new, y_new, tab_ou )
           allocate( x_new(1:long_new),              &
                     y_new(1:larg_new),              &  ! nouvelles abscisses
                    tab_ou(1:long_new, 1:larg_new) )

         do i = 1, long_new
            x_new(i) = hhx*real(i-1, kind=R8)
         enddo

         do j = 1, larg_new
            y_new(j) = hhy*real(j-1, kind=R8)
         enddo

         inbvx = 1
         inbvy = 1
         iloy  = 1

         ! calcul des hauteurs de la surface interpolée
!~          !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NB_PROCS) IF(MULTI_PROCS2)
!~          !$OMP DO SCHEDULE (STATIC,larg_new/NB_PROCS) PRIVATE(i, iflag, val) FIRSTPRIVATE(inbvx, inbvy, iloy)
            do j = 1, larg_new
            do i = 1, long_new
               call db2val(xval = x_new(i),              &
                           yval = y_new(j),              &
                            idx = idx,                   &
                            idy = idy,                   &
                             tx = tx(1:(long +kx)),      &
                             ty = ty(1:(larg +ky)),      &
                             nx = long,                  &
                             ny = larg,                  &
                             kx = kx,                    &
                             ky = ky,                    &
                          bcoef = coeff(1:long, 1:larg), &
                              f = val,                   &
                          iflag = iflag,                 &
                          inbvx = inbvx,                 &
                          inbvy = inbvy,                 &
                           iloy = iloy)
               if (iflag/=0) error stop 'error calling db2val'
               iflag = 0
               tab_ou(i, j) = val
            enddo
            enddo
!~          !$OMP END DO
!~          !$OMP END PARALLEL

         deallocate( coeff_tmp, tx_tmp, ty_tmp )

      enddo
      if (out_lin) close(unit_out_lin)
      if (out_spl) close(unit_out_spl)
      !.............................................................
      nb_pt = nb_pt -1

      vec_y(1:nb_pt) = vec_s(1:nb_pt)
      call interpolate_asfc( bt = beta(1:nb_beta),    &
                           n_bt = nb_beta,            &
                           n_pt = nb_pt,              &
                            asf = asfc1,              &
                             r2 = asfc2)

      asfc_res = [-1000*asfc1, asfc2]

   return
   contains

      subroutine interpolate_asfc(bt, n_bt, n_pt, asf, r2)
      !================================================================================================
      !< @note Function that fits the data points for the asfc determination.<br/>
      !        \(f\_boltz(x_i)=\beta_2 + \beta_3 \tanh \left( \dfrac{x_i -\beta_1}{\beta_4} \right) \)
      !
      !  @endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      integer(kind=I4), intent(in )                    :: n_bt !! *number of parameters*
      integer(kind=I4), intent(in )                    :: n_pt !! *data vector length*
      real   (kind=R8), intent(out), dimension(1:n_bt) :: bt   !! *vector \(\beta\) of parameters*
      real   (kind=R8), intent(out)                    :: asf  !! *Asfc number*
      real   (kind=R8), intent(out)                    :: r2   !! *correlation number to assess validity of the Asfc calculus*

         real   (kind=R8), dimension(1:n_pt, 1:n_bt) :: jacob
         real   (kind=R8), dimension(1:n_pt)         :: v_x, v_y, res_y, pentes
         integer(kind=I4)                            :: i0, i, j, info
         real   (kind=R8)                            :: delta1, delta2, y_mean

         v_x(1:n_pt) = vec_x(1:n_pt)

         ! smoothing
         v_y(1) = vec_y(1)
         do i = 1 + 1, n_pt - 1

            v_y(i) = 0.25_R8 * ( vec_y(i - 1) + 2 * vec_y(i) + vec_y(i + 1) )

         enddo
         v_y(n_pt) = vec_y(n_pt)

         call init_beta_boltz(   bt = bt(1:n_bt),     & !
                               n_bt = n_bt,           & !
                                v_x = v_x(1:n_pt),    & !
                                v_y = v_y(1:n_pt),    & !
                               n_pt = n_pt)             !

         res_y(1:n_pt)         = 0._R8
         jacob(1:n_pt, 1:n_bt) = 0._R8

         call lmder1(   fcn = lmder1_f,               & !
                          m = n_pt,                   & !
                          n = n_bt,                   & !
                          x = bt(1:n_bt),             & !
                       fvec = res_y(1:n_pt),          & !
                       fjac = jacob(1:n_pt, 1:n_bt),  & !
                     ldfjac = n_pt,                   & !
                        tol = 1.0e-8_R8,              & !
                       info = info)                     !

         delta1 = 0._R8
         delta2 = 0._R8
         y_mean = sum( vec_y(1:n_pt) )/n_pt
         do i = 1, n_pt
            delta1 = delta1 +( vec_y(i) -f_boltz(xi = vec_x(i),   & !
                                               beta = bt(1:n_bt), & !
                                             n_beta = n_bt)       & !
                             )**2
            delta2 = delta2 +( vec_y(i) -y_mean )**2
         enddo
         r2 = 1._R8 -delta1/delta2

   !~       if (r2<0) then
   !~       do i=1,n_pt
   !~       write(99,*) v_x(i), v_y(i)
   !~       enddo
   !~       stop 'error'
   !~       endif

         i0 = locate(n = n_pt,        & !
                    xx = v_x(1:n_pt), & !
                     x = bt(1))         !

         j = 0
         do i = i0 -5, i0 +5
            if (i<1 .or. i>n_pt) cycle
            j = j +1
            pentes(j) = +df_boltz(  xi = v_x(i),     & !
                                  beta = bt(1:n_bt), & !
                                n_beta = n_bt,       & !
                                  ivar = 0)            !
         enddo

         asf = sum( pentes(1:j)/j )

      return
      endsubroutine interpolate_asfc

      subroutine lmder1_f(m, n, x, fvec, fjac, ldfjac, iflag)
      !================================================================================================
      !< @note Function called by [[lmder1]] as part of **minpack**, modified by
      !        [Burkardt](https://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html) <br/>
      !        According *iflag* value it calculates the function [[f_boltz]] at the data points
      !        or the jacobian.
      !
      !  @endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      integer(kind=I4), intent(in   )                           :: m       !! *number of points*
      integer(kind=I4), intent(in   )                           :: n       !! *number of parameters*
      integer(kind=I4), intent(in   )                           :: ldfjac  !! *leading dimension of fjac, which must be at least n*
      integer(kind=I4), intent(in   )                           :: iflag   !! *which calculus to perform*
      real   (kind=R8), intent(out  ), dimension(1:m)           :: fvec    !! *vector of f_boltz(xi)*
      real   (kind=R8), intent(out  ), dimension(1:ldfjac, 1:n) :: fjac    !! *jacobian*
      real   (kind=R8), intent(inout), dimension(1:n)           :: x       !! *parameter values*

         integer(kind=I4) :: i, k

         select case (iflag)
            case(0)
               continue

            case(1)
               do i = 1, m
                  fvec(i) = f_boltz(     xi = vec_x(i),  &  !
                                       beta = x(1:n),    &  !
                                     n_beta = n)         &  !
                            -vec_y(i)
               enddo

            case(2)
               do i = 1, m
               do k = 1, n
                  fjac(i, k) = df_boltz(xi = vec_x(i), & !
                                      beta = x(1:n),   & !
                                      ivar = k,        & !
                                    n_beta = n)          !
               enddo
               enddo

            case default
               write ( *, '(a)' ) ' '
               write ( *, '(a)' ) 'LMDER1_F - Fatal error!'
               write ( *, '(a,i6)' ) '  Called with unexpected value of IFLAG = ', iflag
               stop

         endselect
      return
      endsubroutine lmder1_f

   endsubroutine calcul_asfc_spl_all


   subroutine calcul_asfc_hermite(tab_in, scal, asfc_res, omp)
   !================================================================================================
   !! Return the *asfc* of a surface. The different grids are obtained by Hermite interpolation
   implicit none
   type(SCALE_SURF), intent(in )                                      :: scal     !! *surface characteristics*
   real   (kind=R8), intent(in ), dimension(1:scal%xres, 1:scal%yres) :: tab_in   !! *input surface*
   real   (kind=R8), intent(out), dimension(1:2)                      :: asfc_res !! *result: asfc, adjustment factor*
   logical(kind=I4), intent(in )                                      :: omp      !! *with openmp ?*

      real(kind=R8), allocatable, dimension(:) :: x            ! x points in original grid
      real(kind=R8), allocatable, dimension(:) :: y            ! y points in original grid

      integer(kind=I4) :: long_new                             ! number of points in x dimension for new grid
      integer(kind=I4) :: larg_new                             ! number of points in y dimension for new grid

      real(kind=R8), allocatable, dimension(:)   :: x_new      ! new grid x points
      real(kind=R8), allocatable, dimension(:)   :: y_new      ! new grid y points
      real(kind=R8), allocatable, dimension(:,:) :: tab_ou     ! new grid function evaluations


      real(kind=R8), allocatable, dimension(:,:) :: tab_in_dx  ! new grid function evaluations
      real(kind=R8), allocatable, dimension(:,:) :: tab_in_dy  ! new grid function evaluations
      real(kind=R8), allocatable, dimension(:,:) :: tab_in_xy  ! new grid function evaluations

      real(kind=R8), allocatable, dimension(:) :: gx
      real(kind=R8), allocatable, dimension(:) :: gy
      real(kind=R8), allocatable, dimension(:) :: gw
      real(kind=R8), allocatable, dimension(:) :: tab_dnq

      real(kind=R8)    :: rr
      integer(kind=I4) :: i, ii, j, jj, k, long_tmp, larg_tmp

      real(kind=R8), dimension(1:nb_beta) :: beta
      real(kind=R8), dimension(1:npp)     :: vec_s

      real(kind=R8), dimension(1:npp) :: vec_x  ! points coordinates
      real(kind=R8), dimension(1:npp) :: vec_y  ! points coordinates

      real(kind=R8) :: asfc1, asfc2, aire, hx, hy, hhx, hhy, width, height


      real   (kind=R8) :: xi, yi, eps_x
      integer(kind=I4) :: it, ng, nbpt

      integer(kind=I4) :: long, larg

      integer(kind=I4) :: nb_th

      long = scal%xres
      larg = scal%yres

      if (out_her) call get_unit(unit_out_her)
      if (out_her) open(unit=unit_out_her, file="out/asfc_her_her_all.txt")

      width  = scal%lx
      height = scal%ly

      hx =  width/(long -1)
      hy = height/(larg -1)

      eps_x = min(hx/10**3, hy/10**3)

      ! définition d'abscisses pour l'interpolation par splines
      allocate( x(1:long), y(1:larg) )
      do i = 1, long
         x(i) = hx*real(i-1, kind=R8)
      enddo
      do j = 1, larg
         y(j) = hy*real(j-1, kind=R8)
      enddo


      allocate( x_new(1:long),               &  !
                y_new(1:larg),               &  !
               tab_ou(1:long, 1:larg),       &  !
               tab_in_dx(1:long, 1:larg),    &  !
               tab_in_dy(1:long, 1:larg),    &  !
               tab_in_xy(1:long, 1:larg) )      !

      rr =          (real(long0, kind=R8)/long)**(UN/npp)    ! facteur de réduction pour aller du maillage initial au maillage minimal avec npp points
      rr = max( rr, (real(larg0, kind=R8)/larg)**(UN/npp) )  ! facteur de réduction pour aller du maillage initial au maillage minimal avec npp points

      x_new(1:long)          = x(1:long)
      y_new(1:larg)          = y(1:larg)
      tab_ou(1:long, 1:larg) = tab_in(1:long, 1:larg)


      call init_aire_hermite(gx = gx, gy = gy, gw = gw, tab_dnq = tab_dnq, ng = ng)

      ! pour chaque réduction de maillage, calcul du maillage résultant et de l'aire relative associée
      !.............................................................
      long_new = long
      larg_new = larg
      hhx      = hx
      hhy      = hy


      call calcul_aire_hermite(   tab_in = tab_ou(1:long_new, 1:larg_new), &  !
                                    long = long_new,                       &  !
                                    larg = larg_new,                       &  !
                                      gw = gw(1:ng),                       &  !
                                 tab_dnq = tab_dnq(1:ng),                  &  !
                                      ng = ng,                             &  !
                                      hx = hhx,                            &  !
                                      hy = hhy,                            &  !
                                   width = width,                          &  !
                                  height = height,                         &  !
                                    aire = aire)                              !


      vec_s(1) = log( aire )
      vec_x(1) = log( (hhx*1e6)*(hhy*1e6)/2 )

      call calcul_tabd_hermite( tab_in = tab_in   (1:long, 1:larg),  &  !
                                tab_dx = tab_in_dx(1:long, 1:larg),  &  !
                                tab_dy = tab_in_dy(1:long, 1:larg),  &  !
                                tab_xy = tab_in_xy(1:long, 1:larg),  &  !
                                  long = long,                       &  !
                                  larg = larg,                       &  !
                                    hx = hx,                         &  !
                                    hy = hy )                           !

      nb_th = 1
      if (omp) then
         nb_th = omp_get_num_procs()
      endif


      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nb_th) IF (omp)
      !$OMP DO ORDERED SCHEDULE (STATIC,1) PRIVATE(it, long_tmp, larg_tmp, long_new, larg_new, hhx, hhy, i, x_new, j, y_new, yi, jj, xi, ii, tab_ou, aire)


      do it = 1, npp - 1

         long_tmp = nint(long*(rr**it))   ! nb points en suite géométrique
         larg_tmp = nint(larg*(rr**it))

         if ( long_new == long_tmp .or. larg_new == larg_tmp) then
            vec_s(it + 1) = 0
            cycle ! à découper trop fin, on peut tomber sur les mêmes entiers
         endif

         long_new = long_tmp
         larg_new = larg_tmp

         hhx =  width/(long_new -1)
         hhy = height/(larg_new -1)

         vec_x(it + 1) = log( (hhx*1e6)*(hhy*1e6)/2 )


         deallocate( x_new, y_new, tab_ou )

         allocate( x_new(1:long_new),                 &  !
                   y_new(1:larg_new),                 &  ! nouvelles abscisses
                  tab_ou(1:long_new, 1:larg_new) )       !

         do i = 1, long_new
            x_new(i) = hhx*real(i-1, kind=R8)
         enddo

         do j = 1, larg_new
            y_new(j) = hhy*real(j-1, kind=R8)
         enddo

         do j = 1, larg_new

            yi = y_new(j)
            jj = locate2(n = larg, xx = y(1:larg), x = yi, eps = eps_x)
            yi = (yi -(y(jj)+y(jj+1))/2)/(hy/2)

            do i = 1, long_new

               xi = x_new(i)
               ii = locate2(n = long, xx = x(1:long), x = xi, eps = eps_x)
               xi = (xi -(x(ii)+x(ii+1))/2)/(hx/2)

               tab_ou(i, j) =                                              &  !
                     nq_i(xi, yi, 1, 1)*tab_in(ii   , jj   ) +             &  !  u1
                     nq_i(xi, yi, 2, 1)*tab_in(ii +1, jj   ) +             &  !  u2
                     nq_i(xi, yi, 3, 1)*tab_in(ii +1, jj +1) +             &  !  u3
                     nq_i(xi, yi, 4, 1)*tab_in(ii   , jj +1) +             &  !  u4

                     nq_i(xi, yi, 1, 2)*tab_in_dx(ii   , jj   )*hx/2 +     &  ! du1/dx
                     nq_i(xi, yi, 2, 2)*tab_in_dx(ii +1, jj   )*hx/2 +     &  ! du2/dx
                     nq_i(xi, yi, 3, 2)*tab_in_dx(ii +1, jj +1)*hx/2 +     &  ! du3/dx
                     nq_i(xi, yi, 4, 2)*tab_in_dx(ii   , jj +1)*hx/2 +     &  ! du4/dx

                     nq_i(xi, yi, 1, 3)*tab_in_dy(ii   , jj   )*hy/2 +     &  ! du1/dy
                     nq_i(xi, yi, 2, 3)*tab_in_dy(ii +1, jj   )*hy/2 +     &  ! du2/dy
                     nq_i(xi, yi, 3, 3)*tab_in_dy(ii +1, jj +1)*hy/2 +     &  ! du3/dy
                     nq_i(xi, yi, 4, 3)*tab_in_dy(ii   , jj +1)*hy/2 +     &  ! du4/dy

                     nq_i(xi, yi, 1, 4)*tab_in_xy(ii   , jj   )*hx*hy/4 +  &  ! du1/dxdy
                     nq_i(xi, yi, 2, 4)*tab_in_xy(ii +1, jj   )*hx*hy/4 +  &  ! du2/dxdy
                     nq_i(xi, yi, 3, 4)*tab_in_xy(ii +1, jj +1)*hx*hy/4 +  &  ! du3/dxdy
                     nq_i(xi, yi, 4, 4)*tab_in_xy(ii   , jj +1)*hx*hy/4       ! du4/dxdy

            enddo
         enddo

         call calcul_aire_hermite(   tab_in = tab_ou(1:long_new, 1:larg_new), &  !
                                       long = long_new,                       &  !
                                       larg = larg_new,                       &  !
                                         gw = gw(1:ng),                       &  !
                                    tab_dnq = tab_dnq(1:ng),                  &  !
                                         ng = ng,                             &  !
                                         hx = hhx,                            &  !
                                         hy = hhy,                            &  !
                                      width = width,                          &  !
                                     height = height,                         &  !
                                       aire = aire)                              !

         vec_s(it + 1) = log( aire )

!~          write(unit_out_her, *) vec_x(it+1), vec_s(it+1)

      enddo

      !$OMP END DO
      !$OMP END PARALLEL

      deallocate( x_new, y_new, x, y, tab_ou, tab_in_dx, tab_in_dy, tab_in_xy, tab_dnq, gx, gy, gw )

      if (out_lin) close(unit_out_lin)
      if (out_her) close(unit_out_her)
      !.............................................................

      k = 1
      do i = 1, npp
         if ( abs(vec_s(i)) > EPS_R8 ) then
            vec_y(k) = vec_s(i)
            vec_x(k) = vec_x(i)
            k = k + 1
         endif
      enddo
      nbpt = k - 1

      call interpolate_asfc( bt = beta(1:nb_beta),    &  !
                           n_bt = nb_beta,            &  !
                           n_pt = nbpt,               &  !
                            asf = asfc1,              &  !
                             r2 = asfc2 )                !

      asfc_res = [-1000*asfc1, asfc2]

   return

   contains

      subroutine interpolate_asfc(bt, n_bt, n_pt, asf, r2)
      !================================================================================================
      !< @note Function that fits the data points for the asfc determination.<br/>
      !        \(f\_boltz(x_i)=\beta_2 + \beta_3 \tanh \left( \dfrac{x_i -\beta_1}{\beta_4} \right) \)
      !
      !  @endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      integer(kind=I4), intent(in )                    :: n_bt !! *number of parameters*
      integer(kind=I4), intent(in )                    :: n_pt !! *data vector length*
      real   (kind=R8), intent(out), dimension(1:n_bt) :: bt   !! *vector \(\beta\) of parameters*
      real   (kind=R8), intent(out)                    :: asf  !! *Asfc number*
      real   (kind=R8), intent(out)                    :: r2   !! *correlation number to assess validity of the Asfc calculus*

         real   (kind=R8), dimension(1:n_pt, 1:n_bt) :: jacob
         real   (kind=R8), dimension(1:n_pt)         :: v_x, v_y, res_y, pentes
         integer(kind=I4)                            :: i0, i, j, info
         real   (kind=R8)                            :: delta1, delta2, y_mean

         v_x(1:n_pt) = vec_x(1:n_pt)

         ! smoothing
         v_y(1) = vec_y(1)
         do i = 1 + 1, n_pt - 1

            v_y(i) = 0.25_R8 * ( vec_y(i - 1) + 2 * vec_y(i) + vec_y(i + 1) )

         enddo
         v_y(n_pt) = vec_y(n_pt)

         call init_beta_boltz(   bt = bt(1:n_bt),     & !
                               n_bt = n_bt,           & !
                                v_x = v_x(1:n_pt),    & !
                                v_y = v_y(1:n_pt),    & !
                               n_pt = n_pt)             !

         res_y(1:n_pt)         = 0._R8
         jacob(1:n_pt, 1:n_bt) = 0._R8

         call lmder1(   fcn = lmder1_f,               & !
                          m = n_pt,                   & !
                          n = n_bt,                   & !
                          x = bt(1:n_bt),             & !
                       fvec = res_y(1:n_pt),          & !
                       fjac = jacob(1:n_pt, 1:n_bt),  & !
                     ldfjac = n_pt,                   & !
                        tol = 1.0e-8_R8,              & !
                       info = info)                     !

         delta1 = 0._R8
         delta2 = 0._R8
         y_mean = sum( vec_y(1:n_pt) )/n_pt
         do i = 1, n_pt
            delta1 = delta1 +( vec_y(i) -f_boltz(xi = vec_x(i),   & !
                                               beta = bt(1:n_bt), & !
                                             n_beta = n_bt)       & !
                             )**2
            delta2 = delta2 +( vec_y(i) -y_mean )**2
         enddo
         r2 = 1._R8 -delta1/delta2

   !~       if (r2<0) then
   !~       do i=1,n_pt
   !~       write(99,*) v_x(i), v_y(i)
   !~       enddo
   !~       stop 'error'
   !~       endif

         i0 = locate(n = n_pt,        & !
                    xx = v_x(1:n_pt), & !
                     x = bt(1))         !

         j = 0
         do i = i0 -5, i0 +5
            if (i<1 .or. i>n_pt) cycle
            j = j +1
            pentes(j) = +df_boltz(  xi = v_x(i),     & !
                                  beta = bt(1:n_bt), & !
                                n_beta = n_bt,       & !
                                  ivar = 0)            !
         enddo

         asf = sum( pentes(1:j)/j )

      return
      endsubroutine interpolate_asfc

      subroutine lmder1_f(m, n, x, fvec, fjac, ldfjac, iflag)
      !================================================================================================
      !< @note Function called by [[lmder1]] as part of **minpack**, modified by
      !        [Burkardt](https://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html) <br/>
      !        According *iflag* value it calculates the function [[f_boltz]] at the data points
      !        or the jacobian.
      !
      !  @endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      integer(kind=I4), intent(in   )                           :: m       !! *number of points*
      integer(kind=I4), intent(in   )                           :: n       !! *number of parameters*
      integer(kind=I4), intent(in   )                           :: ldfjac  !! *leading dimension of fjac, which must be at least n*
      integer(kind=I4), intent(in   )                           :: iflag   !! *which calculus to perform*
      real   (kind=R8), intent(out  ), dimension(1:m)           :: fvec    !! *vector of f_boltz(xi)*
      real   (kind=R8), intent(out  ), dimension(1:ldfjac, 1:n) :: fjac    !! *jacobian*
      real   (kind=R8), intent(inout), dimension(1:n)           :: x       !! *parameter values*

         integer(kind=I4) :: i, k

         select case (iflag)
            case(0)
               continue

            case(1)
               do i = 1, m
                  fvec(i) = f_boltz(     xi = vec_x(i),  &  !
                                       beta = x(1:n),    &  !
                                     n_beta = n)         &  !
                            -vec_y(i)
               enddo

            case(2)
               do i = 1, m
               do k = 1, n
                  fjac(i, k) = df_boltz(xi = vec_x(i), & !
                                      beta = x(1:n),   & !
                                      ivar = k,        & !
                                    n_beta = n)          !
               enddo
               enddo

            case default
               write ( *, '(a)' ) ' '
               write ( *, '(a)' ) 'LMDER1_F - Fatal error!'
               write ( *, '(a,i6)' ) '  Called with unexpected value of IFLAG = ', iflag
               stop

         endselect
      return
      endsubroutine lmder1_f

      subroutine init_aire_hermite(gx, gy, gw, tab_dnq, ng)
      !================================================================================================
      !< @note
      !
      !  @endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      real   (kind=R8), intent(out), allocatable, dimension(:) :: gx
      real   (kind=R8), intent(out), allocatable, dimension(:) :: gy
      real   (kind=R8), intent(out), allocatable, dimension(:) :: gw
      real   (kind=R8), intent(out), allocatable, dimension(:) :: tab_dnq
      integer(kind=I4), intent(out)                            :: ng

         real   (kind=R8) :: x1, x2, y1, y2
         integer(kind=I4) :: i, k, nb_gauss_1d

         nb_gauss_1d = 2
         ng          = nb_gauss_1d**2

         allocate( gx(1:ng), gy(1:ng), gw(1:ng) )

         select case(nb_gauss_1d)
            case(1)
               x1 = 0._R8
               y1 = 2._R8

               gx(1:ng) = [   x1]
               gy(1:ng) = [   x1]
               gw(1:ng) = [y1*y1]
            case(2)
               x1 = sqrt(1._R8/3._R8)
               y1 = UN

               gx(1:ng) = [   -x1,    -x1,    +x1,    +x1]
               gy(1:ng) = [   -x1,    +x1,    -x1,    +x1]
               gw(1:ng) = [ y1*y1,  y1*y1,  y1*y1,  y1*y1]
            case(3)
               x1 = sqrt(3._R8/5.0_R8)
               x2 = 0._R8
               y1 = 5._R8/9._R8
               y2 = 8._R8/9._R8

               gx(1:ng) = [   -x1,   -x1,   -x1,    x2,    x2,    x2,   +x1,    +x1,   +x1]
               gy(1:ng) = [   -x1,    x2,   +x1,   -x1,    x2,   +x1,   -x1,     x2,   +x1]
               gw(1:ng) = [ y1*y1, y1*y2, y1*y1, y1*y2, y2*y2, y1*y2, y1*y1,  y1*y2, y1*y1]
         endselect

         allocate( tab_dnq(1:32*ng) )
         i = 1
         do k = 1, ng

            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 1, 1) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 2, 1) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 3, 1) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 4, 1) ; i = i +1

            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 1, 2) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 2, 2) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 3, 2) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 4, 2) ; i = i +1

            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 1, 3) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 2, 3) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 3, 3) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 4, 3) ; i = i +1

            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 1, 4) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 2, 4) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 3, 4) ; i = i +1
            tab_dnq(i) = dnq_xi_i(gx(k), gy(k), 4, 4) ; i = i +1
            !!
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 1, 1) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 2, 1) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 3, 1) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 4, 1) ; i = i +1

            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 1, 2) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 2, 2) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 3, 2) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 4, 2) ; i = i +1

            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 1, 3) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 2, 3) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 3, 3) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 4, 3) ; i = i +1

            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 1, 4) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 2, 4) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 3, 4) ; i = i +1
            tab_dnq(i) = dnq_et_i(gx(k), gy(k), 4, 4) ; i = i +1

         enddo

      return
      endsubroutine init_aire_hermite

      subroutine calcul_tabd_hermite( tab_in, tab_dx, tab_dy, tab_xy, long, larg, hx, hy)
      !================================================================================================
      !< @note
      !
      !  @endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      integer(kind=I4), intent(in )                            :: long
      integer(kind=I4), intent(in )                            :: larg
      real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in
      real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_dx
      real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_dy
      real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_xy
      real   (kind=R8), intent(in )                            :: hx
      real   (kind=R8), intent(in )                            :: hy

         integer(kind=I4) :: i, im, ip, j, jm, jp
         real   (kind=R8) :: ui, uim, uip, ujm, ujp, upp, ump, upm, umm

         do j = 1, larg

            jm = max(j -1,    1)
            jp = min(j +1, larg)

            do i = 1, long
               im = max(i -1,    1)
               ip = min(i +1, long)

               ui  = tab_in(i , j )
               uim = tab_in(im, j )
               uip = tab_in(ip, j )
               ujm = tab_in(i , jm)
               ujp = tab_in(i , jp)
               upp = tab_in(ip, jp)
               ump = tab_in(im, jp)
               upm = tab_in(ip, jm)
               umm = tab_in(im, jm)

               tab_dx(i, j) = (uip -uim)/(2*hx)
               tab_dy(i, j) = (ujp -ujm)/(2*hy)
               tab_xy(i, j) = (upp -ump -upm +umm)/(4*hx*hy)
            enddo
         enddo

         tab_dx(   1, 1:larg) = ( tab_in(      2, 1:larg) -      &  !
                                  tab_in(      1, 1:larg) )/hx      !
         tab_dx(long, 1:larg) = ( tab_in(long   , 1:larg) -      &  !
                                  tab_in(long -1, 1:larg) )/hx      !

         tab_dy(1:long,    1) = ( tab_in(1:long,       2) -      &  !
                                  tab_in(1:long,       1) )/hy      !
         tab_dy(1:long, larg) = ( tab_in(1:long, larg   ) -      &  !
                                  tab_in(1:long, larg -1) )/hy      !

         tab_xy(   1, 1:larg) = ( tab_dy(      2, 1:larg) -      &  !
                                  tab_dy(      1, 1:larg) )/hx      !
         tab_xy(long, 1:larg) = ( tab_dy(long   , 1:larg) -      &  !
                                  tab_dy(long -1, 1:larg) )/hx      !

         tab_xy(1:long,    1) = ( tab_dx(1:long,       2) -      &  !
                                  tab_dx(1:long,       1) )/hy      !
         tab_xy(1:long, larg) = ( tab_dx(1:long, larg   ) -      &  !
                                  tab_dx(1:long, larg -1) )/hy      !

      return
      endsubroutine calcul_tabd_hermite

      subroutine calcul_aire_hermite(tab_in, long, larg, gw, tab_dnq, ng, hx, hy, width, height, aire)
      !================================================================================================
      !< @note
      !
      !  @endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      integer(kind=I4), intent(in )                            :: long
      integer(kind=I4), intent(in )                            :: larg
      integer(kind=I4), intent(in )                            :: ng
      real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in
      real   (kind=R8), intent(in ), dimension(1:ng)           :: gw
      real   (kind=R8), intent(in ), dimension(1:32*ng)        :: tab_dnq
      real   (kind=R8), intent(in )                            :: hx
      real   (kind=R8), intent(in )                            :: hy
      real   (kind=R8), intent(in )                            :: width
      real   (kind=R8), intent(in )                            :: height
      real   (kind=R8), intent(out)                            :: aire

         integer(kind=I4) :: i, j, k, i1, i2, j1, j2
         real   (kind=R8) :: aire_tmp

         real(kind=R8), allocatable, dimension(:) :: dfx
         real(kind=R8), allocatable, dimension(:) :: dfy

         real(kind=R8), allocatable, dimension(:,:) :: tab_dx     ! new grid function evaluations
         real(kind=R8), allocatable, dimension(:,:) :: tab_dy     ! new grid function evaluations
         real(kind=R8), allocatable, dimension(:,:) :: tab_xy     ! new grid function evaluations

         allocate( dfx(1:ng),    &  !
                   dfy(1:ng) )      !

         allocate( tab_dx   (1:long, 1:larg),   &  !
                   tab_dy   (1:long, 1:larg),   &  !
                   tab_xy   (1:long, 1:larg) )     !

         call calcul_tabd_hermite( tab_in = tab_in(1:long, 1:larg),  &  !
                                   tab_dx = tab_dx(1:long, 1:larg),  &  !
                                   tab_dy = tab_dy(1:long, 1:larg),  &  !
                                   tab_xy = tab_xy(1:long, 1:larg),  &  !
                                     long = long,                    &  !
                                     larg = larg,                    &  !
                                       hx = hx,                      &  !
                                       hy = hy )                        !


         aire_tmp = 0._R8
         do j = 1, larg -1

            j1 = j ; j2 = j +1

            do i = 1, long -1

               i1 = i ; i2 = i +1

               do k = 1, ng

                  dfx(k) = (2._R8/hx)*(                              &  !
                  tab_dnq(32*(k -1) +01)*tab_in(i1, j1) +            &  !  u1
                  tab_dnq(32*(k -1) +02)*tab_in(i2, j1) +            &  !  u2
                  tab_dnq(32*(k -1) +03)*tab_in(i2, j2) +            &  !  u3
                  tab_dnq(32*(k -1) +04)*tab_in(i1, j2) +            &  !  u4

                  tab_dnq(32*(k -1) +05)*tab_dx(i1, j1)*hx/2 +       &  ! du1/dx
                  tab_dnq(32*(k -1) +06)*tab_dx(i2, j1)*hx/2 +       &  ! du2/dx
                  tab_dnq(32*(k -1) +07)*tab_dx(i2, j2)*hx/2 +       &  ! du3/dx
                  tab_dnq(32*(k -1) +08)*tab_dx(i1, j2)*hx/2 +       &  ! du4/dx

                  tab_dnq(32*(k -1) +09)*tab_dy(i1, j1)*hy/2 +       &  ! du1/dy
                  tab_dnq(32*(k -1) +10)*tab_dy(i2, j1)*hy/2 +       &  ! du2/dy
                  tab_dnq(32*(k -1) +11)*tab_dy(i2, j2)*hy/2 +       &  ! du3/dy
                  tab_dnq(32*(k -1) +12)*tab_dy(i1, j2)*hy/2 +       &  ! du4/dy

                  tab_dnq(32*(k -1) +13)*tab_xy(i1, j1)*hx*hy/4 +    &  ! du1/dxdy
                  tab_dnq(32*(k -1) +14)*tab_xy(i2, j1)*hx*hy/4 +    &  ! du2/dxdy
                  tab_dnq(32*(k -1) +15)*tab_xy(i2, j2)*hx*hy/4 +    &  ! du3/dxdy
                  tab_dnq(32*(k -1) +16)*tab_xy(i1, j2)*hx*hy/4 )       ! du4/dxdy

                  dfy(k) = (2._R8/hy)*(                              &  !
                  tab_dnq(32*(k -1) +17)*tab_in(i1, j1) +            &  !  u1
                  tab_dnq(32*(k -1) +18)*tab_in(i2, j1) +            &  !  u2
                  tab_dnq(32*(k -1) +19)*tab_in(i2, j2) +            &  !  u3
                  tab_dnq(32*(k -1) +20)*tab_in(i1, j2) +            &  !  u4

                  tab_dnq(32*(k -1) +21)*tab_dx(i1, j1)*hx/2 +       &  ! du1/dx
                  tab_dnq(32*(k -1) +22)*tab_dx(i2, j1)*hx/2 +       &  ! du2/dx
                  tab_dnq(32*(k -1) +23)*tab_dx(i2, j2)*hx/2 +       &  ! du3/dx
                  tab_dnq(32*(k -1) +24)*tab_dx(i1, j2)*hx/2 +       &  ! du4/dx

                  tab_dnq(32*(k -1) +25)*tab_dy(i1, j1)*hy/2 +       &  ! du1/dy
                  tab_dnq(32*(k -1) +26)*tab_dy(i2, j1)*hy/2 +       &  ! du2/dy
                  tab_dnq(32*(k -1) +27)*tab_dy(i2, j2)*hy/2 +       &  ! du3/dy
                  tab_dnq(32*(k -1) +28)*tab_dy(i1, j2)*hy/2 +       &  ! du4/dy

                  tab_dnq(32*(k -1) +29)*tab_xy(i1, j1)*hx*hy/4 +    &  ! du1/dxdy
                  tab_dnq(32*(k -1) +30)*tab_xy(i2, j1)*hx*hy/4 +    &  ! du2/dxdy
                  tab_dnq(32*(k -1) +31)*tab_xy(i2, j2)*hx*hy/4 +    &  ! du3/dxdy
                  tab_dnq(32*(k -1) +32)*tab_xy(i1, j2)*hx*hy/4 )       ! du4/dxdy

               enddo

               do k = 1, ng
                  aire_tmp = aire_tmp + gw(k) * sqrt( UN + dfx(k)**2 + dfy(k)**2 )
               enddo

            enddo

         enddo

         aire = aire_tmp*( hx/2 )*( hy/2 )
         aire = aire/(width*height)

         deallocate( tab_dx   ,  &  !
                     tab_dy   ,  &  !
                     tab_xy )       !

         deallocate( dfx, dfy )

      return
      endsubroutine calcul_aire_hermite

   endsubroutine calcul_asfc_hermite


   subroutine calcul_asfc(tab_in, scal, asfc_res, omp)
   !================================================================================================
   !! Return the *asfc* of a surface regarding the default parameter *method_asfc*
   implicit none
   type(SCALE_SURF), intent(in )                                      :: scal     !! *surface characteristics*
   real   (kind=R8), intent(in ), dimension(1:scal%xres, 1:scal%yres) :: tab_in   !! *input surface*
   real   (kind=R8), intent(out), dimension(1:2)                      :: asfc_res !! *result: asfc, adjustment factor*
   logical(kind=I4), intent(in )                                      :: omp      !! *with openmp ?*

      if (out_lin) call get_unit(unit_out_lin)
      if (out_spl) call get_unit(unit_out_spl)

      select case(method_asfc)
         case(lin_all)
            call calcul_asfc_lin_all(tab_in, scal, asfc_res)
         case(spl_all)
            call calcul_asfc_spl_all(tab_in, scal, asfc_res)
         case(hermite)
            call calcul_asfc_hermite(tab_in, scal, asfc_res, omp)
         case default
            stop 'no valid method'
      endselect

   return
   endsubroutine calcul_asfc


   integer(kind=I4) function locate(n, xx, x)
   !================================================================================================
   !< @note Function that returns the location of an element in a vector.
   !<
   !< Given an array xx(1: n) , and given a value x , it returns a value j such that x is between
   !< xx( j ) and xx( j + 1 ).
   !<
   !< xx must be monotonic, either increasing or decreasing. j = 0 or j = n is returned to indicate
   !< that x is out of range.
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in)                 :: n  !! *vector length*
   real   (kind=R8), intent(in)                 :: x  !! *value to locate*
   real   (kind=R8), intent(in), dimension(1:n) :: xx !! *vector*

      integer(kind=I4) :: jl, jm, ju
      logical(kind=I4) :: ascnd
      ascnd = (xx(n) >= xx(1))               !  true if ascending order of table, false otherwise.

      jl = 0                                 ! initialize lower
      ju = n +1                              !  and upper limits.
      do
         if (ju -jl <= 1) exit               ! repeat until this condition is satisfied.
         jm = (ju +jl)/2                     ! compute a midpoint,
         if (ascnd .eqv. (x >= xx(jm))) then
            jl = jm                          ! and replace either the lower limit
         else
            ju = jm                          ! or the upper limit, as appropriate.
         endif
      enddo
      if ( x <= xx(1) ) then                 ! then set the output, being careful with the endpoints.
         locate = 1
      elseif ( x >= xx(n) ) then
         locate = n -1
      else
         locate = jl
      endif
   return
   endfunction locate


   integer(kind=I4) function locate2(n, xx, x, eps)
   !================================================================================================
   !< @note Function that returns the location of an element in a vector.
   !<
   !< Given an array xx(1: n) , and given a value x , it returns a value j such that x is between
   !< xx( j ) and xx( j + 1 ).
   !<
   !< xx must be monotonic, either increasing or decreasing. j = 0 or j = n is returned to indicate
   !< that x is out of range.
   !<
   !< The difference with [[locate]] is the use of *eps* for comparisons of reals.
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in)                 :: n     !! *vector length*
   real   (kind=R8), intent(in)                 :: x     !! *value to locate*
   real   (kind=R8), intent(in)                 :: eps   !! *small value for comparisons of reals*
   real   (kind=R8), intent(in), dimension(1:n) :: xx    !! *vector*

      integer(kind=I4) :: jl, jm, ju

      jl = 0                                 ! initialize lower
      ju = n +1                              !  and upper limits.
      do
         if (ju -jl <= 1) exit               ! repeat until this condition is satisfied.
         jm = (ju +jl)/2                     ! compute a midpoint,
         if ( x -xx(jm) > -eps ) then
            jl = jm                          ! and replace either the lower limit
         else
            ju = jm                          ! or the upper limit, as appropriate.
         endif
      enddo
      if ( abs(x -xx(1))<eps ) then          ! then set the output, being careful with the endpoints.
         locate2 = 1
      elseif ( abs(x -xx(n))<eps ) then
         locate2 = n -1
      else
         locate2 = jl
      endif
   return
   endfunction locate2


   subroutine init_beta_boltz(bt, n_bt, v_x, v_y, n_pt)
   !================================================================================================
   !< @note Function that initializes the fitting *tanh* function [[f_boltz]] parameters.
   !<
   !<  \[ f_{boltz}(x_i)=\beta_2 + \beta_3 \tanh \left( \dfrac{x_i -\beta_1}{\beta_4} \right) \]
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in)                     :: n_bt !! *number of parameters*
   integer(kind=I4), intent(in)                     :: n_pt !! *data vector length*
   real   (kind=R8), intent(out), dimension(1:n_bt) :: bt   !! *vector \(\beta\) of parameters*
   real   (kind=R8), intent(in ), dimension(1:n_pt) :: v_x  !! *x data to fit*
   real   (kind=R8), intent(in ), dimension(1:n_pt) :: v_y  !! *y data to fit*

      real(kind=R8) :: a, pente

      bt(1) = 0.!v_x(1) +(v_x(n_pt) -v_x(1))/3

      bt(2) = ( sum( v_y(n_pt - 9 : n_pt) ) + sum( v_y(1:10) ) )/ (2 * 10)

      bt(3) = ( sum( v_y(n_pt - 1 : n_pt) ) - sum( v_y(1:02) ) )/ (2 * 2 )   !;  a = bt(3)

      !pente = (v_y(n_pt/2) -v_y(n_pt/4))/(v_x(n_pt/2) -v_x(n_pt/4))

      bt(4) = 1.!a/pente
   return
   endsubroutine init_beta_boltz


   real(kind=R8) function f_boltz(xi, beta, n_beta)
   !================================================================================================
   !< @note Fitting function.
   !<
   !<  \[ f_{boltz}(x_i)=\beta_2 + \beta_3 \tanh \left( \dfrac{x_i -\beta_1}{\beta_4} \right) \]
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                      :: n_beta !! *number of parameters*
   real   (kind=R8), intent(in   )                      :: xi     !! *data point*
   real   (kind=R8), intent(inout), dimension(1:n_beta) :: beta   !! *parameter vector*

      real(kind=R8) :: x, x0, y0, a, b

      x0 = beta(1) ; y0 = beta(2) ; a = beta(3) ; b = beta(4) !; c = beta(5) ; d = beta(6)
      x  = xi -x0
      f_boltz = y0 +a*tanh( x/b ) !+c*x +d*x*x
   return
   endfunction f_boltz


   real(kind=R8) function df_boltz(xi, beta, n_beta, ivar)
   !================================================================================================
   !! Fitting function partial derivatives.
   implicit none
   integer(kind=I4), intent(in   )                       :: n_beta   !! *number of parameters*
   integer(kind=I4), intent(in   )                       :: ivar     !! *parameter number*
   real   (kind=R8), intent(in   )                       :: xi       !! *data point*
   real   (kind=R8), intent(inout), dimension(1:n_beta)  :: beta     !! *parameter vector*

      real(kind=R8) :: x0, y0, a, b, x, th

      x0 = beta(1) ; y0 = beta(2) ; a = beta(3) ; b = beta(4) !; c = beta(5) ; d = beta(6)
      x  = xi -x0 ; th = tanh( x/b )
      select case(ivar)
         case(0)
            df_boltz = +(a/b)*(UN -th**2) !+c +2*d*x  ! special case: derivative regarding xi
         case(1)
            df_boltz = -(a/b)*(UN -th**2)
         case(2)
            df_boltz = UN
         case(3)
            df_boltz = th
         case(4)
            df_boltz = -(a*x/b**2)*(UN -th**2)
         !case(5)
         !   df_boltz = x
         !case(6)
         !   df_boltz = x*x
         case default
            stop 'df_boltz, bad choice'
      endselect
   return
   endfunction df_boltz


   real(kind=R8) function n_i(xi, i)
   !================================================================================================
   !! 1D shape function, quadratic case.
   implicit none
   real   (kind=R8), intent(in) :: xi
   integer(kind=I4), intent(in) :: i

      select case(i)
         case(1)
            n_i = 0.250_R8 *( (+UN -xi)**2 )*(2._R8 +xi)
         case(2)
            n_i = 0.250_R8 *( (+UN -xi)**2 )*(   UN +xi)
         case(3)
            n_i = 0.250_R8 *( (+UN +xi)**2 )*(2._R8 -xi)
         case(4)
            n_i = 0.250_R8 *( (+UN +xi)**2 )*(  -UN +xi)
         case default
            stop 'n_i, bad node ddl'
      endselect
   return
   endfunction n_i


   real(kind=R8) function dn_i(xi, i)
   !================================================================================================
   !! 1D shape function derivative, quadratic case.
   implicit none
   real   (kind=R8), intent(in) :: xi
   integer(kind=I4), intent(in) :: i

      select case(i)
         case(1)
            dn_i = 0.250_R8 *( -3*(+UN -xi**2)            )
         case(2)
            dn_i = 0.250_R8 *(    (+UN -xi   )*(-UN -3*xi))
         case(3)
            dn_i = 0.250_R8 *( +3*(+UN -xi**2)            )
         case(4)
            dn_i = 0.250_R8 *(    (+UN +xi   )*(-UN +3*xi))
         case default
            stop 'dn_i, bad node ddl'
      endselect
   return
   endfunction dn_i


   real(kind=R8) function nq_i(xi, et, i, j)
   !================================================================================================
   !! 2D shape function, quadratic case.
   implicit none
   real   (kind=R8), intent(in) :: xi, et
   integer(kind=I4), intent(in) :: i, j

      integer(kind=I4) :: k

      k = 4*(i -1) +j
      select case(k)
         case( 1)
            nq_i = n_i(xi, 1)*n_i(et, 1) ! node 1 : u
         case( 2)
            nq_i = n_i(xi, 2)*n_i(et, 1) ! node 1 : du/dx
         case( 3)
            nq_i = n_i(xi, 1)*n_i(et, 2) ! node 1 : du/dy
         case( 4)
            nq_i = n_i(xi, 2)*n_i(et, 2) ! node 1 : d2u/dxdy

         case( 5)
            nq_i = n_i(xi, 3)*n_i(et, 1) ! node 2 : u
         case( 6)
            nq_i = n_i(xi, 4)*n_i(et, 1) ! node 2 : du/dx
         case( 7)
            nq_i = n_i(xi, 3)*n_i(et, 2) ! node 2 : du/dy
         case( 8)
            nq_i = n_i(xi, 4)*n_i(et, 2) ! node 2 : d2u/dxdy

         case( 9)
            nq_i = n_i(xi, 3)*n_i(et, 3) ! node 3 : u
         case(10)
            nq_i = n_i(xi, 4)*n_i(et, 3) ! node 3 : du/dx
         case(11)
            nq_i = n_i(xi, 3)*n_i(et, 4) ! node 3 : du/dy
         case(12)
            nq_i = n_i(xi, 4)*n_i(et, 4) ! node 3 : d2u/dxdy

         case(13)
            nq_i = n_i(xi, 1)*n_i(et, 3) ! node 4 : u
         case(14)
            nq_i = n_i(xi, 2)*n_i(et, 3) ! node 4 : du/dx
         case(15)
            nq_i = n_i(xi, 1)*n_i(et, 4) ! node 4 : du/dy
         case(16)
            nq_i = n_i(xi, 2)*n_i(et, 4) ! node 4 : d2u/dxdy
         case default
            stop 'nq_i, bad node ddl'
      endselect
   return
   endfunction nq_i


   real(kind=R8) function dnq_xi_i(xi, et, i, j)
   !================================================================================================
   !! 2D shape function \(\xi\) derivative, quadratic case.
   implicit none
   real   (kind=R8), intent(in) :: xi, et
   integer(kind=I4), intent(in) :: i, j

      integer(kind=I4) :: k

      k = 4*(i -1) +j
      select case(k)
         case( 1)
            dnq_xi_i = dn_i(xi, 1)*n_i(et, 1)
         case( 2)
            dnq_xi_i = dn_i(xi, 2)*n_i(et, 1)
         case( 3)
            dnq_xi_i = dn_i(xi, 1)*n_i(et, 2)
         case( 4)
            dnq_xi_i = dn_i(xi, 2)*n_i(et, 2)

         case( 5)
            dnq_xi_i = dn_i(xi, 3)*n_i(et, 1)
         case( 6)
            dnq_xi_i = dn_i(xi, 4)*n_i(et, 1)
         case( 7)
            dnq_xi_i = dn_i(xi, 3)*n_i(et, 2)
         case( 8)
            dnq_xi_i = dn_i(xi, 4)*n_i(et, 2)

         case( 9)
            dnq_xi_i = dn_i(xi, 3)*n_i(et, 3)
         case(10)
            dnq_xi_i = dn_i(xi, 4)*n_i(et, 3)
         case(11)
            dnq_xi_i = dn_i(xi, 3)*n_i(et, 4)
         case(12)
            dnq_xi_i = dn_i(xi, 4)*n_i(et, 4)

         case(13)
            dnq_xi_i = dn_i(xi, 1)*n_i(et, 3)
         case(14)
            dnq_xi_i = dn_i(xi, 2)*n_i(et, 3)
         case(15)
            dnq_xi_i = dn_i(xi, 1)*n_i(et, 4)
         case(16)
            dnq_xi_i = dn_i(xi, 2)*n_i(et, 4)
         case default
            stop 'dnq_xi_i, bad node ddl'
      endselect
   return
   endfunction dnq_xi_i


   real(kind=R8) function dnq_et_i(xi, et, i, j)
   !================================================================================================
   !! 2D shape function \(\eta\) derivative, quadratic case.
   implicit none
   real   (kind=R8), intent(in) :: xi, et
   integer(kind=I4), intent(in) :: i, j

      integer(kind=I4) :: k

      k = 4*(i -1) +j
      select case(k)
         case( 1)
            dnq_et_i = n_i(xi, 1)*dn_i(et, 1)
         case( 2)
            dnq_et_i = n_i(xi, 2)*dn_i(et, 1)
         case( 3)
            dnq_et_i = n_i(xi, 1)*dn_i(et, 2)
         case( 4)
            dnq_et_i = n_i(xi, 2)*dn_i(et, 2)

         case( 5)
            dnq_et_i = n_i(xi, 3)*dn_i(et, 1)
         case( 6)
            dnq_et_i = n_i(xi, 4)*dn_i(et, 1)
         case( 7)
            dnq_et_i = n_i(xi, 3)*dn_i(et, 2)
         case( 8)
            dnq_et_i = n_i(xi, 4)*dn_i(et, 2)

         case( 9)
            dnq_et_i = n_i(xi, 3)*dn_i(et, 3)
         case(10)
            dnq_et_i = n_i(xi, 4)*dn_i(et, 3)
         case(11)
            dnq_et_i = n_i(xi, 3)*dn_i(et, 4)
         case(12)
            dnq_et_i = n_i(xi, 4)*dn_i(et, 4)

         case(13)
            dnq_et_i = n_i(xi, 1)*dn_i(et, 3)
         case(14)
            dnq_et_i = n_i(xi, 2)*dn_i(et, 3)
         case(15)
            dnq_et_i = n_i(xi, 1)*dn_i(et, 4)
         case(16)
            dnq_et_i = n_i(xi, 2)*dn_i(et, 4)
         case default
            stop 'dnq_xi_i, bad node ddl'
      endselect
   return
   endfunction dnq_et_i


   real(kind=R8) function dnq_xi_et_i(xi, et, i, j)
   !================================================================================================
   !! 2D shape function \(\xi\), \(\eta\) derivative, quadratic case.
   implicit none
   real   (kind=R8), intent(in) :: xi, et
   integer(kind=I4), intent(in) :: i, j

      integer(kind=I4) :: k

      k = 4*(i -1) +j
      select case(k)
         case( 1)
            dnq_xi_et_i = dn_i(xi, 1)*dn_i(et, 1)
         case( 2)
            dnq_xi_et_i = dn_i(xi, 2)*dn_i(et, 1)
         case( 3)
            dnq_xi_et_i = dn_i(xi, 1)*dn_i(et, 2)
         case( 4)
            dnq_xi_et_i = dn_i(xi, 2)*dn_i(et, 2)

         case( 5)
            dnq_xi_et_i = dn_i(xi, 3)*dn_i(et, 1)
         case( 6)
            dnq_xi_et_i = dn_i(xi, 4)*dn_i(et, 1)
         case( 7)
            dnq_xi_et_i = dn_i(xi, 3)*dn_i(et, 2)
         case( 8)
            dnq_xi_et_i = dn_i(xi, 4)*dn_i(et, 2)

         case( 9)
            dnq_xi_et_i = dn_i(xi, 3)*dn_i(et, 3)
         case(10)
            dnq_xi_et_i = dn_i(xi, 4)*dn_i(et, 3)
         case(11)
            dnq_xi_et_i = dn_i(xi, 3)*dn_i(et, 4)
         case(12)
            dnq_xi_et_i = dn_i(xi, 4)*dn_i(et, 4)

         case(13)
            dnq_xi_et_i = dn_i(xi, 1)*dn_i(et, 3)
         case(14)
            dnq_xi_et_i = dn_i(xi, 2)*dn_i(et, 3)
         case(15)
            dnq_xi_et_i = dn_i(xi, 1)*dn_i(et, 4)
         case(16)
            dnq_xi_et_i = dn_i(xi, 2)*dn_i(et, 4)
         case default
            stop 'dnq_xi_i, bad node ddl'
      endselect
   return
   endfunction dnq_xi_et_i


   subroutine indice_fractal(tab_in, long, larg, indf)
   !================================================================================================
   !! Function that returns the fractal dimension with the box counting method
   implicit none
   integer(kind=I4), intent(in )                            :: long   !! *surface array length*
   integer(kind=I4), intent(in )                            :: larg   !! *surface array width*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in !! *surface array*
   real   (kind=R8), intent(out), dimension(3)              :: indf   !! *result: ordinate at origin, slope, R2*

      integer(kind=I4) :: i, j, k, ib, nbhmax, nni, nb, dec_i, dec_j, ri, rj
      real   (kind=R8) :: lboite, hhmax, hhmin, ddh, t1, t2, t3, t4, hmax, hmin
      real   (kind=R8) :: yibarr, yi_m_yichap, yi_m_yibarr, r2adj

      real(kind=R8), dimension(1:2) :: var

      real   (kind=R8), allocatable, dimension(:,:) :: tab, Jf
      real   (kind=R8), allocatable, dimension(:)   :: tab_nni

      integer(kind=I4), allocatable, dimension(:,:) :: pas_i, pas_j
      integer(kind=I4), allocatable, dimension(:)   :: nnb

      type(moment_stat)                             :: mom

      call calc_moments(   tab = reshape( tab_in(1:long, 1:larg), [long * larg] ),  &  !
                            mx = mom,                                               &  !
                        nb_mom = 2 )                                                   !

      allocate(tab(1:long, 1:larg))
      tab(1:long, 1:larg) = (tab_in(1:long, 1:larg) -mom%mu)/mom%si

      hhmin = minval(tab(1:long, 1:larg))     ! hauteur min de la surface
      tab(1:long, 1:larg) = tab(1:long, 1:larg) +hhmin +1._R8

      hhmax = maxval(tab(1:long, 1:larg))     ! hauteur max de la surface
      hhmin = minval(tab(1:long, 1:larg))     ! hauteur min de la surface
      ddh   = hhmax -hhmin                    ! amplitude

      nbhmax = nint( log( UN*min(long-1,larg-1) )/log(2.) )      ! nbre max de comptages, attention il y a n-1 intervalles pour n points
      if (2**nbhmax > min(long-1, larg-1)) nbhmax = nbhmax -1

      allocate(      Jf(1:nbhmax, 1:2) )
      allocate( tab_nni(1:nbhmax) )

      nb = 2**nbhmax ! nombre de boîtes dans une direction
      allocate( pas_i(1:nb, 1:nbhmax), pas_j(1:nb, 1:nbhmax) )
      allocate(   nnb(1:nb) )

      pas_i = 0
      pas_j = 0

      pas_i(1:nb, nbhmax) = (long -1)/nb ! longueur moyenne d'UN pas selon x (le pas n'est pas forcément constant : nbpts /= 2**n +1, par exemple)
      pas_j(1:nb, nbhmax) = (larg -1)/nb ! id selon y

      ri = mod(long -1, nb) ! si la division au-dessus ne tombe pas juste, c'est le reste selon x
      rj = mod(larg -1, nb) ! id selon y

      if ( ri>0 ) then ! s'il y a un résidu, répartition régulière de ce résidu sur le découpage le plus fin
         do i = 1, ri
            pas_i( (i-1)*(nb/ri) +1, nbhmax ) = pas_i( (i-1)*(nb/ri) +1, nbhmax ) +1
         enddo
      endif
      if ( rj>0 ) then
         do j = 1, rj
            pas_j( (j-1)*(nb/rj) +1, nbhmax ) = pas_j( (j-1)*(nb/rj) +1, nbhmax ) +1
         enddo
      endif

      do ib = nbhmax-1, 1, -1 ! agglomération des pas 2 à 2 pour former le pas de la boîte englobante
         do k = 1, 2**ib
            pas_i(k, ib) = pas_i(2*(k-1)+1, ib+1) +pas_i(2*k, ib+1)
            pas_j(k, ib) = pas_j(2*(k-1)+1, ib+1) +pas_j(2*k, ib+1)
         enddo
      enddo

      do ib = 1, nbhmax               ! niveau de découpage
         nb = 2**ib                   ! nombre de boîtes dans une direction pour ce découpage
         lboite = (ddh -2*EPS_R8)/nb  ! taille z de la boîte

         nni   = 0
         dec_i = 1
         do i = 1, nb ! numéro i de la colonne de boîtes

            dec_j = 1
            do j = 1, nb ! numéro j de la colonne de boîtes
               ! on considère le plan résultant de la coupe de la colonne par la surface
               t1 = tab(dec_i              , dec_j              )
               t2 = tab(dec_i +pas_i(i, ib), dec_j              )
               t3 = tab(dec_i +pas_i(i, ib), dec_j +pas_j(j, ib))
               t4 = tab(dec_i              , dec_j +pas_j(j, ib))
               ! ce plan traverse plusieurs boîtes, il suffit de les compter
               hmax  = max(t1, t2, t3, t4)
               hmin  = min(t1, t2, t3, t4)
               nni   = nni + ceiling( hmax/lboite ) &
                           -   floor( hmin/lboite )

               dec_j = dec_j +pas_j(j, ib)
            enddo

            dec_i = dec_i +pas_i(i, ib)
         enddo

         tab_nni(ib) = log(UN*nni)
         Jf(ib, 1)   = UN
         Jf(ib, 2)   = log(UN/lboite)

      enddo

      call moindres_carres_lineaire( nb_var = 2,                   & !
                                     nb_pts = nbhmax,              & !
                                        hij = tab_nni(1:nbhmax),   & !
                                       beta = var(1:2),            & !
                                         Jf = Jf(1:nbhmax, 1:2) )    !

      yibarr = sum( tab_nni(1:nbhmax) )/(nbhmax)
      yi_m_yichap = 0
      yi_m_yibarr = 0
      do i = 1, nbhmax
         yi_m_yichap = yi_m_yichap +( tab_nni(i) -( Jf(i, 1)*var(1) +Jf(i, 2)*var(2) ) )**2
         yi_m_yibarr = yi_m_yibarr +( tab_nni(i) -yibarr )**2
      enddo
      r2adj = UN -(yi_m_yichap/(nbhmax -2))/(yi_m_yibarr/(nbhmax -1))

      indf = [var(2), var(1), r2adj] ! fractal index first

      deallocate( tab, tab_nni, Jf, pas_i, pas_j, nnb )

   return
   endsubroutine indice_fractal


endmodule asfc
