!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: april, 08 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Subroutines for anisotropy detection in surfaces**
!<  </span>
module anisotropy
use data_arch,     only : I4, R4, R8, PI_R8, EPS_R8, UN
use miscellaneous, only : trans_corner2center, trans_center2corner, get_unit
use sort_arrays,   only : sort_array2
use tchebychev,    only : least_squares_tcheby
use fftw3,         only : apod, fftw_plan_with_nthreads, PAD_FFT, extend,     &  !
                          calc_fftw3, tab_calc_fftw3, FORWARD, BACKWARD
use stat_mom,      only : calc_moments, MOMENT_STAT
use filter,        only : fft_filter
use surfile,       only : init_scal, write_surf, SCALE_SURF
!$ use omp_lib
implicit none

   real(kind=R8)    :: PAD_FFT_ANI !! *dimension multiplier for 0-padding*

   character(len=6) :: APO_FFT_ANI !! *dimension multiplier for 0-padding*

!> {!ANISO/src/inc_doc/analyses.md!}
!> {!css/button.html!}

contains


   subroutine simple_anisotropy(tabin, long, larg, vec_len, vec_pks, vec_slp, scale_xy, multi_fft)
   !================================================================================================
   !< @note Function that returns some anisotropy parameters calculated on a polar representation,
   !< for each angle from 0 to 179
   !<
   !< + the excess of length
   !< + the RMS slope
   !< + the mean peak width
   !<
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                             :: long       !! *surface length*
   integer(kind=I4), intent(in )                             :: larg       !! *surface width*
   logical(kind=I4), intent(in )                             :: multi_fft  !! *parallel ffts?*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg)  :: tabin      !! *surface acf array*
   real   (kind=R8), intent(in ), dimension(1:2)             :: scale_xy   !! *lag along x, y*
   real   (kind=R8), intent(out), dimension(0:179)           :: vec_len    !! *vector containing path length ratios*
   real   (kind=R8), intent(out), dimension(0:179)           :: vec_pks    !! *vector containing peak mean width*
   real   (kind=R8), intent(out), dimension(0:179)           :: vec_slp    !! *vector containing RMS slopes*

      integer(kind=I4) :: lo, la, ll, nx, ny, ird
      integer(kind=I4) :: p, q, k, nb_cross, nb_peak
      real   (kind=R8) :: scx, scy, dr, s, c
      real   (kind=R8) :: mx, my, reg_a, reg_b
      real   (kind=R8) :: h1, h2, h3, h4, hh
      real   (kind=R8) :: inc_a, theta, x, y, xb, yb, xm, ym, xp, yp, tmp

      real   (kind=R8), allocatable, dimension(:,:) :: height_disc

      real   (kind=R8), allocatable, dimension(:)   :: vec_tmp, slope

      integer(kind=I4), allocatable, dimension(:)   :: cross, peaks

      ! lags in micron
      scx = scale_xy(1)  ! x lag
      scy = scale_xy(2)  ! y lag

      ! define the surface center location
      if ( long == 2 * (long / 2) ) then ; lo = long/2 + 1 ; else ; lo = long/2 ; endif
      if ( larg == 2 * (larg / 2) ) then ; la = larg/2 + 1 ; else ; la = larg/2 ; endif

      ! define the square length embedded in the surface
      ll = min(lo, la) - 1

      allocate( height_disc(0:ll, 0:359) )

      allocate( vec_tmp(0:2*ll) )
      allocate( cross(1:2*ll) )
      allocate( slope(1:2*ll) )
      allocate( peaks(1:2*ll) )

      ! angle increment
      inc_a = 2 * PI_R8 / 360

      ! determination of heights on a diameter of the rosette, obtained by linear interpolation:
      ! a point "falls" into a rectangular element [h1,h2,h3,h4], its height is determined by linear interpolation.
      !* Calculate heights along the radius for each angular position.

      do p = 0, ll                ! point number on a radius
         ird = p                  ! point radius

         do q = 0, 359            ! point angular location (°)
            theta = q * inc_a     ! angular location (rad)

            ! projection on x and y of the point determined by its radius and angle
            ! by taking the floor, we get the number of the lower row and the left column of the rectangle
            ! the remainder (x-nx) thus represents the abscissa of the point in the rectangle with sides=1

            x = lo + ird * cos(theta)

            if ( abs(x-nint(x)) < 1.e-3_R8) then

               nx = nint(x)

            else

               nx = floor(x)

            endif
            xb = x - nx

            y = la + ird * sin(theta)

            if ( abs(y-nint(y)) < 1.e-3_R8) then

               ny = nint(y)

            else

               ny = floor(y)

            endif
            yb = y - ny

            xm = UN - xb ; xp = xb
            ym = UN - yb ; yp = yb

            if ( nx+1 <= long .and. ny+1 <= larg .and.      &  !
                 nx   >= 1    .and. ny   >= 1) then
               ! note ird may be greater than lo or la
               h1 = tabin(nx    , ny    )
               h2 = tabin(nx + 1, ny    )
               h3 = tabin(nx + 1, ny + 1)
               h4 = tabin(nx    , ny + 1)

               hh = h1 * xm * ym +  &  !
                    h2 * xp * ym +  &  !
                    h3 * xp * yp +  &  !
                    h4 * xm * yp       !

               height_disc(p, q) = hh

            endif

         enddo ! q = 0, 359

      enddo ! p = 0, ll

      ! centering
      tmp = sum( height_disc(0:ll, 0:359) ) / ( (ll + 1) * 360 )

      height_disc(0:ll, 0:359) = height_disc(0:ll, 0:359) - tmp

      ! for each degree, determine the mean peak width and the path length ratio
      do q = 0, 359 - 180

         do p = ll, 1, -1 ; vec_tmp(ll - p) = height_disc(p, q + 180) ; enddo ! left part of the height vector
         do p = 0, ll, +1 ; vec_tmp(ll + p) = height_disc(p, q      ) ; enddo ! right part

         theta = q * inc_a ; s = sin(theta) ; c = cos(theta)

         ! element length, usually scx=scy
         dr = sqrt( (s * scx)**2 + (c * scy)**2 )

         ! path relative length  -----------------------------------------------------------------------------------------------
         !
         tmp = 0.
         do p = 1, 2*ll

            tmp = tmp + sqrt( dr**2 + (vec_tmp(p) - vec_tmp(p - 1))**2 )

         enddo

         vec_len(q) = tmp / (2*ll*dr) - 1.

         ! path absolute slope  ------------------------------------------------------------------------------------------------
         !
         do p = 1, 2*ll

            slope(p) = ( vec_tmp(p) - vec_tmp(p - 1) ) / dr

         enddo
         tmp = sum( slope(1:2*ll) ) / (2*ll)
         slope(1:2*ll) = slope(1:2*ll) - tmp

         vec_slp(q) = sqrt( sum( slope(1:2*ll)**2 ) ) / (2*ll)

         ! find each abscissa point where the height crosses z=0 -----------------------------------------------------------------
         !

         !------- subtract mean profile (least square)
         mx = (ll + 1) * dr                       ! x mean x = 0, dr, 2*dr, ..., 2*ll*dr
         my = sum( vec_tmp(0:2*ll) ) / (2*ll + 1) ! y mean

         reg_a = 0
         tmp   = 0
         do p = 0, 2*ll

            reg_a = reg_a + (p * dr - mx) * (vec_tmp(p) - my)
            tmp   = tmp   + (p * dr - mx)**2

         enddo
         reg_a = reg_a / tmp        ! first regressor

         reg_b = my - reg_a * mx    ! second regressor

         do p = 0, 2*ll

            vec_tmp(p) = vec_tmp(p) - (reg_a * p * dr + reg_b)

         enddo

         !------- find peaks
         cross(1:2*ll) = 0
         peaks(1:2*ll) = 0

         k = 1
         cross(k) = 0
         do p = 1, 2*ll

            if ( vec_tmp(p - 1) * vec_tmp(p) < 0. ) then ! the height crosses z=0

               k = k + 1
               cross(k) = p

            endif

         enddo
         k = k + 1
         cross(k) = 2*ll
         nb_cross = k

         ! determine the peak width
         !
         k = 0
         do p = 1, nb_cross - 1

            if ( vec_tmp( (cross(p + 1) + cross(p))/2 ) > 0 ) then ! it must be a peak

               k = k + 1

               peaks(k) = cross(k + 1) - cross(k)

            endif

         enddo
         nb_peak = k

         ! mean peak width
         !
         if ( nb_peak > 0 ) then

            vec_pks(q) = dr * sum( peaks(1:nb_peak) ) / nb_peak

         else

            vec_pks(q) = dr * 2*ll

         endif

      enddo

      deallocate( height_disc, vec_tmp, cross, peaks, slope )

   return
   endsubroutine simple_anisotropy


   subroutine multiple_anisotropy(tabin, long, larg, scale_xy, multi_fft, vec_ani)
   !================================================================================================
   !< @note Function that returns simple_anisotropy min, max and max/min for different Gaussian filter
   !< cutoff
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                                  :: long       !! *surface length*
   integer(kind=I4), intent(in )                                  :: larg       !! *surface width*
   logical(kind=I4), intent(in )                                  :: multi_fft  !! *parallel ffts?*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg)       :: tabin      !! *surface acf array*
   real   (kind=R8), intent(in ), dimension(1:2)                  :: scale_xy   !! *lag along x and y*
   real   (kind=R8), intent(out), dimension(1:9)                  :: vec_ani    !! *anisotropy parameters*

      integer(kind=I4) :: icut, i
      real   (kind=R8) :: dx, dy, fft_cutoff,length

      real   (kind=R8) :: max_mean_pks, min_mean_pks, rat_mean_pks, max_mean_len
      real   (kind=R8) :: min_mean_len, rat_mean_len, max_mean_slp, min_mean_slp, rat_mean_slp

      character(len=2) :: str

      real   (kind=R8), dimension(0:179) :: mean_len, mean_pks, mean_slp

      real   (kind=R8), allocatable, dimension(:,:) :: bf_tab

      real   (kind=R8), allocatable, dimension(:,:) :: mat_len    !! *array containing anisotropy outputs*
      real   (kind=R8), allocatable, dimension(:,:) :: mat_pks    !! *array containing anisotropy outputs*
      real   (kind=R8), allocatable, dimension(:,:) :: mat_slp    !! *array containing anisotropy outputs*

      type(SCALE_SURF) :: scal_surf

      allocate( bf_tab(1:long, 1:larg) )

      dx = scale_xy(1) ! x lag
      dy = scale_xy(2) ! y lag

      length = dx / 1.e-6 ! base cutoff = 1 µm

      allocate( mat_len(0:179, 1:10) )
      allocate( mat_pks(0:179, 1:10) )
      allocate( mat_slp(0:179, 1:10) )

      call init_scal( scal   = scal_surf,       &  !
                      nx     = long,            &  !
                      ny     = larg,            &  !
                      lx     = long * dx,       &  !
                      ly     = larg * dy,       &  !
                      unit_z = 'nm')               !

      do icut = 1, 10

         fft_cutoff = length / icut

         call fft_filter(tab       = tabin(1:long, 1:larg),    &  ! in
                         long      = long,                     &  ! in
                         larg      = larg,                     &  ! in
                         cutoff    = fft_cutoff,               &  ! in
                         bf_tab    = bf_tab(1:long, 1:larg),   &  ! out
                         multi_fft = multi_fft)                   ! in

         write(str,'(i2.2)') icut

         if (.false.) then
            call write_surf( nom_fic = "out/bf_tab"//str//".sur",    &  !
                             tab_s   = bf_tab(1:long, 1:larg),       &  !
                             scal    = scal_surf )                      !
         endif

         call simple_anisotropy( tabin     = bf_tab(1:long, 1:larg),     &  ! in
                                 long      = long,                       &  ! in
                                 larg      = larg,                       &  ! in
                                 vec_len   = mat_len(0:179, icut),       &  ! out   path length ratios
                                 vec_pks   = mat_pks(0:179, icut),       &  ! out   peak mean width
                                 vec_slp   = mat_slp(0:179, icut),       &  ! out   RMS slopes
                                 scale_xy  = scale_xy(1:2),              &  ! in
                                 multi_fft = multi_fft)                     ! in

      enddo

      do i = 0, 179

        mean_pks(i) = sum( mat_pks(i, 06:10) ) / 5. ! mean on smoother profiles
        mean_len(i) = sum( mat_len(i, 01:05) ) / 5.
        mean_slp(i) = sum( mat_slp(i, 01:05) ) / 5.

      enddo

      deallocate( bf_tab, mat_len, mat_pks, mat_slp )

      max_mean_pks = maxval( mean_pks(0:179) )     ! robust maximum of the peak mean width
      min_mean_pks = minval( mean_pks(0:179) )     !        minimum
      rat_mean_pks = max_mean_pks / min_mean_pks   ! ratio

      max_mean_len = maxval( mean_len(0:179) )     ! robust maximum of the path length ratio
      min_mean_len = minval( mean_len(0:179) )     !        minimum
      rat_mean_len = max_mean_len / min_mean_len   ! ratio

      max_mean_slp = maxval( mean_slp(0:179) )     ! robust maximum of the RMS slopes
      min_mean_slp = minval( mean_slp(0:179) )     !        minimum
      rat_mean_slp = max_mean_slp / min_mean_slp   ! ratio

      vec_ani(1:9) = [ max_mean_pks, min_mean_pks, rat_mean_pks,  &  !
                       max_mean_len, min_mean_len, rat_mean_len,  &  !
                       max_mean_slp, min_mean_slp, rat_mean_slp ]    !

   return
   endsubroutine multiple_anisotropy


   subroutine ellipse_acf(tabin, long, larg, p_acv, cut, scale_xy, omp)
   !================================================================================================
   !< @note Function that returns p_acv which contains parameters on anisotropy.
   !<
   !<  - p_acv(1) = axe_a,                                ellipsis big axis
   !<  - p_acv(2) = axe_b,                                ellipsis small axis
   !<  - p_acv(3) = axe_a/axe_b                           another anisotropy factor
   !<  - p_acv(4) = nint(angle/inc_a),                    main texture orientation
   !<  - p_acv(5) = ray_pente,                            radius of greatest slope
   !<  - p_acv(6) = max_pente,                            greatest slope
   !<  - p_acv(7) = max_pente/min_pente                   slope anisotropy factor
   !<  - p_acv(8) = highest curvature/smallest curvature, curvature anisotropy factor
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                             :: long       !! *surface length*
   integer(kind=I4), intent(in )                             :: larg       !! *surface width*
   logical(kind=I4), intent(in )                             :: omp        !! *multithreaded ?*
   real   (kind=R8), intent(in ), optional                   :: cut        !! *cut height*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg)  :: tabin      !! *surface acf array*
   real   (kind=R8), intent(in ), dimension(1:2)             :: scale_xy   !! *lag along x and y in micrometers*
   real   (kind=R8), intent(out), dimension(1:8)             :: p_acv      !! *vector containing anisotropy outputs*

      integer(kind=I4) :: i, k, ll, p, q, qp, qm, nb_p, nx, ny, lo2, la2, funit, nb_th
      logical(kind=I4) :: verif
      real   (kind=R8) :: x, y, xb, yb, xm, ym, xp, yp, inc_a
      real   (kind=R8) :: r, r1, rc, theta, h0, hh, h1, h2, h3, h4, coupe, pente_locale, pente_locale_precede
      real   (kind=R8) :: angle, axe_a, axe_b, ech_x, ech_y, ech_z, ech_r
      real   (kind=R8) :: min_pente, max_pente, ray_pente
      real   (kind=R8) :: c, s

      integer(kind=I4), dimension(0:179) :: p_min, e_angle
      integer(kind=I4), dimension(1:2)   :: loc_max
      real   (kind=R8), dimension(0:179) :: pente_max
      real   (kind=R8), dimension(0:359) :: ellipse

      real   (kind=R8), allocatable, dimension(:,:) :: tabou, tab_tmp
      real   (kind=R8), allocatable, dimension(:)   :: courbure

      ech_x = scale_xy(1)  !SCALE_IMG%dx * unit2IUf(SCALE_IMG%dx_unit) / 1.e-6  ! x lag in micron
      ech_y = scale_xy(2)  !SCALE_IMG%dy * unit2IUf(SCALE_IMG%dy_unit) / 1.e-6  ! y lag in micron
      ech_z = 1.                                                                ! acf(0,0) = 1

      loc_max(1:2) = maxloc( tabin(1:long, 1:larg) )
      lo2 = loc_max(1)
      la2 = loc_max(2)

      !  on prend une surface carrée inscrite
      ll = min(lo2, la2) - 1

      allocate(    tabou(0:ll, 0:359) )
      allocate( courbure(      0:359) )
      tabou    = EPS_R8
      courbure = EPS_R8

      allocate( tab_tmp(1:long, 1:larg) )

      tab_tmp(1:long, 1:larg) = tabin(1:long, 1:larg)

      verif = .false.

      if ( verif ) then ! sortie xyz
         call get_unit(funit)
         open(funit,file='out/test_pol.dat')
      endif

      !  angle increment for polar representation of the rosette
      inc_a = 2*PI_R8/360

      !  determination of heights on a rosette diameter, obtained by linear interpolation :
      !  a point “falls” into a rectangular element [h1,h2,h3,h4], its height is determined by lin. interp.
      nb_p = 0

      nb_th = 1
      if (omp) nb_th = omp_get_num_procs()

      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nb_th) IF (omp)
      !$OMP DO SCHEDULE (STATIC,(ll+1)/nb_th) PRIVATE(p,r,q,theta,x,y,xm,ym,h1,h2,h3,h4,hh)

      do p = 0, ll                !  identifying a point on the diameter
         r = p                    !  corresponding algebraic radius
         do q = 0, 359            !  angular increment identification
            theta = q*inc_a

            !  projection on x and y of the point marked by its radius and angle, taking the lower integer,
            !  gives the number of the bottom line and left-hand column of the rectangle
            !  the remainder (x-nx) represents the abscissa of the point in the rectangle with sides=1
            !  the 0.9999 coefficient is used to avoid falling right on an existing point
            x = lo2 + r * cos(theta) * 0.9999_R8 ; nx = floor(x) ; xb = x -nx
            y = la2 + r * sin(theta) * 0.9999_R8 ; ny = floor(y) ; yb = y -ny

            xm = UN -xb ; xp = xb
            ym = UN -yb ; yp = yb

            if ( nx+1 <= long .and. ny+1 <= larg .and. &
                 nx   >= 1    .and. ny   >= 1) then
               ! attention r may be greater than lo2 or la2
               h1 = tab_tmp(nx   , ny   )
               h2 = tab_tmp(nx +1, ny   )
               h3 = tab_tmp(nx +1, ny +1)
               h4 = tab_tmp(nx   , ny +1)

               hh = h1*xm*ym + h2*xp*ym + h3*xp*yp + h4*xm*yp
               tabou(p, q) = hh
               nb_p = nb_p +1
            else
               hh = 0.
            endif
            if ( verif ) write(funit,*) real(x,kind=R4), real(y,kind=R4), real(hh,kind=R4)
         enddo
      enddo

      !$OMP END DO
      !$OMP END PARALLEL

      deallocate( tab_tmp )

      if ( verif ) close(funit)

      if ( present(cut) ) then
         coupe = cut
      else
         coupe  = 0.5
      endif

      do q = 0, 359        !  angle en deg

         theta = q*inc_a   !  angle en rad
         ech_r = sqrt( (ech_x*cos(theta))**2 + &   !
                       (ech_y*sin(theta))**2 )     !  unit according to angle q
         ellipse(q)  = (ll -1)*ech_r               !  max value of ellipse radius

         do p = 0, ll -1                           !  identifying a point on the diameter

            r1 = p                                 !  algebraic radius
            h1 = tabou( p   , q )
            h2 = tabou( p +1, q )
            if (abs(h2) < 10*EPS_R8) exit          !  useful for ll = floor(sqrt(UN*(lo2**2 +la2**2))) -1

            if (h1 > coupe .and. h2 < coupe) then
               rc = r1 +(h1 -coupe)/(h1 -h2)
               ellipse(q)  = rc*ech_r
               exit  ! if you don't pass here: no intersection with the cutting plane -> max value taken by default
            endif

         enddo

         ! curvature averaged over 3 points
         rc = 0
         do p = 1, 3

            h0 = tabou( p - 1, q )
            h1 = tabou( p    , q )
            h2 = tabou( p + 1, q )
            if (abs(h2) < 10*EPS_R8) exit          !  useful for ll = floor(sqrt(UN*(lo2**2 +la2**2))) -1

            rc = rc + ( h0 - 2 * h1 + h2 )

         enddo
         courbure(q) = - ( rc / ( 2 * ech_r**2 ) ) / 3.

      enddo

      do q = 0, 179
         e_angle(q) = 2 * q ! angle doubled for the continuity of sin and cos
      enddo

      call sort_array2( tab_inout = ellipse(0:179),                   &  !
                             tab1 = e_angle(0:179),                   &  !
                                n = 180 )                                !

      axe_b   = sum(ellipse(  0:  2))/3.
      axe_a   = sum(ellipse(177:179))/3.

      c = sum( cos(e_angle(175:179) * PI_R8 / 180) ) / 5.
      s = sum( sin(e_angle(175:179) * PI_R8 / 180) ) / 5.

      angle   = atan2(s, c) * 180 / PI_R8 / 2     ! angle halfed (see above)

      p_acv(1) = axe_a
      p_acv(2) = axe_b
      p_acv(3) = axe_a/axe_b
      p_acv(4) = nint(angle)

      !-----------------------------------------------

      p_min(0:179) = ll -1
      do q = 0, 179                                !  incrément angulaire

         qp = q +1 ; if (qp>179) qp = 0
         qm = q -1 ; if (qm<0  ) qm = 179
         theta = q*inc_a
         ech_r = sqrt( (ech_x*cos(theta))**2 +  &  !
                       (ech_y*sin(theta))**2 )     !  unit according to angle q

         pente_locale_precede = 0._R8
         do p = 1, ll-1

            pente_locale = ( (tabou(p+1, q ) -tabou(p-1, q )) /  2              )**2 +   &  !
                           ( (tabou(p  , qp) -tabou(p  , qm)) / (2 * p * inc_a) )**2        !

            pente_locale = sqrt( pente_locale ) * ech_z / ech_r

            if ( abs( pente_locale ) < abs( pente_locale_precede ) ) then
               p_min(q) = p                                             !  the steepest slope distance is recorded
               exit
            else
               pente_locale_precede = pente_locale
            endif

         enddo

      enddo

      call sort_array2( tab_inout = p_min(0:179),  &  !
                        n         = 180 )             !

      k = int( sum( p_min(0:4) ) / 5., kind = I4 )          !  the minimum distance

      do q = 0, 179                                         !  angular increment

         qp = q +1 ; if (qp > 179) qp = 0
         qm = q -1 ; if (qm < 0  ) qm = 179
         theta = q * inc_a

         ech_r = sqrt( ( ech_x * cos(theta) )**2 + &  !
                        (ech_y * sin(theta) )**2 )    !  unit according to angle q

         pente_locale = ( (tabou(k+1, q ) - tabou(k-1, q )) /  2              )**2 +   &  !
                        ( (tabou(k  , qp) - tabou(k  , qm)) / (2 * k * inc_a) )**2        !

         pente_locale = sqrt( pente_locale ) * ech_z / ech_r
         pente_max(q) =  abs( pente_locale )                               !  for angle q, greater slope

      enddo

      do q = 0, 179
         e_angle(q) = q
      enddo

      call sort_array2( tab_inout = pente_max(0:179),                   &  !
                             tab1 =   e_angle(0:179),                   &  !
                                n = 180 )                                  !

      angle = sum( e_angle(175:179) ) / 5.

      theta = angle*inc_a
      ech_r = sqrt( ( ech_x * cos(theta) )**2 +    &  !
                    ( ech_y * sin(theta) )**2 )       !  unit according to angle q

      min_pente = sum( pente_max(  0:  2) ) / 3.
      max_pente = sum( pente_max(177:179) ) / 3.
      ray_pente = k * ech_r

      p_acv(5) = ray_pente
      p_acv(6) = max_pente
      p_acv(7) = max_pente / min_pente         !  anisotropy indicator
      p_acv(8) = maxval( courbure(0:359) ) / minval( courbure(0:359) )

      deallocate( tabou, courbure )

   return
   endsubroutine ellipse_acf


   subroutine acv(tab_in, tab_out, long, larg, sub_samp)
   !================================================================================================
   !< @note Function that returns the *acf* of an array.
   !<
   !< \[
   !< \begin{align*}
   !<    acf(i,j) &= (z \ast z)(i,j) = \sum_{k,l}^{n,n} z(k+1-i,l+1-j)z(k,l)  \\
   !<    TF(acf)  &= ACF = Z \cdot Z                                          \\
   !<    acf      &= TF^{-1}(ACF) = TF^{-1}(Z^2)
   !< \end{align*}
   !< \]
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long        !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg        !! *2D array width*
   logical(kind=I4), intent(in )                            :: sub_samp    !! *sampling?*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in      !! *input array*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_out     !! *acf of the input array*

      integer(kind=I4) :: nx2, ny2, iex, iey, ibx, iby, i, j, lo2, la2
      real   (kind=R8) :: tmp

      integer(kind=I4), dimension(1:2) :: loc_max

      complex(kind=R8), dimension(:,:), allocatable :: cmpl1, cmpl2
      real   (kind=R8), dimension(:,:), allocatable :: tab_ext1, tab_ext2

      type(SCALE_SURF)  :: scal_surf

      type(MOMENT_STAT) :: m_res

      ! 0-padding
      if ( PAD_FFT_ANI < 0 ) then

         nx2 = long
         ny2 = larg

         ibx = 1 ; iex = long
         iby = 1 ; iey = larg

      else

         nx2 = 2 * ( nint(PAD_FFT_ANI * long) / 2 )
         ny2 = 2 * ( nint(PAD_FFT_ANI * larg) / 2 )

         ibx = max( ceiling( (nx2 - long)/2. ), 1 ) ; iex = ibx + long - 1
         iby = max( ceiling( (ny2 - larg)/2. ), 1 ) ; iey = iby + larg - 1

      endif

      allocate( tab_ext1(1:nx2, 1:ny2),   &  !
                tab_ext2(1:nx2, 1:ny2) )     !

      allocate(    cmpl1(1:nx2, 1:ny2),   &  !
                   cmpl2(1:nx2, 1:ny2) )     !

      tab_ext1(1:nx2, 1:ny2) = 0

      call calc_moments(   tab = tab_in(1:long, 1:larg),    &  !
                            mx = m_res,                     &  !
                        nb_mom = 2 )                           !

      tab_ext1(ibx:iex, iby:iey) = ( tab_in(1:long, 1:larg) - m_res%mu ) / m_res%si

      if ( APO_FFT_ANI /= "no_apo" ) then

         call apod( tab_in = tab_ext1(1:nx2, 1:ny2),  &  !
                   tab_out = tab_ext2(1:nx2, 1:ny2),  &  !
                      long = nx2,                     &  !
                      larg = ny2,                     &  !
                  type_apo = trim(APO_FFT_ANI),       &  !
                     param = 0.1_R8 )                    !
      else

         tab_ext2(1:nx2, 1:ny2) = tab_ext1(1:nx2, 1:ny2)

      endif

      !----------------

      cmpl1(1:nx2, 1:ny2) = cmplx( tab_ext2(1:nx2, 1:ny2), 0, kind = R8 )

      !----------------

      if (sub_samp) then

         call tab_calc_fftw3(   sens = FORWARD,                  &  !
                              tab_in = cmpl1(1:nx2, 1:ny2),      &  !
                              tab_ou = cmpl2(1:nx2, 1:ny2),      &  !
                                long = nx2,                      &  !
                                larg = ny2)                         !

      else

         call calc_fftw3(   sens = FORWARD,                    &  !
                          tab_in = cmpl1(1:nx2, 1:ny2),        &  !
                          tab_ou = cmpl2(1:nx2, 1:ny2),        &  !
                            long = nx2,                        &  !
                            larg = ny2)                           !

      endif

      cmpl1(1:nx2, 1:ny2) = cmplx( abs( cmpl2(1:nx2, 1:ny2) )**2, 0, kind = R8 )

      ! théorème de wiener

      if (sub_samp) then

         call tab_calc_fftw3(   sens = BACKWARD,               &  !
                              tab_in = cmpl1(1:nx2, 1:ny2),    &  !
                              tab_ou = cmpl2(1:nx2, 1:ny2),    &  !
                                long = nx2,                    &  !
                                larg = ny2)                       !

      else

         call calc_fftw3(   sens = BACKWARD,                   &  !
                          tab_in = cmpl1(1:nx2, 1:ny2),        &  !
                          tab_ou = cmpl2(1:nx2, 1:ny2),        &  !
                            long = nx2,                        &  !
                            larg = ny2)                           !

      endif

      tab_ext1(1:nx2, 1:ny2) = real(cmpl2(1:nx2, 1:ny2), kind=R8)

      call trans_corner2center(  tab_in  = tab_ext1(1:nx2, 1:ny2),  &  !
                                 tab_out = tab_ext2(1:nx2, 1:ny2),  &  !
                                 long    = nx2,                     &  !
                                 larg    = ny2  )                      !

      tab_out(1:long, 1:larg) = tab_ext2(ibx:iex, iby:iey)

      ! normalisation
      loc_max(1:2) = maxloc( tab_out(1:long, 1:larg) )
      lo2 = loc_max(1)
      la2 = loc_max(2)

      tmp = tab_out(lo2, la2)

      tab_out(1:long, 1:larg) = tab_out(1:long, 1:larg) / tmp

      if (.false.) then
         call init_scal( scal   = scal_surf,    &  !
                         nx     = long,         &  !
                         ny     = larg,         &  !
                         lx     = 1._R8,        &  !
                         ly     = 1._R8,        &  !
                         unit_z = 'm ')            !

         call write_surf( nom_fic = "test_acv.sur",            &  !
                          tab_s   = tab_out(1:long, 1:larg),   &  !
                          scal    = scal_surf )                   !
      endif

      deallocate(cmpl1, cmpl2)
      deallocate(tab_ext1, tab_ext2)

   return
   endsubroutine acv


   subroutine correlation_parameters(tab, long, larg, res, cut, sub_plane, sub_sampl, scale_xy, omp)
   !================================================================================================
   !< @note Function that returns [[ellipse_acf]] parameters calculated on the autocorrelation
   !< function. But prior to the acf calculation, the mean plane is subtracted.
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                             :: long       !! *2D array length*
   integer(kind=I4), intent(in )                             :: larg       !! *2D array height*
   logical(kind=I4), intent(in )                             :: sub_plane  !! *subtract least square plane?*
   logical(kind=I4), intent(in )                             :: sub_sampl  !! *subsampling?*
   logical(kind=I4), intent(in )                             :: omp        !! *multithreaded ?*
   real   (kind=R8), intent(in ), optional                   :: cut        !! *cut height*
   real   (kind=R8), intent(in ), dimension(1:2)             :: scale_xy   !! *lag along x and y in micrometers*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg)  :: tab        !! *2D array in*
   real   (kind=R8), intent(out), dimension(1:8)             :: res        !! *correlation parameters*

      real(kind=R8), dimension(1:long, 1:larg) :: tab_tmp1, tab_tmp2

      if ( sub_plane ) then

         ! mean plane subtracted
         call least_squares_tcheby(tab_in =      tab(1:long, 1:larg),  &  ! IN
                                  tab_out = tab_tmp1(1:long, 1:larg),  &  ! OUT
                                    long1 = long,                      &  ! IN
                                    long2 = larg,                      &  ! IN
                                    nvarx = 1,                         &  ! IN
                                    nvary = 1)                            ! IN

         tab_tmp1(1:long, 1:larg) = tab(1:long, 1:larg) - tab_tmp1(1:long, 1:larg)

      else

         tab_tmp1(1:long, 1:larg) = tab(1:long, 1:larg)

      endif

      call acv( tab_in    = tab_tmp1(1:long, 1:larg),  &  ! IN
                tab_out   = tab_tmp2(1:long, 1:larg),  &  ! OUT
                long      = long,                      &  ! IN
                larg      = larg,                      &  ! IN
                sub_samp  = sub_sampl )                   ! IN

      call ellipse_acf( tabin = tab_tmp2(1:long, 1:larg),   &  ! IN
                         long = long,                       &  ! IN
                         larg = larg,                       &  ! IN
                        p_acv = res(1:8),                   &  ! OUT
                          cut = cut,                        &  ! IN
                     scale_xy = scale_xy,                   &  ! IN
                          omp = omp )                          ! IN

   return
   endsubroutine correlation_parameters

endmodule anisotropy
