!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: april, 06 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Surface smoothers**
!<  </span>

module filter
use data_arch,     only: I4, R8, HIG_R8, PI_R8
use miscellaneous, only: trans_center2corner
use stat_mom,      only: calc_median, moment_stat, calc_moments
use sort_arrays,   only: sort_array2
use fftw3,         only: calc_fftw3, tab_calc_fftw3, extend, &  !
                         PAD_FFT, FORWARD, BACKWARD
!$ use omp_lib
implicit none

real(kind=R8) :: PAD_FFT_FILTER = PAD_FFT
! with 1.5, the results for gaussian filtering are the same as Mountains'

private

public :: median_filter, median_smooth, morpho_filter, soften, fft_filter, PAD_FFT_FILTER

!> {!FILTR/src/inc_doc/comp_filtr.md!}
!> {!css/button.html!}

contains

   subroutine morpho_filter(tabin, tabou, long, larg, scale_xyz, ray, omp, nb_div, mtype)
   !================================================================================================
   !< @note
   !<
   !< Morphological filter: uses combinations of [[roll_smooth]] to provide all kind of transformation :
   !<
   !< + closing
   !< + opening
   !< + dilation
   !< + erosion
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long        !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg        !! *2D array width*
   integer(kind=I4), intent(in )                            :: nb_div      !! *number of macro elements along an axis*
   logical(kind=I4), intent(in )                            :: omp         !! *if multithreading*
   real   (kind=R8), intent(in )                            :: ray         !! *roll radius*
   character(len=*), intent(in )                            :: mtype       !! *closing, opening, dilation or erosion*
   real   (kind=R8), intent(in ), dimension(1:3)            :: scale_xyz   !! *lag along x, y and scale z*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tabin       !! *2D array in*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tabou       !! *2D array out*

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      logical(kind=I4) :: op1, op2
      integer(kind=I4) :: k

      select case ( mtype(1:7) )

         case("closing")

            op1 = .true. ; op2 = .true. ; k = +1

         case("opening")

            op1 = .true. ; op2 = .true. ; k = -1

         case("dilatio")

            op1 = .true. ; op2 = .false. ; k = +1

         case("erosion")

            op1 = .true. ; op2 = .false. ; k = -1

         case default

            stop "bad choice, morpho_filter"

      endselect

      if (op1) then

         call roll_smooth( tabin       = tabin(1:long, 1:larg),      &  ! IN
                           tabou       = tabou(1:long, 1:larg),      &  ! OUT
                           long        = long,                       &  ! IN
                           larg        = larg,                       &  ! IN
                           scale_xyz   = scale_xyz,                  &  ! IN
                           sgn         = +k,                         &  ! IN
                           ray         = ray,                        &  ! IN
                           omp         = omp,                        &  ! IN
                           nb_div      = nb_div )                       ! IN

      endif

      if (op2) then

         allocate( tab_tmp(1:long, 1:larg) )

         call roll_smooth( tabin       = tabou(1:long, 1:larg),      &  ! IN
                           tabou       = tab_tmp(1:long, 1:larg),    &  ! OUT
                           long        = long,                       &  ! IN
                           larg        = larg,                       &  ! IN
                           scale_xyz   = scale_xyz,                  &  ! IN
                           sgn         = -k,                         &  ! IN
                           ray         = ray,                        &  ! IN
                           omp         = omp,                        &  ! IN
                           nb_div      = nb_div )                       ! IN

         tabou(1:long, 1:larg) = tab_tmp(1:long, 1:larg)

         deallocate( tab_tmp )

      endif

   return
   endsubroutine morpho_filter


   subroutine roll_smooth(tabin, tabou, long, larg, scale_xyz, sgn, ray, omp, nb_div)
   !================================================================================================
   !< @note
   !<
   !< A ball of radius "ray" rolls on / below the surface, hence defining a closing or an opening enveloppe.
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long        !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg        !! *2D array width*
   integer(kind=I4), intent(in )                            :: sgn         !! *+ 1: dilation, -1:erosion*
   integer(kind=I4), intent(in )                            :: nb_div      !! *number of macro elements along an axis*
   logical(kind=I4), intent(in )                            :: omp         !! *if multithreading*
   real   (kind=R8), intent(in )                            :: ray         !! *roll radius*
   real   (kind=R8), intent(in ), dimension(1:3)            :: scale_xyz   !! *lag along x, y and scale z*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tabin       !! *2D array in*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tabou       !! *2D array out*

      integer(kind=I4) :: i, j, ii, jj, nb_th
      integer(kind=I4) :: hw

      integer(kind=I4) :: ik, jk, idiv, jdiv, ista, iend, jsta, jend

      real(kind=R8) :: h1, h2, ht, tmp, delta_h_max
      real(kind=R8) :: ech_x, ech_y

      real(kind=R8), allocatable, dimension(:,:) :: elem, tab_tmp

      integer(kind=I4), allocatable, dimension(:,:) :: thw

      ech_x = scale_xyz(1) !SCALE_IMG%dx * unit2IUf(SCALE_IMG%dx_unit)
      ech_y = scale_xyz(2) !SCALE_IMG%dy * unit2IUf(SCALE_IMG%dy_unit)

      h1 = minval( tabin(1:long, 1:larg) )
      h2 = maxval( tabin(1:long, 1:larg) )

      delta_h_max = h2 - h1

      ht = min( -abs(h1), -abs(h2) )
      h2 = max( +abs(h1), +abs(h2) )
      h1 = ht

      ! the normal width of the ball is : hw = int( ray / ech_x )
      ! However the ball curvature makes the ball height sometimes higher than the surfaces heights
      hw = int( sqrt( 2*delta_h_max*ray - delta_h_max**2 ) / ech_x ) !+ 1

      allocate( elem(   -hw    :hw,        -hw    :hw       ) )
      allocate( tab_tmp(-hw + 1:hw + long, -hw + 1:hw + larg) ) ! surface is extended

      tab_tmp(1:long, 1:larg) = tabin(1:long, 1:larg)           ! original surface

      do i = 1, long

         tab_tmp(i,  -hw + 1:         0) = tab_tmp(i,    1)
         tab_tmp(i, larg + 1: larg + hw) = tab_tmp(i, larg)

      enddo

      do j = -hw + 1, larg + hw

         tab_tmp( -hw + 1:        0, j) = tab_tmp(1,    j)
         tab_tmp(long + 1:long + hw, j) = tab_tmp(long, j)

      enddo

      nb_th = 1
      if (omp) then
         nb_th = omp_get_num_procs()
      endif

      allocate( thw(1:nb_div, 1:nb_div) )    ! number of macro squares on which the ball active width is determined

      idiv = long / nb_div
      jdiv = larg / nb_div

      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nb_th) IF (omp)
      !$OMP DO  SCHEDULE (STATIC, max(nb_div/nb_th, 1)) PRIVATE(jk, ik, delta_h_max, ista, iend, jsta, jend)

      do jk = 1, nb_div

      if (jk==nb_div) then ; jend = larg ; else ; jend = jk * jdiv ; endif ; jsta = 1 + jdiv * (jk-1)

      do ik = 1, nb_div

      if (ik==nb_div) then ; iend = long ; else ; iend = ik * idiv ; endif ; ista = 1 + idiv * (ik-1)

            delta_h_max = maxval( tab_tmp(-hw + ista:hw + iend, -hw + jsta:hw + jend) ) - & !
                          minval( tab_tmp(-hw + ista:hw + iend, -hw + jsta:hw + jend) )

            ! beyond the present width, the ball height is greater than the surface height
            thw(ik, jk) = int( sqrt( 2*delta_h_max*ray - delta_h_max**2 ) / ech_x )

      enddo

      enddo

      !$OMP END DO
      !$OMP END PARALLEL

      if (sgn == +1) then

         do jj = -hw, hw
         do ii = -hw, hw
            tmp = ray**2 - (ii * ech_x)**2 - (jj * ech_y)**2
            if ( tmp < 0. ) then
               elem(ii, jj) = ray               + 1.1 * abs(h2)
            else
               ! the ball location is a little above the surface
               elem(ii, jj) = ray - sqrt( tmp )
            endif
         enddo
         enddo

      else

         do jj = -hw, hw
         do ii = -hw, hw
            tmp = ray**2 - (ii * ech_x)**2 - (jj * ech_y)**2
            if ( tmp < 0. ) then
               elem(ii, jj) = -ray               - 1.1 * abs(h1)
            else
               ! the ball location is a little below the surface
               elem(ii, jj) = -ray + sqrt( tmp )
            endif
         enddo
         enddo

      endif

      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nb_th) IF (omp)
      !$OMP DO SCHEDULE (STATIC,larg/nb_th) PRIVATE(i, j, jk, ik, hw)

      do j = 1, larg ; jk = 1 + j/(jdiv + 1)
      do i = 1, long ; ik = 1 + i/(idiv + 1)

         hw = thw(ik, jk)

         tabou(i, j) = - minval( sgn * ( elem(-hw    :hw    , -hw    :hw    ) - & !
                                      tab_tmp(-hw + i:hw + i, -hw + j:hw + j) ) )

      enddo
      enddo

      !$OMP END DO
      !$OMP END PARALLEL

      if (sgn == +1) then

         tabou(1:long, 1:larg) = + tabou(1:long, 1:larg) + ray

      else

         tabou(1:long, 1:larg) = - tabou(1:long, 1:larg) - ray

      endif

      deallocate( elem, tab_tmp, thw )

   return
   endsubroutine roll_smooth


   subroutine median_smooth(tab, long, larg, kernel, omp)
   !================================================================================================
   !! Very classical smoothing
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                            :: long   !! *2D array length*
   integer(kind=I4), intent(in   )                            :: larg   !! *2D array width*
   integer(kind=I4), intent(in   )                            :: kernel !! *kernel size*
   logical(kind=I4), intent(in   )                            :: omp    !! *if multithreading*
   real   (kind=R8), intent(inout), dimension(1:long, 1:larg) :: tab    !! *2D array*

      integer(kind=I4) :: i, j, k, ii, jj, nt, nk, nb_th
      real(kind=R8)    :: md

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp, t
      real(kind=R8), allocatable, dimension(:)   :: vt

      allocate( tab_tmp(1:long, 1:larg) ) ; tab_tmp = HIG_R8

      k  = kernel
      nt = ( 2*k + 1 )*( 2*k + 1 )

      allocate( t(-k:k, -k:k), vt(1:(2*k+1)*(2*k+1)) )

      nb_th = 1
      if (omp) then
         nb_th = omp_get_num_procs()
      endif

      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nb_th) IF (omp)
      !$OMP DO SCHEDULE (STATIC,larg/nb_th) PRIVATE(i, t, nk, ii, jj, vt, md)

      do j = 1, larg
      do i = 1, long

         t(-k:k, -k:k) = -HIG_R8/10
         nk = 0

         do jj = -k, +k

            if (j +jj < 1 .or. j +jj > larg) cycle

            do ii = -k, +k
               if (i +ii < 1 .or. i +ii > long) cycle
               nk = nk +1
               t(ii, jj) = tab(i +ii, j +jj)
            enddo

         enddo

         vt(1:nt) = reshape(t(-k:k, -k:k), [nt])
         call sort_array2(tab_inout = vt(1:nt), n = nt)

         call calc_median(tab = vt(nt -nk +1:nt), md = md)

         tab_tmp(i, j) = md

      enddo
      enddo

      !$OMP END DO
      !$OMP END PARALLEL

      tab(1:long, 1:larg) = tab_tmp(1:long, 1:larg)

      deallocate( tab_tmp, t, vt )

   return
   endsubroutine median_smooth


   subroutine median_filter(tab, long, larg, snb, kernel, sig, omp)
   !================================================================================================
   !! A bit more complex filter: the overall height standard deviation is taken into account
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                            :: long   !! *2D array length*
   integer(kind=I4), intent(in   )                            :: larg   !! *2D array width*
   integer(kind=I4), intent(in   )                            :: snb    !! *patch number along a direction*
   integer(kind=I4), intent(in   )                            :: kernel !! *kernel size*
   logical(kind=I4), intent(in   )                            :: omp    !! *if multithreading*
   real   (kind=R8), intent(in   )                            :: sig    !! *error std*
   real   (kind=R8), intent(inout), dimension(1:long, 1:larg) :: tab    !! *2D array*

      integer(kind=I4) :: i, j, k
      real(kind=R8)    :: md

      integer(kind=I4), dimension(1:snb + 1) :: li, lj

      real(kind=R8), dimension(1:snb * snb)  :: ect

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp1, tab_tmp2

      type( moment_stat ) :: mx_smooth

      allocate( tab_tmp1(1:long, 1:larg) )
      allocate( tab_tmp2(1:long, 1:larg) )

      ! first determine the difference between the input surface and a median filtered one
      tab_tmp1(1:long, 1:larg) = tab(1:long, 1:larg)
      call median_smooth(  tab    = tab_tmp1(1:long, 1:larg),  &  !
                           kernel = kernel,                    &  !
                           long   = long,                      &  !
                           larg   = larg,                      &  !
                           omp    = omp )                         !
      tab_tmp2(1:long, 1:larg) = tab(1:long, 1:larg) - tab_tmp1(1:long, 1:larg)

      ! bounds when patching domain
      li(1) = 1 ; li(snb + 1) = long
      lj(1) = 1 ; lj(snb + 1) = larg

      do i = 2, snb
         li(i) = li(i - 1) + int( real(long, kind=R8)/snb, kind=I4 )
         lj(i) = lj(i - 1) + int( real(larg, kind=R8)/snb, kind=I4 )
      enddo

      k  = 0
      do j = 1, snb
      do i = 1, snb

         call calc_moments(   tab = reshape( tab_tmp2( li(i):li(i + 1) - 1, lj(j):lj(j + 1) - 1 ),    &  !
                                             [( li(i) - li(i + 1) ) * ( lj(j) - lj(j + 1) )] ),       &  !
                               mx = mx_smooth,                                                        &  !
                           nb_mom = 2 )                                                                  !
         k = k +1
         ect(k) = mx_smooth%si

      enddo
      enddo

      call calc_median( tab = ect( 1:snb*snb ), &  !
                         md = md )                 !

      call calc_moments(   tab = reshape( tab_tmp2(1:long, 1:larg),  &  !
                                          [long * larg] ),           &  !
                            mx = mx_smooth,                          &  !
                        nb_mom = 2 )                                    !

      where( abs(tab_tmp2(1:long, 1:larg) -mx_smooth%mu) > sig*md ) tab(1:long, 1:larg) = tab_tmp1(1:long, 1:larg)

      deallocate( tab_tmp1, tab_tmp2 )

   return
   endsubroutine median_filter


   subroutine soften(tabin, mask, tabout, long, larg)
   !================================================================================================
   !< @note Function to smooth out a 2D array: each point is replaced by a weighted mean of its neighbors.
   !<
   !< \[
   !<   h_{i,j} = \frac{1}{16} \left( 4 h_{i, j} + 2 h_{i + 1, j    } + 2 h_{i - 1, j    } + 2 h_{i    , j + 1} + 2 h_{i    , j - 1}
   !<                                              + h_{i + 1, j - 1} +   h_{i - 1, j + 1} +   h_{i - 1, j + 1} +   h_{i + 1, j - 1} \right)
   !< \]
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                                      :: long       !! *2D array width*
   integer(kind=I4), intent(in )                                      :: larg       !! *2D array height*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg)           :: tabin      !! *2D array in*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg)           :: tabout     !! *2D array out*
   integer(kind=I4), intent(in ), dimension(1:long, 1:larg), optional :: mask       !! *mask*

      integer(kind=I4) :: i, j

      tabout(1:long, 1:larg) = tabin(1:long, 1:larg)

      if ( present(mask) ) then

         do j = 1 +1, larg -1
         do i = 1 +1, long -1

            if ( sum(mask(i-1:i+1, j-1:j+1)) < 9 ) then
               cycle
            else
               tabout(i, j) = ( 2*tabin(i, j) +tabin(i +1, j   ) +tabin(i -1, j   ) +                                          &  !
                                               tabin(i   , j +1) +tabin(i   , j -1) + ( tabin(i +1, j -1) +tabin(i -1, j -1) + &  !
                                                                                        tabin(i -1, j +1) +tabin(i +1, j +1) ) / 2._R8 ) / 8
            endif

         enddo
         enddo

      else

         do j = 1 +1, larg -1
         do i = 1 +1, long -1
            tabout(i, j) = ( 2*tabin(i, j) +tabin(i +1, j   ) +tabin(i -1, j   ) +                                          &  !
                                            tabin(i   , j +1) +tabin(i   , j -1) + ( tabin(i +1, j -1) +tabin(i -1, j -1) + &  !
                                                                                     tabin(i -1, j +1) +tabin(i +1, j +1) ) / 2._R8 ) / 8
         enddo
         enddo

      endif

   return
   endsubroutine soften


   subroutine fft_filter(tab, long, larg, cutoff, bf_tab, multi_fft, pad, ext, type_apo, shift)
   !================================================================================================
   !! Classical Gaussian filter
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in ) :: long                                !! *2D array width*
   integer(kind=I4), intent(in ) :: larg                                !! *2D array height*
   real   (kind=R8), intent(in ) :: cutoff                              !! *cut-off wavelength*
   logical(kind=I4), intent(in ) :: multi_fft                           !! *multiple fft at once ?*
   real   (kind=R8), intent(in ), optional :: pad                       !! *fft padding*
   character(len=*), intent(in ), optional :: ext                       !! *extension*
   character(len=*), intent(in ), optional :: type_apo                  !! *apodization type*
   real   (kind=R8), intent(in ), dimension(2), optional    :: shift    !! *surface shift fraction along x and y*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab      !! *2D array in*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: bf_tab   !! *2D array out*

      integer(kind=I4) :: i, j, nx2, ny2, iex, iey, ibx, iby
      real   (kind=R8) :: o_pad
      logical(kind=I4) :: with_pad

      character(len=:), allocatable :: o_ext
      character(len=:), allocatable :: o_type_apo

      complex(kind=R8) :: cdx, cdy

      real   (kind=R8), dimension(:,:), allocatable :: tab_ext, gauss_tab, top_tab

      complex(kind=R8), dimension(:,:), allocatable :: cmpl1, cmpl2, shift_tab

      with_pad = .true.

      if ( .not.present( ext ) ) then

         o_ext = 'symmetry'

      else

         o_ext = ext

      endif

      if ( .not.present( type_apo ) ) then

         o_type_apo = 'no_apo'

      else

         o_type_apo = type_apo

      endif

      if ( .not.present( pad ) ) then

         ! default padding
         o_pad = PAD_FFT_FILTER

      else

         ! no padding
         if ( pad < 0. ) then

            o_pad = 1

            nx2 = long
            ny2 = larg

            with_pad = .false.

         else

            ! padding with argument
            o_pad = pad

         endif

      endif

      if ( with_pad ) then

         ! make even surface dimensions

         nx2 = 2 * ( nint(o_pad * long)/2 )
         ny2 = 2 * ( nint(o_pad * larg)/2 )

      endif

      allocate( shift_tab(1:nx2, 1:ny2) )

      shift_tab(1:nx2, 1:ny2) = 1

      if ( present( shift ) ) then

         ! shift surface with frequencies rotation
         call shifting(nx2, ny2, shift, shift_tab)

      endif

      allocate( tab_ext(1:nx2, 1:ny2)  )     !

      allocate( cmpl1(1:nx2, 1:ny2),      &  !
                cmpl2(1:nx2, 1:ny2) )        !

      if ( nx2 > long ) then

         ibx = ceiling( (nx2 - long)/2. ) ; iex = ibx + long - 1
         iby = ceiling( (ny2 - larg)/2. ) ; iey = iby + larg - 1

         call extend(   tab_in = tab(1:long, 1:larg),       &  !
                       tab_out = tab_ext(1:nx2, 1:ny2),     &  !
                            nx = long,                      &  !
                            ny = larg,                      &  !
                           nx2 = nx2,                       &  !
                           ny2 = ny2,                       &  !
                           ext = o_ext,                     &  !
                      type_apo = o_type_apo )                  !

      else

         tab_ext(1:nx2, 1:ny2) = tab(1:long, 1:larg)

      endif

      cmpl1(1:nx2, 1:ny2) = cmplx( tab_ext(1:nx2, 1:ny2), 0, kind = R8 )

      if (multi_fft) then

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

      if ( cutoff < 0. ) then

         ! TOP HAT FILTER
         allocate( top_tab(1:nx2, 1:ny2) )

         call top_hat_filter(  long     = nx2,                     &  !
                               larg     = ny2,                     &  !
                               xc       = -cutoff,                 &  !
                               top_filt = top_tab(1:nx2, 1:ny2) )     !

         cmpl1(1:nx2, 1:ny2) = cmpl2(1:nx2, 1:ny2) * top_tab(1:nx2, 1:ny2) * shift_tab(1:nx2, 1:ny2)

         deallocate( top_tab )

      else

         ! GAUSSIAN FILTER
         allocate( gauss_tab(1:nx2, 1:ny2) )

         call gaussian_filter( long       = nx2,                     &  !
                               larg       = ny2,                     &  !
                               xc         = cutoff,                  &  !
                               gauss_filt = gauss_tab(1:nx2, 1:ny2) )   !

         cmpl1(1:nx2, 1:ny2) = cmpl2(1:nx2, 1:ny2) * gauss_tab(1:nx2, 1:ny2)

         deallocate( gauss_tab )

      endif

      if (multi_fft) then

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

      tab_ext(1:nx2, 1:ny2) = real(cmpl2(1:nx2, 1:ny2), kind=R8)

      if ( nx2 > long ) then

         bf_tab(1:long, 1:larg) = tab_ext(ibx:iex, iby:iey)

      else

         bf_tab(1:long, 1:larg) = tab_ext(1:nx2, 1:ny2)

      endif

      deallocate(cmpl1, cmpl2, shift_tab)
      deallocate(tab_ext)

   return
   endsubroutine fft_filter


   subroutine shifting(long, larg, shift, shift_tab)
   !================================================================================================
   !! Rotation matrix for frequencies -> results in surface shifting
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long        !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg        !! *2D array width*
   real   (kind=R8), intent(in ), dimension(2), optional    :: shift       !! *surface shift fraction along x and y*
   complex(kind=R8), intent(out), dimension(1:long, 1:larg) :: shift_tab   !! *2D array out*

      integer(kind=I4) :: i, j
      real   (kind=R8) :: x, y

      complex(kind=R8) :: eixi, tmp

      x = - shift(1) * ( 2 * pi_r8 / (long - 1) )
      y = - shift(2) * ( 2 * pi_r8 / (larg - 1) )

      do i = 1, long

         eixi = cmplx( cos((i - 1) * x), sin((i - 1) * x), kind=r8 )

         do j = 1, larg

            shift_tab(i, j) = eixi * cmplx( cos((j - 1) * y), sin((j - 1) * y), kind=r8 )

         enddo

      enddo

   return
   endsubroutine shifting


   subroutine gaussian_filter(long, larg, xc, gauss_filt)
   !================================================================================================
   !! Gaussian kernel
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long        !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg        !! *2D array width*
   real   (kind=R8), intent(in )                            :: xc          !! *the cut-off wavelength*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: gauss_filt  !! *2D array out*

      integer(kind=I4) :: i, j
      real   (kind=R8) :: tmp, xi, xj

      real   (kind=R8), parameter :: const = sqrt( log(2._R8)/PI_R8 )

      do j = 2, larg/2 +1
      do i = 2, long/2 +1
         xi = (i-1) ; xj = (j-1)
         xi = xi/(long -1) ; xj = xj/(larg -1)
         tmp = gaussian_function(xi, xj, xc)
         gauss_filt(       +i,        +j) = tmp
         gauss_filt(long+2 -i,        +j) = tmp
         gauss_filt(       +i, larg+2 -j) = tmp
         gauss_filt(long+2 -i, larg+2 -j) = tmp
      enddo
      enddo
      do j = 2, larg/2 +1
         i = 1
         xi = (i-1) ; xj = (j-1)
         xi = xi/(long -1) ; xj = xj/(larg -1)
         tmp = gaussian_function(xi, xj, xc)
         gauss_filt(i,         j) = tmp
         gauss_filt(i, larg+2 -j) = tmp
      enddo
      do i = 2, long/2 +1
         j = 1
         xi = (i-1) ; xj = (j-1)
         xi = xi/(long -1) ; xj = xj/(larg -1)
         tmp = gaussian_function(xi, xj, xc)
         gauss_filt(        i, j) = tmp
         gauss_filt(long+2 -i, j) = tmp
      enddo
      i = 1
      j = 1
      xi = (i-1) ; xj = (j-1)
      xi = xi/(long -1) ; xj = xj/(larg -1)
      gauss_filt(i, j) = gaussian_function(xi, xj, xc)

   contains
      !-----------------------------------------
      real(kind=R8) function gaussian_function(xi, xj, xc)
      implicit none
      real(kind=R8), intent(in) :: xi
      real(kind=R8), intent(in) :: xj
      real(kind=R8), intent(in) :: xc ! fréquence de coupure, plus exactement proportion : (freq coup) / (nb points)

         gaussian_function = exp( -PI_R8 * const**2 * (xi /  xc                           ) **2 ) *   &  ! const = sqrt(ln(2)/pi)=0.47
                             exp( -PI_R8 * const**2 * (xj / (xc * (long - 1) / (larg - 1))) **2 )        !

      return
      endfunction gaussian_function
      !-----------------------------------------
   endsubroutine gaussian_filter


   subroutine top_hat_filter(long, larg, xc, top_filt)
   !================================================================================================
   !! Top-hat kernel
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long        !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg        !! *2D array width*
   real   (kind=R8), intent(in )                            :: xc          !! *the cut-off wavelength*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: top_filt    !! *2D array out*

      integer(kind=I4) :: i, j
      real   (kind=R8) :: tmp, xi, xj

      real   (kind=R8), parameter :: const = sqrt( log(2._R8)/PI_R8 )

      do j = 2, larg/2 +1
      do i = 2, long/2 +1
         xi = (i-1) ; xj = (j-1)
         xi = xi/(long -1) ; xj = xj/(larg -1)
         tmp = top_hat_function(xi, xj, xc)
         top_filt(       +i,        +j) = tmp
         top_filt(long+2 -i,        +j) = tmp
         top_filt(       +i, larg+2 -j) = tmp
         top_filt(long+2 -i, larg+2 -j) = tmp
      enddo
      enddo
      do j = 2, larg/2 +1
         i = 1
         xi = (i-1) ; xj = (j-1)
         xi = xi/(long -1) ; xj = xj/(larg -1)
         tmp = top_hat_function(xi, xj, xc)
         top_filt(i,         j) = tmp
         top_filt(i, larg+2 -j) = tmp
      enddo
      do i = 2, long/2 +1
         j = 1
         xi = (i-1) ; xj = (j-1)
         xi = xi/(long -1) ; xj = xj/(larg -1)
         tmp = top_hat_function(xi, xj, xc)
         top_filt(        i, j) = tmp
         top_filt(long+2 -i, j) = tmp
      enddo
      i = 1
      j = 1
      xi = (i-1) ; xj = (j-1)
      xi = xi/(long -1) ; xj = xj/(larg -1)
      top_filt(i, j) = top_hat_function(xi, xj, xc)

   contains
      !-----------------------------------------
      real(kind=R8) function top_hat_function(xi, xj, xc)
      implicit none
      real(kind=R8), intent(in) :: xi
      real(kind=R8), intent(in) :: xj
      real(kind=R8), intent(in) :: xc ! fréquence de coupure, plus exactement proportion : (freq coup) / (nb points)

         real(kind=R8) :: val

         val = (xi / xc) **2 + (xj / (xc * (long - 1) / (larg - 1))) **2

         top_hat_function = 0
         if ( val < 1. ) top_hat_function = 1

      return
      endfunction top_hat_function
      !-----------------------------------------
   endsubroutine top_hat_filter



endmodule filter
