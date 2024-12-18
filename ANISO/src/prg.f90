!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: november, 01 2024
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<     **Anisotropy detection. Example of use**
!<  </span>
program test_anisotropy
use data_arch,     only : I4, R8, PI_R8
use miscellaneous, only : get_unit
use surfile,       only : init_scal, read_surf, write_surf, SCALE_SURF, unit2IUf
use fftw3
use anisotropy
!$ use omp_lib
implicit none

integer(kind=I4), parameter :: nx = 800, ny = 400
real   (kind=R8), allocatable, dimension(:,:)  :: array, array_tmp

real(kind=R8), dimension(1:3)  :: ell_param
real(kind=R8), dimension(1:2)  :: ech_param
real(kind=R8), dimension(1:8)  :: param_acv
real(kind=R8), dimension(1:9)  :: param_ani
real(kind=R8), dimension(1:2)  :: scale_xy

real(kind=R8) :: lx, dx, dy

type(SCALE_SURF) :: scal_surf

integer(kind=I4) :: nnx, nny, n_th, nx2, ny2

!-----------------------------------------------------------------------------------------------------------------

call read_surf( nom_fic = "sur/test1.sur",         &  !
                tab_s   = array,                   &  !
                scal    = scal_surf )                 !

nnx = scal_surf%xres
nny = scal_surf%yres

scale_xy = [ scal_surf%dx * unit2IUf(scal_surf%dx_unit),   &  !
             scal_surf%dy * unit2IUf(scal_surf%dy_unit) ]     !

n_th = omp_get_max_threads()
call fftw_plan_with_nthreads(nthreads = n_th)

nx2 = 2 * ( nint(PAD_FFT * nnx)/2 )
ny2 = 2 * ( nint(PAD_FFT * nny)/2 )

call init_fftw3(long = nx2,  & !
                larg = ny2 )   !


call multiple_anisotropy( tabin     = array(1:nnx,1:nny),    &  !
                          long      = nnx,                   &  !
                          larg      = nny,                   &  !
                          scale_xy  = scale_xy,              &  !
                          multi_fft = .false.,               &  !
                          vec_ani   = param_ani(1:9) )          !

call end_fftw3()

write(*, '(e12.4, T20, a)' ) param_ani(1) , 'bmp, maximum over [0,179°] of the peaks mean width'
write(*, '(e12.4, T20, a)' ) param_ani(2) , 'smp, minimum over [0,179°] of the peaks mean width'
write(*, '(e12.4, T20, a)' ) param_ani(3) , 'rmp, ratio bmp/smp'

write(*, '(e12.4, T20, a)' ) param_ani(4) , 'bml, maximum over [0,179°] of the path length'
write(*, '(e12.4, T20, a)' ) param_ani(5) , 'sml, minimum over [0,179°] of the path length'
write(*, '(e12.4, T20, a)' ) param_ani(6) , 'rml, ratio bml/sml'

write(*, '(e12.4, T20, a)' ) param_ani(7) , 'bms, maximum over [0,179°] of the standard deviation of slope'
write(*, '(e12.4, T20, a)' ) param_ani(8) , 'sms, minimum over [0,179°] of the standard deviation of slope'
write(*, '(e12.4, T20, a)' ) param_ani(9) , 'rms, ratio bms/sms'

deallocate( array )

!-----------------------------------------------------------------------------------------------------------------

ell_param = [ 40.e-6_R8, 10.e-6_R8, 25._R8 ]    ! big semi-axis, small semi-axis, angle (°)

ech_param = [ 0.25e-6_R8, 0.25e-6_R8 ]          ! scale x (m/pix), scale_y (m/pix)
                                                ! domain is (800 * 0.25e-6) X (400 * 0.25e-6) = 200 µm X 100 µm

call fake_acv( acv_array = array,               &  !
               long      = nx,                  &  !
               larg      = ny,                  &  !
               param     = ell_param(1:3),      &  !
               ech       = ech_param(1:2) )        !

call init_scal( scal   = scal_surf,             &  !
                nx     = nx,                    &  !
                ny     = ny,                    &  !
                lx     = nx * ech_param(1),     &  !
                ly     = ny * ech_param(2),     &  !
                unit_z = 'm ')                     !

call write_surf( nom_fic = "out/fake_acv.sur",  &  !
                 tab_s   = array(1:nx, 1:ny),   &  !
                 scal    = scal_surf )             !

call ellipse_acf( tabin    = array(1:nx, 1:ny), &  !
                  long     = nx,                &  !
                  larg     = ny,                &  !
                  p_acv    = param_acv(1:8),    &  !
                  cut      = 0.5_R8,            &  !
                  scale_xy = ech_param(1:2),    &  !
                  omp      = .true. )              !

write(*, '(    a, T20, a)' ) '************', '*******************************************************************'
write(*, '(e12.4, T20, a)' ) param_acv(1) , 'axe_a,                                ellipsis big axis'
write(*, '(e12.4, T20, a)' ) param_acv(2) , 'axe_b,                                ellipsis small axis'
write(*, '(e12.4, T20, a)' ) param_acv(3) , 'axe_a/axe_b                           another anisotropy factor'
write(*, '(e12.4, T20, a)' ) param_acv(4) , 'nint(angle/inc_a),                    main texture orientation'
write(*, '(e12.4, T20, a)' ) param_acv(5) , 'ray_pente,                            radius of greatest slope'
write(*, '(e12.4, T20, a)' ) param_acv(6) , 'max_pente,                            greatest slope'
write(*, '(e12.4, T20, a)' ) param_acv(7) , 'max_pente/min_pente                   slope anisotropy factor'
write(*, '(e12.4, T20, a)' ) param_acv(8) , 'highest curvature/smallest curvature, curvature anisotropy factor'

deallocate( array )

!-----------------------------------------------------------------------------------------------------------------

allocate( array    (1:nx, 1:ny) )
allocate( array_tmp(1:nx, 1:ny) )

array(1:nx, 1:ny) = 1
call apod( tab_in = array    (1:nx, 1:ny),  &  !
          tab_out = array_tmp(1:nx, 1:ny),  &  !
             long = nx,                     &  !
             larg = ny,                     &  !
         type_apo = 'tuckey' )                 !

call write_surf( nom_fic = "out/apod_tuckey.sur",   &  !
                 tab_s   = array_tmp(1:nx, 1:ny),   &  !
                 scal    = scal_surf )                 !

call apod( tab_in = array    (1:nx, 1:ny),  &  !
          tab_out = array_tmp(1:nx, 1:ny),  &  !
             long = nx,                     &  !
             larg = ny,                     &  !
         type_apo = 'blackm' )                 !

call write_surf( nom_fic = "out/apod_blackman.sur",   &  !
                 tab_s   = array_tmp(1:nx, 1:ny),     &  !
                 scal    = scal_surf )                   !

call apod( tab_in = array    (1:nx, 1:ny),  &  !
          tab_out = array_tmp(1:nx, 1:ny),  &  !
             long = nx,                     &  !
             larg = ny,                     &  !
         type_apo = 'hann__' )                 !

call write_surf( nom_fic = "out/apod_hann.sur",   &  !
                 tab_s   = array_tmp(1:nx, 1:ny),     &  !
                 scal    = scal_surf )                   !

deallocate( array, array_tmp )

!-----------------------------------------------------------------------------------------------------------------

allocate( array_tmp(1:nx, 1:ny) )

call read_surf( nom_fic = "sur/verif_corr.sur",    &  !
                tab_s   = array,                   &  !
                scal    = scal_surf )                 !

nnx = scal_surf%xres
nny = scal_surf%yres

! let suppose that the length along x is 100 µm
lx = 100.e-6_R8 ; dx = lx / nnx ; dy = dx

n_th = omp_get_max_threads()
call fftw_plan_with_nthreads(nthreads = n_th)

nx2 = 2 * ( nint(PAD_FFT * nnx) / 2)
ny2 = 2 * ( nint(PAD_FFT * nny) / 2)

call init_fftw3(long = nx2,  & !
                larg = ny2 )   !

! results gwyddion : correlation length between 3.7 and 4.2 µm
!                    no particular angle
call correlation_parameters( tab       = array(1:nnx, 1:nny),  &  !
                             long      = nnx,                  &  !
                             larg      = nny,                  &  !
                             res       = param_acv(1:8),       &  !
                             cut       = 0.5_R8,               &  !
                             sub_plane = .true.,               &  !
                             scale_xy  = [ dx, dy ],           &  !
                             omp       = .true.  )                !

write(*, '(    a, T20, a)' ) '************', '*******************************************************************'
write(*, '(e12.4, T20, a)' )  param_acv(1) , 'axe_a,                                ellipsis big axis'
write(*, '(e12.4, T20, a)' )  param_acv(2) , 'axe_b,                                ellipsis small axis'
write(*, '(e12.4, T20, a)' )  param_acv(3) , 'axe_a/axe_b                           another anisotropy factor'
write(*, '(e12.4, T20, a)' )  param_acv(4) , 'nint(angle/inc_a),                    main texture orientation'
write(*, '(e12.4, T20, a)' )  param_acv(5) , 'ray_pente,                            radius of greatest slope'
write(*, '(e12.4, T20, a)' )  param_acv(6) , 'max_pente,                            greatest slope'
write(*, '(e12.4, T20, a)' )  param_acv(7) , 'max_pente/min_pente                   slope anisotropy factor'
write(*, '(e12.4, T20, a)' )  param_acv(8) , 'highest curvature/smallest curvature, curvature anisotropy factor'

call end_fftw3()

deallocate( array  )

stop
contains

subroutine fake_acv(acv_array, long, larg, param, ech)
implicit none
real   (kind=R8), intent(out), allocatable, dimension(:,:)  :: acv_array
integer(kind=I4), intent(in)                                :: long, larg
real   (kind=R8), intent(in), dimension(1:3)                :: param
real   (kind=R8), intent(in), dimension(1:2)                :: ech

   integer(kind=I4) :: i, j, ia, ib, i0, j0
   real   (kind=R8) :: a, b, ang, c, s, echx, echy

   allocate( acv_array(1:long, 1:larg) )

   a   = param(1)
   b   = param(2)
   ang = param(3)

   echx = ech(1)
   echy = ech(2)

   ia = int( a / echx )
   ib = int( b / echy )

   i0 = long / 2 + 1
   j0 = larg / 2 + 1

   c = cos( ang * PI_R8 / 180 )
   s = sin( ang * PI_R8 / 180 )

   do j = 1, larg
   do i = 1, long

      acv_array(i, j) = exp( log(1./2.) * sqrt( ( ( +c * (i - i0) + s * (j - j0) ) / ia )**2 +   &  !
                                                ( ( -s * (i - i0) + c * (j - j0) ) / ib )**2 ) )    !

   enddo
   enddo

return
endsubroutine fake_acv

endprogram test_anisotropy

!~ RESULTS ---------------------------------------

!~   0.1997E-03       bmp, maximum over [0,179°] of the peaks mean width
!~   0.1085E-04       smp, minimum over [0,179°] of the peaks mean width
!~   0.1841E+02       rmp, ratio bmp/smp
!~   0.2810E-02       bml, maximum over [0,179°] of the path length
!~   0.2428E-03       sml, minimum over [0,179°] of the path length
!~   0.1157E+02       rml, ratio bml/sml
!~   0.1283E-02       bms, maximum over [0,179°] of the standard deviation of slope
!~   0.3437E-03       sms, minimum over [0,179°] of the standard deviation of slope
!~   0.3734E+01       rms, ratio bms/sms
!~   0.3994E-04       axe_a,                                ellipsis big axis
!~   0.1000E-04       axe_b,                                ellipsis small axis
!~   0.3993E+01       axe_a/axe_b                           another anisotropy factor
!~   0.2500E+02       nint(angle/inc_a),                    main texture orientation
!~   0.5000E-06       ray_pente,                            radius of greatest slope
!~   0.6709E+05       max_pente,                            greatest slope
!~   0.3318E+01       max_pente/min_pente                   slope anisotropy factor
!~   0.3258E-01       highest curvature/smallest curvature, curvature anisotropy factor
!~ ************       *******************************************************************
!~   0.4284E-05       axe_a,                                ellipsis big axis
!~   0.3831E-05       axe_b,                                ellipsis small axis
!~   0.1118E+01       axe_a/axe_b                           another anisotropy factor
!~   0.7700E+02       nint(angle/inc_a),                    main texture orientation
!~   0.3906E-06       ray_pente,                            radius of greatest slope
!~   0.1893E+06       max_pente,                            greatest slope
!~   0.2458E+01       max_pente/min_pente                   slope anisotropy factor
!~  -0.2528E+02       highest curvature/smallest curvature, curvature anisotropy factor

