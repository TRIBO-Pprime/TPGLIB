program test_smooth
use data_arch,     only : I4, R8
use miscellaneous, only : get_unit
use surfile,       only : read_surf, write_surf, init_scal, SCALE_SURF, unit2IUf
use filter,        only : median_smooth, median_filter, soften, morpho_filter, fft_filter, PAD_FFT_FILTER
use fftw3,         only : fftw_plan_with_nthreads, init_fftw3, end_fftw3, extend
!$ use omp_lib

implicit none

real(kind=R8), allocatable, dimension(:,:) :: heights, heights_out, heights_copy, mask, ext_heights

real(kind=R8) :: dx, dy, dz, fft_cutoff, pad

integer(kind=I4) :: nx, ny, nx2, ny2, n_th

type(SCALE_SURF) :: scal_surf, scal_mask

!========================================================================= padding / windowing

call read_surf( nom_fic = "sur/AB-Zinv-ART-8-21-200x200.sur", &  ! IN
                  tab_s = heights,                            &  ! OUT
                   scal = scal_surf )                            ! OUT

nx = scal_surf%xres
ny = scal_surf%yres

pad = 1.5!PAD_FFT_FILTER

nx2 = 2 * ( nint(pad * nx)/2 )
ny2 = 2 * ( nint(pad * ny)/2 )

allocate( ext_heights(1:nx2, 1:ny2) )

call extend(   tab_in = heights(1:nx, 1:ny),          &  !
              tab_out = ext_heights(1:nx2, 1:ny2),    &  !
                   nx = nx,                           &  !
                   ny = ny,                           &  !
                  nx2 = nx2,                          &  !
                  ny2 = ny2,                          &  !
                  ext = 'constant',                   &  !
             type_apo = 'hann__')                        !

scal_surf%xres = nx2
scal_surf%yres = ny2

call write_surf( nom_fic = "out/test_extend.sur",        &  !
                 tab_s   = ext_heights(1:nx2, 1:ny2),    &  !
                 scal    = scal_surf )                      !

deallocate( ext_heights )

!========================================================================= fft_filter 5 µm

call read_surf( nom_fic = "sur/AB-Zinv-ART-8-21-200x200.sur", &  ! IN
                  tab_s = heights,                            &  ! OUT
                   scal = scal_surf )                            ! OUT

nx = scal_surf%xres
ny = scal_surf%yres

dx = scal_surf%dx * unit2IUf(scal_surf%dx_unit)
dy = scal_surf%dy * unit2IUf(scal_surf%dy_unit)
dz = 0

write(*,*) 'nx, ny = ', nx, ny
write(*,*) 'dx, dy = ', dx, dy

allocate( heights_out(1:nx, 1:ny) )


n_th = omp_get_max_threads()
call fftw_plan_with_nthreads(nthreads = n_th)

nx2 = 2 * ( nint(PAD_FFT_FILTER * nx)/2 )
ny2 = 2 * ( nint(PAD_FFT_FILTER * ny)/2 )

fft_cutoff = dx / 5.e-6

call init_fftw3 (long      = nx2, & !
                 larg      = ny2)   !

call fft_filter(tab       = heights(1:nx, 1:ny),      & !
                long      = nx,                       & !
                larg      = ny,                       & !
                cutoff    = fft_cutoff,               & !
                bf_tab    = heights_out(1:nx, 1:ny),  & !
                multi_fft = .FALSE.)                    !

call write_surf( nom_fic = "out/test_fft_filter_005µm.sur",    &  !
                 tab_s   = heights_out(1:nx, 1:ny),            &  !
                 scal    = scal_surf )                            !

write(*,*) 'gaussian filter cutoff = 5 µm'

!========================================================================= fft_filter 80 µm

fft_cutoff = dx / 80.e-6

call fft_filter(tab       = heights(1:nx, 1:ny),      & !
                long      = nx,                       & !
                larg      = ny,                       & !
                cutoff    = fft_cutoff,               & !
                bf_tab    = heights_out(1:nx, 1:ny),  & !
                multi_fft = .FALSE.)                    !

call write_surf( nom_fic = "out/test_fft_filter_080µm.sur",    &  !
                 tab_s   = heights_out(1:nx, 1:ny),            &  !
                 scal    = scal_surf )                            !

write(*,*) 'gaussian filter cutoff = 80 µm'

!========================================================================= fft_filter 150 µm

fft_cutoff = dx / 150.e-6

call fft_filter(tab       = heights(1:nx, 1:ny),      & !
                long      = nx,                       & !
                larg      = ny,                       & !
                cutoff    = fft_cutoff,               & !
                bf_tab    = heights_out(1:nx, 1:ny),  & !
                multi_fft = .FALSE.)                    !

call write_surf( nom_fic = "out/test_fft_filter_150µm.sur",    &  !
                 tab_s   = heights_out(1:nx, 1:ny),            &  !
                 scal    = scal_surf )                            !

write(*,*) 'gaussian filter cutoff = 150 µm'

call end_fftw3()

!========================================================================= ROLL +, RAY = 10 µm

call morpho_filter( mtype       = "dilation",                 &  ! IN
                    tabin       = heights(1:nx, 1:ny),        &  ! IN
                    tabou       = heights_out(1:nx, 1:ny),    &  ! OUT
                    long        = nx,                         &  ! IN
                    larg        = ny,                         &  ! IN
                    scale_xyz   = [dx, dy, dz],               &  ! IN
                    ray         = 10.e-6_R8,                  &  ! IN
                    omp         = .true.,                     &  ! IN
                    nb_div      = 50 )                           ! IN

call write_surf( nom_fic = "out/test_roll_smooth_dilation_10µm.sur",   &  !
                 tab_s   = heights_out(1:nx, 1:ny),                    &  !
                 scal    = scal_surf )                                    !

write(*,*) 'morpho filter: dilation disk 10 µm'

!------------------------------------------------------------------------- ROLL -, RAY = 10 µm

call morpho_filter( mtype       = "erosion",                  &  ! IN
                    tabin       = heights(1:nx, 1:ny),        &  ! IN
                    tabou       = heights_out(1:nx, 1:ny),    &  ! OUT
                    long        = nx,                         &  ! IN
                    larg        = ny,                         &  ! IN
                    scale_xyz   = [dx, dy, dz],               &  ! IN
                    ray         = 10.e-6_R8,                  &  ! IN
                    omp         = .true.,                     &  ! IN
                    nb_div      = 50 )                           ! IN

call write_surf( nom_fic = "out/test_roll_smooth_erosion_10µm.sur",   &  !
                 tab_s   = heights_out(1:nx, 1:ny),                   &  !
                 scal    = scal_surf )                                   !

write(*,*) 'morpho filter: erosion disk 10 µm'

!------------------------------------------------------------------------- CLOSING, RAY = 10 µm

call morpho_filter( mtype       = "closing",                  &  ! IN
                    tabin       = heights(1:nx, 1:ny),        &  ! IN
                    tabou       = heights_out(1:nx, 1:ny),    &  ! OUT
                    long        = nx,                         &  ! IN
                    larg        = ny,                         &  ! IN
                    scale_xyz   = [dx, dy, dz],               &  ! IN
                    ray         = 10.e-6_R8,                  &  ! IN
                    omp         = .true.,                     &  ! IN
                    nb_div      = 50 )                           ! IN

call write_surf( nom_fic = "out/test_roll_smooth_closing_10µm.sur",   &  !
                 tab_s   = heights_out(1:nx, 1:ny),                   &  !
                 scal    = scal_surf )                                   !

write(*,*) 'morpho filter: closing disk 10 µm'

!------------------------------------------------------------------------- OPENING, RAY = 10 µm

call morpho_filter( mtype       = "opening",                  &  ! IN
                    tabin       = heights(1:nx, 1:ny),        &  ! IN
                    tabou       = heights_out(1:nx, 1:ny),    &  ! OUT
                    long        = nx,                         &  ! IN
                    larg        = ny,                         &  ! IN
                    scale_xyz   = [dx, dy, dz],               &  ! IN
                    ray         = 10.e-6_R8,                  &  ! IN
                    omp         = .true.,                     &  ! IN
                    nb_div      = 50 )                           ! IN

call write_surf( nom_fic = "out/test_roll_smooth_opening_10µm.sur",   &  !
                 tab_s   = heights_out(1:nx, 1:ny),                   &  !
                 scal    = scal_surf )                                   !

write(*,*) 'morpho filter: opening disk 10 µm'

!========================================================================= KERNEL SIZE = 3

allocate( heights_copy(1:nx, 1:ny) )

heights_copy(1:nx, 1:ny) = heights(1:nx, 1:ny)

call median_smooth( tab    = heights(1:nx, 1:ny),        &  ! INOUT
                    long   = nx,                         &  ! IN
                    larg   = ny,                         &  ! IN
                    kernel = 3,                          &  ! IN
                    omp    = .true. )                       ! IN

call write_surf( nom_fic = "out/test_med_smooth_3.sur",  &  !
                 tab_s   = heights(1:nx, 1:ny),          &  !
                 scal    = scal_surf )                      !

write(*,*) 'median filter kernel: 3'

!========================================================================= KERNEL SIZE = 5

heights(1:nx, 1:ny) = heights_copy(1:nx, 1:ny)

call median_smooth( tab    = heights(1:nx, 1:ny),        &  ! INOUT
                    long   = nx,                         &  ! IN
                    larg   = ny,                         &  ! IN
                    kernel = 5,                          &  ! IN
                    omp    = .true. )                       ! IN

call write_surf( nom_fic = "out/test_med_smooth_5.sur",  &  !
                 tab_s   = heights(1:nx, 1:ny),          &  !
                 scal    = scal_surf )                      !

write(*,*) 'median filter kernel: 5'

!========================================================================= KERNEL SIZE = 9

heights(1:nx, 1:ny) = heights_copy(1:nx, 1:ny)

call median_smooth( tab    = heights(1:nx, 1:ny),        &  ! INOUT
                    long   = nx,                         &  ! IN
                    larg   = ny,                         &  ! IN
                    kernel = 9,                          &  ! IN
                    omp    = .true. )                       ! IN

call write_surf( nom_fic = "out/test_med_smooth_9.sur",  &  !
                 tab_s   = heights(1:nx, 1:ny),          &  !
                 scal    = scal_surf )                      !

write(*,*) 'median filter kernel: 9'

!========================================================================= SOFTEN NO MASK

heights(1:nx, 1:ny) = heights_copy(1:nx, 1:ny)

call soften( tabin  = heights(1:nx, 1:ny),         &  ! IN
             tabout = heights_out(1:nx, 1:ny),     &  ! OUT
             long   = nx,                          &  ! IN
             larg   = ny )                            ! IN

call write_surf( nom_fic = "out/test_soften_no_mask.sur",   &  !
                 tab_s   = heights_out(1:nx, 1:ny),         &  !
                 scal    = scal_surf )                         !

write(*,*) 'simple smooth with no mask'

!========================================================================= SOFTEN WITH MASK

call read_surf( nom_fic = "sur/test_mask.sur",  &  ! IN
                  tab_s = mask,                 &  ! OUT
                   scal = scal_mask )              ! OUT

call soften( tabin  = heights(1:nx, 1:ny),         &  ! IN
             tabout = heights_out(1:nx, 1:ny),     &  ! OUT
             mask   = nint(mask(1:nx, 1:ny)),      &  ! IN
             long   = nx,                          &  ! IN
             larg   = ny )                            ! IN

call write_surf( nom_fic = "out/test_soften_with_mask.sur",   &  !
                 tab_s   = heights_out(1:nx, 1:ny),           &  !
                 scal    = scal_surf )                           !

write(*,*) 'simple smooth in upper mask'

!========================================================================= MEDIAN FILTER SNB = 10, KERN = 15

call median_filter( tab    = heights(1:nx, 1:ny),        &  ! INOUT
                    long   = nx,                         &  ! IN
                    larg   = ny,                         &  ! IN
                    kernel = 15,                         &  ! IN
                    snb    = 10,                         &  ! IN
                    sig    = 3._R8,                      &  ! IN
                    omp    = .true. )                       ! IN

call write_surf( nom_fic = "out/test_med_filt_k15_s10.sur", &  !
                 tab_s   = heights(1:nx, 1:ny),             &  !
                 scal    = scal_surf )                         !

deallocate( heights, heights_copy, heights_out, mask )

stop
endprogram test_smooth
