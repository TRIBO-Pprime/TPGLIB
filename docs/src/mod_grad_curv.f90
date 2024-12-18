!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: may, 3 2024
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Routines to calculate surface gradients and curvatures**
!<  </span>

module grad_curv
use data_arch,   only : I4, R8, UN, PI_R8
use sort_arrays, only : sort_array2
use surfile,     only : read_surf, write_surf, SCALE_SURF, unit2IUf, init_scal
use filter,      only : fft_filter
use fftw3,       only : fftw_plan_with_nthreads, init_fftw3, end_fftw3, PAD_FFT
!$ use omp_lib

implicit none

private

public :: gradient, gauss_curv, curv2, curvature, label_surf_summits, peaks_and_pits_curvatures
public :: test_labelize_point, test_label_surf_summits, test_peaks_and_pits_curvatures

contains

   subroutine gradient(tab, nx, ny, dx, dy, gradx, grady)
   !================================================================================================
   !< @note Function to calculate the gradient of a 2D array
   !<
   !<  It implements the details given in ISO 25178.
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in ) :: nx                            !! *number of pixels along x*
   integer(kind=I4), intent(in ) :: ny                            !! *number of pixels along x*
   real   (kind=R8), intent(in ) :: dx                            !! *x lag*
   real   (kind=R8), intent(in ) :: dy                            !! *y lag*
   real   (kind=R8), intent(in ), dimension(1:nx, 1:ny) :: tab    !! *Input 2D array*
   real   (kind=R8), intent(out), dimension(1:nx, 1:ny) :: gradx  !! *derivative along x 2D array*
   real   (kind=R8), intent(out), dimension(1:nx, 1:ny) :: grady  !! *derivative along y 2D array*

      integer(kind=I4) :: i, j

      !------------------------------------------------------------------ GRADX
      i = 1

      gradx(1, 1:ny) = (UN / (60 * dx) ) * ( -147 * tab(i+0, 1:ny)   &  !
                                             +360 * tab(i+1, 1:ny)   &  !
                                             -450 * tab(i+2, 1:ny)   &  !
                                             +400 * tab(i+3, 1:ny)   &  !
                                             -225 * tab(i+4, 1:ny)   &  !
                                             +072 * tab(i+5, 1:ny)   &  !
                                             -010 * tab(i+6, 1:ny) )    !

      i = 1

      gradx(2, 1:ny) = (UN / (60 * dx) ) * ( -010 * tab(i+0, 1:ny)   &  !
                                             -077 * tab(i+1, 1:ny)   &  !
                                             +150 * tab(i+2, 1:ny)   &  !
                                             -100 * tab(i+3, 1:ny)   &  !
                                             +050 * tab(i+4, 1:ny)   &  !
                                             -015 * tab(i+5, 1:ny)   &  !
                                             +002 * tab(i+6, 1:ny) )    !

      i = 1

      gradx(3, 1:ny) = (UN / (60 * dx) ) * ( +002 * tab(i+0, 1:ny)   &  !
                                             -024 * tab(i+1, 1:ny)   &  !
                                             -035 * tab(i+2, 1:ny)   &  !
                                             +080 * tab(i+3, 1:ny)   &  !
                                             -030 * tab(i+4, 1:ny)   &  !
                                             +008 * tab(i+5, 1:ny)   &  !
                                             -001 * tab(i+6, 1:ny) )    !

      do i = 4, nx - 3
         gradx(i, 1:ny) = (UN / (60 * dx) ) * ( +01 * tab(i+3, 1:ny)   &  !
                                                -09 * tab(i+2, 1:ny)   &  !
                                                +45 * tab(i+1, 1:ny)   &  !
                                                -45 * tab(i-1, 1:ny)   &  !
                                                +09 * tab(i-2, 1:ny)   &  !
                                                -01 * tab(i-3, 1:ny) )    !
      enddo

      i = nx

      gradx(nx  , 1:ny) = -(UN / (60 * dx) ) * ( -147 * tab(i-0, 1:ny)   &  !
                                                 +360 * tab(i-1, 1:ny)   &  !
                                                 -450 * tab(i-2, 1:ny)   &  !
                                                 +400 * tab(i-3, 1:ny)   &  !
                                                 -225 * tab(i-4, 1:ny)   &  !
                                                 +072 * tab(i-5, 1:ny)   &  !
                                                 -010 * tab(i-6, 1:ny) )    !

      i = nx

      gradx(nx-1, 1:ny) = -(UN / (60 * dx) ) * ( -010 * tab(i-0, 1:ny)   &  !
                                                 -077 * tab(i-1, 1:ny)   &  !
                                                 +150 * tab(i-2, 1:ny)   &  !
                                                 -100 * tab(i-3, 1:ny)   &  !
                                                 +050 * tab(i-4, 1:ny)   &  !
                                                 -015 * tab(i-5, 1:ny)   &  !
                                                 +002 * tab(i-6, 1:ny) )    !

      i = nx

      gradx(nx-2, 1:ny) = -(UN / (60 * dx) ) * ( +002 * tab(i-0, 1:ny)   &  !
                                                 -024 * tab(i-1, 1:ny)   &  !
                                                 -035 * tab(i-2, 1:ny)   &  !
                                                 +080 * tab(i-3, 1:ny)   &  !
                                                 -030 * tab(i-4, 1:ny)   &  !
                                                 +008 * tab(i-5, 1:ny)   &  !
                                                 -001 * tab(i-6, 1:ny) )    !

      !------------------------------------------------------------------ GRADY

      j = 1

      grady(1:nx, 1) = (UN / (60 * dy) ) * ( -147 * tab(1:nx, j+0)   &  !
                                             +360 * tab(1:nx, j+1)   &  !
                                             -450 * tab(1:nx, j+2)   &  !
                                             +400 * tab(1:nx, j+3)   &  !
                                             -225 * tab(1:nx, j+4)   &  !
                                             +072 * tab(1:nx, j+5)   &  !
                                             -010 * tab(1:nx, j+6) )    !

      j = 1

      grady(1:nx, 2) = (UN / (60 * dy) ) * ( -010 * tab(1:nx, j+0)   &  !
                                             -077 * tab(1:nx, j+1)   &  !
                                             +150 * tab(1:nx, j+2)   &  !
                                             -100 * tab(1:nx, j+3)   &  !
                                             +050 * tab(1:nx, j+4)   &  !
                                             -015 * tab(1:nx, j+5)   &  !
                                             +002 * tab(1:nx, j+6) )    !

      j = 1

      grady(1:nx, 3) = (UN / (60 * dy) ) * ( +002 * tab(1:nx, j+0)   &  !
                                             -024 * tab(1:nx, j+1)   &  !
                                             -035 * tab(1:nx, j+2)   &  !
                                             +080 * tab(1:nx, j+3)   &  !
                                             -030 * tab(1:nx, j+4)   &  !
                                             +008 * tab(1:nx, j+5)   &  !
                                             -001 * tab(1:nx, j+6) )    !

      do j = 4, ny - 3
         grady(1:nx, j) = (UN / (60 * dy) ) * ( +01 * tab(1:nx, j+3)   &  !
                                                -09 * tab(1:nx, j+2)   &  !
                                                +45 * tab(1:nx, j+1)   &  !
                                                -45 * tab(1:nx, j-1)   &  !
                                                +09 * tab(1:nx, j-2)   &  !
                                                -01 * tab(1:nx, j-3) )    !
      enddo

      j = ny

      grady(1:nx, ny  ) = -(UN / (60 * dy) ) * ( -147 * tab(1:nx, j-0)   &  !
                                                 +360 * tab(1:nx, j-1)   &  !
                                                 -450 * tab(1:nx, j-2)   &  !
                                                 +400 * tab(1:nx, j-3)   &  !
                                                 -225 * tab(1:nx, j-4)   &  !
                                                 +072 * tab(1:nx, j-5)   &  !
                                                 -010 * tab(1:nx, j-6) )    !

      j = ny

      grady(1:nx, ny-1) = -(UN / (60 * dy) ) * ( -010 * tab(1:nx, j-0)   &  !
                                                 -077 * tab(1:nx, j-1)   &  !
                                                 +150 * tab(1:nx, j-2)   &  !
                                                 -100 * tab(1:nx, j-3)   &  !
                                                 +050 * tab(1:nx, j-4)   &  !
                                                 -015 * tab(1:nx, j-5)   &  !
                                                 +002 * tab(1:nx, j-6) )    !

      j = ny

      grady(1:nx, ny-2) = -(UN / (60 * dy) ) * ( +002 * tab(1:nx, j-0)   &  !
                                                 -024 * tab(1:nx, j-1)   &  !
                                                 -035 * tab(1:nx, j-2)   &  !
                                                 +080 * tab(1:nx, j-3)   &  !
                                                 -030 * tab(1:nx, j-4)   &  !
                                                 +008 * tab(1:nx, j-5)   &  !
                                                 -001 * tab(1:nx, j-6) )    !

   return
   endsubroutine gradient


   subroutine gauss_curv(gradx, grady, nx, ny, dx, dy, gradxx, gradyy, gradxy)
   !================================================================================================
   !! Function to calculate the double derivatives of a 2D array
   implicit none
   integer(kind=I4), intent(in ) :: nx                              !! *number of pixels along x*
   integer(kind=I4), intent(in ) :: ny                              !! *number of pixels along x*
   real   (kind=R8), intent(in ) :: dx                              !! *x lag*
   real   (kind=R8), intent(in ) :: dy                              !! *y lag*
   real   (kind=R8), intent(in ), dimension(1:nx, 1:ny) :: gradx    !! *derivative along x 2D array*
   real   (kind=R8), intent(in ), dimension(1:nx, 1:ny) :: grady    !! *derivative along y 2D array*
   real   (kind=R8), intent(out), dimension(1:nx, 1:ny) :: gradxx   !! *double derivative along x, x 2D array*
   real   (kind=R8), intent(out), dimension(1:nx, 1:ny) :: gradyy   !! *double derivative along y, y 2D array*
   real   (kind=R8), intent(out), dimension(1:nx, 1:ny) :: gradxy   !! *double derivative along x, y 2D array*

      call gradient(tab   = gradx (1:nx, 1:ny), &  ! IN
                    gradx = gradxx(1:nx, 1:ny), &  ! OUT
                    grady = gradxy(1:nx, 1:ny), &  ! OUT
                    nx    = nx,                 &  ! IN
                    ny    = ny,                 &  ! IN
                    dx    = dx,                 &  ! IN
                    dy    = dy)                    ! IN

      call gradient(tab   = grady (1:nx, 1:ny), &  !IN
                    gradx = gradxy(1:nx, 1:ny), &  !OUT
                    grady = gradyy(1:nx, 1:ny), &  !OUT
                    nx    = nx,                 &  !IN
                    ny    = ny,                 &  !IN
                    dx    = dx,                 &  !IN
                    dy    = dy)                    !IN

   return
   endsubroutine gauss_curv


   subroutine curv2(tab, nx, ny, dx, dy, gradxx, gradyy)
   !================================================================================================
   !< @note Function to calculate the double derivatives of a 2D array
   !<
   !<  It implements the details given in ISO 25178.
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in ) :: nx                               !! *number of pixels along x*
   integer(kind=I4), intent(in ) :: ny                               !! *number of pixels along x*
   real   (kind=R8), intent(in ) :: dx                               !! *x lag*
   real   (kind=R8), intent(in ) :: dy                               !! *y lag*
   real   (kind=R8), intent(in ), dimension(1:nx, 1:ny) :: tab       !! *input 2D array*
   real   (kind=R8), intent(out), dimension(1:nx, 1:ny) :: gradxx    !! *double derivative along x, x 2D array*
   real   (kind=R8), intent(out), dimension(1:nx, 1:ny) :: gradyy    !! *double derivative along y, y 2D array*

      integer(kind=I4) :: i, j

      !------------------------------------------------------------------ GRADXX
      i = 1

      gradxx(1, 1:ny) = (UN / (180 * dx**2) ) * ( +0812 * tab(i+0, 1:ny)   &  !
                                                  -3132 * tab(i+1, 1:ny)   &  !
                                                  +5265 * tab(i+2, 1:ny)   &  !
                                                  -5080 * tab(i+3, 1:ny)   &  !
                                                  +2970 * tab(i+4, 1:ny)   &  !
                                                  -0972 * tab(i+5, 1:ny)   &  !
                                                  +0137 * tab(i+6, 1:ny) )    !

      i = 1

      gradxx(2, 1:ny) = (UN / (180 * dx**2) ) * ( +0137 * tab(i+0, 1:ny)   &  !
                                                  -0147 * tab(i+1, 1:ny)   &  !
                                                  -0255 * tab(i+2, 1:ny)   &  !
                                                  +0470 * tab(i+3, 1:ny)   &  !
                                                  -0285 * tab(i+4, 1:ny)   &  !
                                                  +0093 * tab(i+5, 1:ny)   &  !
                                                  -0013 * tab(i+6, 1:ny) )    !

      i = 1

      gradxx(3, 1:ny) = (UN / (180 * dx**2) ) * ( -0013 * tab(i+0, 1:ny)   &  !
                                                  +0228 * tab(i+1, 1:ny)   &  !
                                                  -0420 * tab(i+2, 1:ny)   &  !
                                                  +0200 * tab(i+3, 1:ny)   &  !
                                                  +0015 * tab(i+4, 1:ny)   &  !
                                                  -0012 * tab(i+5, 1:ny)   &  !
                                                  +0002 * tab(i+6, 1:ny) )    !

      do i = 4, nx - 3
         gradxx(i, 1:ny) = (UN / (180 * dx**2) ) * ( +002 * tab(i+3, 1:ny)   &  !
                                                     -027 * tab(i+2, 1:ny)   &  !
                                                     +270 * tab(i+1, 1:ny)   &  !
                                                     -490 * tab(i  , 1:ny)   &  !
                                                     +270 * tab(i-1, 1:ny)   &  !
                                                     -027 * tab(i-2, 1:ny)   &  !
                                                     +002 * tab(i-3, 1:ny) )    !
      enddo

      i = nx

      gradxx(nx  , 1:ny) = (UN / (180 * dx**2) ) * ( +0812 * tab(i-0, 1:ny)   &  !
                                                     -3132 * tab(i-1, 1:ny)   &  !
                                                     +5265 * tab(i-2, 1:ny)   &  !
                                                     -5080 * tab(i-3, 1:ny)   &  !
                                                     +2970 * tab(i-4, 1:ny)   &  !
                                                     -0972 * tab(i-5, 1:ny)   &  !
                                                     +0137 * tab(i-6, 1:ny) )    !

      i = nx

      gradxx(nx-1, 1:ny) = (UN / (180 * dx**2) ) * ( +0137 * tab(i-0, 1:ny)   &  !
                                                     -0147 * tab(i-1, 1:ny)   &  !
                                                     -0255 * tab(i-2, 1:ny)   &  !
                                                     +0470 * tab(i-3, 1:ny)   &  !
                                                     -0285 * tab(i-4, 1:ny)   &  !
                                                     +0093 * tab(i-5, 1:ny)   &  !
                                                     -0013 * tab(i-6, 1:ny) )    !

      i = nx

      gradxx(nx-2, 1:ny) = (UN / (180 * dx**2) ) * ( -0013 * tab(i-0, 1:ny)   &  !
                                                     +0228 * tab(i-1, 1:ny)   &  !
                                                     -0420 * tab(i-2, 1:ny)   &  !
                                                     +0200 * tab(i-3, 1:ny)   &  !
                                                     +0015 * tab(i-4, 1:ny)   &  !
                                                     -0012 * tab(i-5, 1:ny)   &  !
                                                     +0002 * tab(i-6, 1:ny) )    !

      !------------------------------------------------------------------ GRADYY
      j = 1

      gradyy(1:nx, 1) = (UN / (180 * dy**2) ) * ( +0812 * tab(1:nx, j+0)   &  !
                                                  -3132 * tab(1:nx, j+1)   &  !
                                                  +5265 * tab(1:nx, j+2)   &  !
                                                  -5080 * tab(1:nx, j+3)   &  !
                                                  +2970 * tab(1:nx, j+4)   &  !
                                                  -0972 * tab(1:nx, j+5)   &  !
                                                  +0137 * tab(1:nx, j+6) )    !

      j = 1

      gradyy(1:nx, 2) = (UN / (180 * dy**2) ) * ( +0137 * tab(1:nx, j+0)   &  !
                                                  -0147 * tab(1:nx, j+1)   &  !
                                                  -0255 * tab(1:nx, j+2)   &  !
                                                  +0470 * tab(1:nx, j+3)   &  !
                                                  -0285 * tab(1:nx, j+4)   &  !
                                                  +0093 * tab(1:nx, j+5)   &  !
                                                  -0013 * tab(1:nx, j+6) )    !

      j = 1

      gradyy(1:nx, 3) = (UN / (180 * dy**2) ) * ( -0013 * tab(1:nx, j+0)   &  !
                                                  +0228 * tab(1:nx, j+1)   &  !
                                                  -0420 * tab(1:nx, j+2)   &  !
                                                  +0200 * tab(1:nx, j+3)   &  !
                                                  +0015 * tab(1:nx, j+4)   &  !
                                                  -0012 * tab(1:nx, j+5)   &  !
                                                  +0002 * tab(1:nx, j+6) )    !

      do j = 4, ny - 3
         gradyy(1:nx, j) = (UN / (180 * dy**2) ) * ( +002 * tab(1:nx, j+3)   &  !
                                                     -027 * tab(1:nx, j+2)   &  !
                                                     +270 * tab(1:nx, j+1)   &  !
                                                     -490 * tab(1:nx, j  )   &  !
                                                     +270 * tab(1:nx, j-1)   &  !
                                                     -027 * tab(1:nx, j-2)   &  !
                                                     +002 * tab(1:nx, j-3) )    !
      enddo

      j = ny

      gradyy(1:nx, ny  ) = (UN / (180 * dy**2) ) * ( +0812 * tab(1:nx, j-0)   &  !
                                                     -3132 * tab(1:nx, j-1)   &  !
                                                     +5265 * tab(1:nx, j-2)   &  !
                                                     -5080 * tab(1:nx, j-3)   &  !
                                                     +2970 * tab(1:nx, j-4)   &  !
                                                     -0972 * tab(1:nx, j-5)   &  !
                                                     +0137 * tab(1:nx, j-6) )    !

      j = ny

      gradyy(1:nx, ny-1) = (UN / (180 * dy**2) ) * ( +0137 * tab(1:nx, j-0)   &  !
                                                     -0147 * tab(1:nx, j-1)   &  !
                                                     -0255 * tab(1:nx, j-2)   &  !
                                                     +0470 * tab(1:nx, j-3)   &  !
                                                     -0285 * tab(1:nx, j-4)   &  !
                                                     +0093 * tab(1:nx, j-5)   &  !
                                                     -0013 * tab(1:nx, j-6) )    !

      j = ny

      gradyy(1:nx, ny-2) = (UN / (180 * dy**2) ) * ( -0013 * tab(1:nx, j-0)   &  !
                                                     +0228 * tab(1:nx, j-1)   &  !
                                                     -0420 * tab(1:nx, j-2)   &  !
                                                     +0200 * tab(1:nx, j-3)   &  !
                                                     +0015 * tab(1:nx, j-4)   &  !
                                                     -0012 * tab(1:nx, j-5)   &  !
                                                     +0002 * tab(1:nx, j-6) )    !


   return
   endsubroutine curv2


   subroutine curvature(tab, nx, ny, dx, dy, S_param_grad, S_param_curv, gcurvt)
   !================================================================================================
   !! Function to calculate the gaussian curvature of a 2D array,
   !!        its mean quadratic value and the gradient mean quadratic value
   implicit none
   integer(kind=I4), intent(in ) :: nx                              !! *number of pixels along x*
   integer(kind=I4), intent(in ) :: ny                              !! *number of pixels along x*
   real   (kind=R8), intent(in ) :: dx                              !! *x lag*
   real   (kind=R8), intent(in ) :: dy                              !! *y lag*
   real   (kind=R8), intent(out) :: S_param_grad                    !! *mean quadratic gradient value*
   real   (kind=R8), intent(out) :: S_param_curv                    !! *mean quadratic curvature value*
   real   (kind=R8), intent(in ), dimension(1:nx, 1:ny) :: tab      !! *input 2D array*
   real   (kind=R8), intent(out), dimension(1:nx, 1:ny) :: gcurvt   !! *gaussian curvature  2D array*

      real(kind=R8), allocatable, dimension(:,:) :: gradx    !! *derivative along x 2D array*
      real(kind=R8), allocatable, dimension(:,:) :: grady    !! *derivative along y 2D array*
      real(kind=R8), allocatable, dimension(:,:) :: gradxx   !! *double derivative along x, x 2D array*
      real(kind=R8), allocatable, dimension(:,:) :: gradyy   !! *double derivative along y, y 2D array*
      real(kind=R8), allocatable, dimension(:,:) :: gradxy   !! *double derivative along x, y 2D array*

      allocate( gradx (1:nx, 1:ny) )
      allocate( grady (1:nx, 1:ny) )
      allocate( gradxx(1:nx, 1:ny) )
      allocate( gradxy(1:nx, 1:ny) )
      allocate( gradyy(1:nx, 1:ny) )

      call gradient(tab   = tab(1:nx, 1:ny),          & ! IN
                    nx    = nx,                       & ! IN
                    ny    = ny,                       & ! IN
                    dx    = dx,                       & ! IN
                    dy    = dy,                       & ! IN
                    gradx = gradx(1:nx, 1:ny),        & ! OUT
                    grady = grady(1:nx, 1:ny))          ! OUT

!~       S_param_grad = sum( sqrt( gradx(1:nx, 1:ny)**2 + grady(1:nx, 1:ny)**2 ) ) / (nx * ny)
      S_param_grad = sqrt( sum( gradx(1:nx, 1:ny)**2 + grady(1:nx, 1:ny)**2 ) ) / (nx * ny)

      call gauss_curv(gradx  = gradx (1:nx, 1:ny),    & ! IN
                      grady  = grady (1:nx, 1:ny),    & ! IN
                      gradxx = gradxx(1:nx, 1:ny),    & ! OUT
                      gradxy = gradxy(1:nx, 1:ny),    & ! OUT
                      gradyy = gradyy(1:nx, 1:ny),    & ! OUT
                      nx     = nx,                    & ! IN
                      ny     = ny,                    & ! IN
                      dx     = dx,                    & ! IN
                      dy     = dy)                      ! IN

      gcurvt(1:nx, 1:ny) = ( gradxx(1:nx, 1:ny) * gradyy(1:nx, 1:ny) - gradxy(1:nx, 1:ny)**2 ) / ( UN + gradx(1:nx, 1:ny)**2 + grady(1:nx, 1:ny)**2 )**2

!~       S_param_curv = sum( sqrt( gradxx(1:nx, 1:ny)**2 + gradyy(1:nx, 1:ny)**2 ) ) / (nx * ny)
      S_param_curv = sqrt( sum( gradxx(1:nx, 1:ny)**2 + gradyy(1:nx, 1:ny)**2 ) ) / (nx * ny)

      deallocate( gradx, grady, gradxy, gradxx, gradyy )

   return
   endsubroutine curvature


   subroutine deriv_N(x, y, mat_d)
   !================================================================================================
   !< @note Function to provide the interpolation functions of a QU9 element, as well as its derivatives
   !<
   !<  It implements the details given in code ASTER r3.01.01.pdf doc
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real(kind=R8), intent(in ) :: x                                !! *abscissa between -1 and +1*
   real(kind=R8), intent(in ) :: y                                !! *ordinate between -1 and +1*
   real(kind=R8), intent(out), dimension(1:9, 1:6) :: mat_d       !! *array containing N, dN/di, d2N/di2*

      real(kind=R8) :: xm1, xp1, xm12, xp12, umx2, dxm1, dxp1, xy
      real(kind=R8) :: ym1, yp1, ym12, yp12, umy2, dyp1, dym1

      xm1 = x - 1 ; xp1 = x + 1 ; xm12 = x - 1/2._R8 ; xp12 = x + 1/2._R8 ; umx2 = 1 - x**2 ; dxm1 = 2*x - 1 ; dxp1 = 2*x + 1 ; xy = x*y
      ym1 = y - 1 ; yp1 = y + 1 ; ym12 = y - 1/2._R8 ; yp12 = y + 1/2._R8 ; umy2 = 1 - y**2 ; dym1 = 2*y - 1 ; dyp1 = 2*y + 1

!     Nodes order:
!
!     4---7---3
!     |   |   |
!     8---9-- 6
!     |   |   |
!     1---5---2

!     mat_d(1:6, .) = [ N(1:9), dNdx(1:9), dNdy(1:9), d2Ndx2(1:9), d2Ndxdy(1:9), d2Ndy2(1:9) ]

      mat_d(1:9, 1) = [ xy*xm1*ym1/4  , &   ! 1
                        xy*xp1*ym1/4  , &   ! 2
                        xy*xp1*yp1/4  , &   ! 3
                        xy*xm1*yp1/4  , &   ! 4
                        y*umx2*ym1/2  , &   ! 5
                        x*umy2*xp1/2  , &   ! 6
                        y*umx2*yp1/2  , &   ! 7
                        x*umy2*xm1/2  , &   ! 8
                          umx2*umy2 ]       ! 9

      mat_d(1:9, 2) = [ dxm1*y*ym1/4  , &   ! 1
                        dxp1*y*ym1/4  , &   ! 2
                        dxp1*y*yp1/4  , &   ! 3
                        dxm1*y*yp1/4  , &   ! 4
                        -xy*ym1       , &   ! 5
                        dxp1*umy2/2   , &   ! 6
                        -xy*yp1       , &   ! 7
                        dxm1*umy2/2   , &   ! 8
                        -2*x*umy2 ]         ! 9

      mat_d(1:9, 3) = [ x*xm1*dym1/4  , &   ! 1
                        x*xp1*dym1/4  , &   ! 2
                        x*xp1*dyp1/4  , &   ! 3
                        x*xm1*dyp1/4  , &   ! 4
                         umx2*dym1/2  , &   ! 5
                        -xy*xp1       , &   ! 6
                         umx2*dyp1/2  , &   ! 7
                        -xy*xm1       , &   ! 8
                        -2*y*umx2 ]         ! 9

      mat_d(1:9, 4) = [ y*ym1/2       , &   ! 1
                        y*ym1/2       , &   ! 2
                        y*yp1/2       , &   ! 3
                        y*yp1/2       , &   ! 4
                        -y*ym1        , &   ! 5
                           umy2       , &   ! 6
                        -y*yp1        , &   ! 7
                           umy2       , &   ! 8
                        -2*umy2 ]           ! 9

      mat_d(1:9, 5) = [ xm12*ym12     , &   ! 1
                        xp12*ym12     , &   ! 2
                        xp12*yp12     , &   ! 3
                        xm12*yp12     , &   ! 4
                        -x*dym1       , &   ! 5
                        -y*dxp1       , &   ! 6
                        -x*dyp1       , &   ! 7
                        -y*dxm1       , &   ! 8
                        +4*xy ]             ! 9

      mat_d(1:9, 6) = [ x*xm1/2       , &   ! 1
                        x*xp1/2       , &   ! 2
                        x*xp1/2       , &   ! 3
                        x*xm1/2       , &   ! 4
                          +umx2       , &   ! 5
                        -x*xp1        , &   ! 6
                          +umx2       , &   ! 7
                        -x*xm1        , &   ! 8
                        -2*umx2 ]           ! 9

   return
   endsubroutine deriv_N


   subroutine gradient_corner(hgt, gdx, gdy)
   !================================================================================================
   !< @note Function that gives the nodal height gradients
   !<
   !< gdx : (2,1)----(2,2)----(2,3) : QU9 notation
   !<
   !< gdy : (1,2)----(2,2)----(3,2) : QU9 notation
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real   (kind=R8), intent(in ), dimension(1:3, 1:3) :: hgt   !! *local height 2D array*
   real   (kind=R8), intent(out), dimension(1:3)      :: gdx   !! *nodal x gradient*
   real   (kind=R8), intent(out), dimension(1:3)      :: gdy   !! *nodal x gradient*

!~       real(kind=R8), dimension(1:9) :: tab

!~       tab(1:9) = [ hgt(1, 1), & 1 !  (3,1)----(3,2)----(3,3)
!~                    hgt(1, 3), & 2 !    |        |        |
!~                    hgt(3, 3), & 3 !    |        |        |
!~                    hgt(3, 1), & 4 !  (2,1)----(2,2)----(2,3)
!~                    hgt(1, 2), & 5 !    |        |        |
!~                    hgt(2, 3), & 6 !    |        |        |
!~                    hgt(3, 2), & 7 !  (1,1)----(1,2)----(1,3)
!~                    hgt(2, 1), & 8 !
!~                    hgt(2, 2) ]  9 !

     ! (3,1)----(3,2)----(3,3)
     !   |        |        |
     !   |        |        |
     ! (2,1)----(2,2)----(2,3)
     !   |        |        |
     !   |        |        |
     ! (1,1)----(1,2)----(1,3)

      gdx(1:3) = [ -1.5_R8 * hgt(2, 1), &  !
                   +0.0_R8 * hgt(2, 2), &  !
                   +1.5_R8 * hgt(2, 3) ]   !

      gdy(1:3) = [ -1.5_R8 * hgt(1, 2), &  !
                   +0.0_R8 * hgt(2, 2), &  !
                   +1.5_R8 * hgt(3, 2) ]   !


!~       gdx(1:9) = [ -1.5_R8 * tab(1), &  !
!~                    +1.5_R8 * tab(2), &  !
!~                    +1.5_R8 * tab(3), &  !
!~                    -1.5_R8 * tab(4), &  !
!~                    +0.0_R8 * tab(5), &  !
!~                    +1.5_R8 * tab(6), &  !
!~                    +0.0_R8 * tab(7), &  !
!~                    -1.5_R8 * tab(8), &  !
!~                    +0.0_R8 * tab(9) ]   !


!~       gdy(1:9) = [ -1.5_R8 * tab(1), &  !
!~                    -1.5_R8 * tab(2), &  !
!~                    +1.5_R8 * tab(3), &  !
!~                    +1.5_R8 * tab(4), &  !
!~                    -1.5_R8 * tab(5), &  !
!~                    +0.0_R8 * tab(6), &  !
!~                    +1.5_R8 * tab(7), &  !
!~                    +0.0_R8 * tab(8), &  !
!~                    +0.0_R8 * tab(9) ]   !

   return
   endsubroutine gradient_corner


   subroutine labelize_point(height, label, x, y)
   !================================================================================================
   !! Function to label a point as: peak, valley, saddle or nothing particular
   implicit none
   real   (kind=R8), intent(in ), dimension(1:3, 1:3) :: height      !! *nodal height values as a 2D array*
   character(len=1), intent(out)                      :: label       !! *kind of point*
   real   (kind=R8), intent(out), optional            :: x, y        !! *coordinates of the extremum found*

      real(kind=R8), dimension(1:9)      :: tab
      real(kind=R8), dimension(1:9, 1:6) :: derivatives

      real(kind=R8) :: h, dhdx, dhdy, d2hdx2, d2hdy2, d2hdxdy, delta, dx, dy
      real(kind=R8) :: xs, ys

      real(kind=R8), parameter :: eps = 1.0e-8_R8

      integer(kind=I4) :: k, status

      ! following QU9 node notation:
      tab(1:9) = [ height(1, 1), &  !  (3,1)----(3,2)----(3,3)
                   height(1, 3), &  !    |        |        |
                   height(3, 3), &  !    |        |        |
                   height(3, 1), &  !  (2,1)----(2,2)----(2,3)
                   height(1, 2), &  !    |        |        |
                   height(2, 3), &  !    |        |        |
                   height(3, 2), &  !  (1,1)----(1,2)----(1,3)
                   height(2, 1), &  !
                   height(2, 2) ]   !

      ! Newton Raphson to locate the null first derivatives, which means an extremum
      ! Intiate with the center of the square [-1,1]X[-1,1] and initalize the counting
      !
      xs = 0 ; ys = 0 ; k = 0
      do

         ! for the current coordinates, what are the surface height and its derivatives:
         call deriv_N(x = xs, y = ys, mat_d = derivatives(1:9, 1:6))

         dhdx    = sum( derivatives(1:9, 2) * tab(1:9) )
         dhdy    = sum( derivatives(1:9, 3) * tab(1:9) )
         d2hdx2  = sum( derivatives(1:9, 4) * tab(1:9) )
         d2hdxdy = sum( derivatives(1:9, 5) * tab(1:9) )
         d2hdy2  = sum( derivatives(1:9, 6) * tab(1:9) )

         ! jacobian denominator
         delta = d2hdx2 * d2hdy2 - d2hdxdy**2

         if ( abs(dhdx) < eps .and. abs(dhdy) < eps ) then             ! converge ok
            status = 0                                                 ! extremum found, whatever its kind
            exit                                                       ! nothing more to do
         endif

         if ( abs(xs) >= 5._R8 .or. abs(ys) >= 5._R8 ) then            ! during the convergence process, the point is far from the square, so exit
            status = 1                                                 ! with the appropriate status
            exit                                                       !
         endif                                                         !

         k = k + 1

         if (k > 1000) then   ! limit the number of iterations
            status = 1
            exit
         endif

         dx = (-1. / delta) * ( + d2hdy2  * dhdx - d2hdxdy * dhdy)
         dy = (-1. / delta) * ( - d2hdxdy * dhdx + d2hdx2  * dhdy)

         xs = xs + 0.9 * dx
         ys = ys + 0.9 * dy


      enddo

      ! outside the square [-0.5,0.5]X[-0.5,0.5] the extremum belongs to another node
      if ( abs(xs) > 0.5_R8 .or. abs(ys) > 0.5_R8 ) status = 1

      ! if derivatives are null, what kind of point is it ?
      !     Nodes order:
      !
      !     4---7---3
      !     |   |   |
      !     8---9-- 6
      !     |   |   |
      !     1---5---2

      if (status == 0) then

         call deriv_N(x = xs, y = ys, mat_d = derivatives(1:9, 1:6))

         ! height at the extremum point xs, ys
         h = sum( derivatives(1:9, 1) * tab(1:9) )

         if (     all(tab(1:8) >= h) ) then ; label = 'V' ! valley or pit
         elseif ( all(tab(1:8) <= h) ) then ; label = 'P' ! hill or peak
         else ;                               label = 'S' ! saddle
         endif

      else

         label = 'N' ! nothing particular

      endif

      if ( present(x) ) then
         x = xs
         y = ys
      endif

   return
   endsubroutine labelize_point


   subroutine label_surf_summits(tab, nx, ny, valleys, peaks, saddles, nb_summits)
   !================================================================================================
   !! Function to output the extrema of a 2D array, as peaks, valleys or saddles.
   implicit none
   integer(kind=I4), intent(in ) :: nx                                     !! *number of pixels along x*
   integer(kind=I4), intent(in ) :: ny                                     !! *number of pixels along y*
   real   (kind=R8), intent(in ), dimension(1:nx, 1:ny) :: tab             !! *input 2D array*
   integer(kind=I4), intent(out), dimension(1:3)        :: nb_summits      !! *number of extrema of each kind*
   integer(kind=I4), intent(out), dimension(:,:), allocatable :: valleys   !! *list of valley coordinates*
   integer(kind=I4), intent(out), dimension(:,:), allocatable :: peaks     !! *list of peaks coordinates*
   integer(kind=I4), intent(out), dimension(:,:), allocatable :: saddles   !! *list of saddles coordinates*

      integer(kind=I4) :: i, j
      integer(kind=I4) :: ip ! peak counter
      integer(kind=I4) :: iv ! valley counter
      integer(kind=I4) :: is ! saddle counter

      character(len=1) :: label

      integer(kind=I4), dimension(1:3, 1:3) :: tp

      integer(kind=I4), dimension(:,:), allocatable :: topo ! point kind 2D array

      real(kind=R8), dimension(1:3, 1:3) :: ht

      real(kind=R8), dimension(1:3) :: gdx, gdy ! gradient vectors

      allocate( topo(1:nx, 1:ny) )

      topo(1:nx, 1:ny) = 0

      ! Loop through each point in the surface
      ip = 1 ; iv = 1 ; is = 1
      do i = 1 + 1, nx - 1

         do j = 1 + 1, ny - 1

            tp(1:3, 1:3) = topo(i-1:i+1, j-1:j+1)
            ht(1:3, 1:3) =  tab(i-1:i+1, j-1:j+1)

            ! if the gradients along x (resp. y) have the same sign, there is no extremum in the middle node
            call gradient_corner( hgt = ht(1:3, 1:3),   &  ! in
                                  gdx = gdx(1:3),       &  ! out
                                  gdy = gdy(1:3) )         ! out

!~             if ( all( [gdx(6) > 0, gdx(8) > 0, gdx(9) > 0] ) .or. all( [gdy(5) > 0, gdy(7) > 0, gdy(9) > 0] ) ) cycle
!~             if ( all( [gdx(6) > 0, gdx(8) > 0, gdx(9) > 0] ) .or. all( [gdy(5) < 0, gdy(7) < 0, gdy(9) < 0] ) ) cycle
!~             if ( all( [gdx(6) < 0, gdx(8) < 0, gdx(9) < 0] ) .or. all( [gdy(5) < 0, gdy(7) < 0, gdy(9) < 0] ) ) cycle
!~             if ( all( [gdx(6) < 0, gdx(8) < 0, gdx(9) < 0] ) .or. all( [gdy(5) > 0, gdy(7) > 0, gdy(9) > 0] ) ) cycle

            if ( all( gdx(1:3) > 0 ) .or. all( gdy(1:3) > 0 ) ) cycle
            if ( all( gdx(1:3) > 0 ) .or. all( gdy(1:3) < 0 ) ) cycle
            if ( all( gdx(1:3) < 0 ) .or. all( gdy(1:3) < 0 ) ) cycle
            if ( all( gdx(1:3) < 0 ) .or. all( gdy(1:3) > 0 ) ) cycle

            ! condition to avoid summits glued to each other: if a summit has been found in the neighborhood, cycle.
            if ( any( [tp(1, 1) > 0, &  !
                       tp(1, 2) > 0, &  !
                       tp(1, 3) > 0, &  !
                       tp(2, 1) > 0] ) ) cycle

            call labelize_point( height = ht(1:3, 1:3), &  ! in
                                 label  = label )          ! out

            selectcase(label)
               case('V') ; topo(i, j) = 1 ; iv = iv + 1 ! one valley more detected
               case('S') ; topo(i, j) = 2 ; is = is + 1 ! one saddle more detected
               case('P') ; topo(i, j) = 3 ; ip = ip + 1 ! one peak   more detected
            endselect

         enddo

      enddo

      nb_summits = [iv - 1, is - 1, ip - 1]

      ! now the number of extrema is known
      allocate( valleys(1:nb_summits(1), 1:2) )
      allocate( saddles(1:nb_summits(2), 1:2) )
      allocate( peaks  (1:nb_summits(3), 1:2) )

      ip = 1 ; iv = 1 ; is = 1
      do i = 1 + 1, nx - 1

         do j = 1 + 1, ny - 1

            selectcase( topo(i, j) )
               case(1) ; valleys(iv, 1) = i ;  valleys(iv, 2) = j ; iv = iv + 1
               case(2) ; saddles(is, 1) = i ;  saddles(is, 2) = j ; is = is + 1
               case(3) ; peaks  (ip, 1) = i ;  peaks  (ip, 2) = j ; ip = ip + 1
            endselect

         enddo

      enddo

      deallocate( topo )

   return
   endsubroutine label_surf_summits


   subroutine peaks_and_pits_curvatures(heights, nx, ny, dx, dy, S_param_grad, S_param_curv, peak_curv, pits_curv)
   !================================================================================================
   !! Function to calculate and output the peaks and pits curvatures as well as then mean quadratic
   !!        gradient value and the mean quadratic curvature value.
   implicit none
   integer(kind=I4), intent(in ) :: nx                              !! *number of pixels along x*
   integer(kind=I4), intent(in ) :: ny                              !! *number of pixels along x*
   real   (kind=R8), intent(in ) :: dx                              !! *x lag*
   real   (kind=R8), intent(in ) :: dy                              !! *y lag*
   real   (kind=R8), intent(out) :: peak_curv                       !! *3 first peaks mean curvature*
   real   (kind=R8), intent(out) :: pits_curv                       !! *3 first pits  mean curvature*
   real   (kind=R8), intent(out) :: S_param_grad                    !! *mean quadratic gradient value*
   real   (kind=R8), intent(out) :: S_param_curv                    !! *mean quadratic curvature value*
   real   (kind=R8), intent(in ), dimension(1:nx, 1:ny) :: heights  !! *input 2D array*

      integer(kind=I4) :: i, npeak, npits

      real(kind=R8) :: spg, spc, adim

      integer(kind=I4), allocatable, dimension(:,:) :: vall, peak, sadd

      integer(kind=I4), dimension(1:3) :: nb_extr

      real(kind=R8), dimension(1:nx*ny) :: tpits
      real(kind=R8), dimension(1:nx*ny) :: tpeak

      real(kind=R8), allocatable, dimension(:,:) :: cvt

      allocate( cvt(1:nx, 1:ny) )

      ! first determine the surface curvature
      call curvature(tab          = heights(1:nx, 1:ny),      & ! in
                     nx           = nx,                       & ! in
                     ny           = ny,                       & ! in
                     dx           = dx,                       & ! in
                     dy           = dy,                       & ! in
                     S_param_grad = spg,                      & ! out
                     S_param_curv = spc,                      & ! out
                     gcurvt       = cvt(1:nx, 1:ny))            ! out

      ! OUTPUT
      S_param_grad = spg
      S_param_curv = spc

      ! no need to carry very high/low values, so normalize curvature
      adim = maxval( abs(cvt(1:nx, 1:ny)) )

      cvt(1:nx, 1:ny) = cvt(1:nx, 1:ny) / adim

      call label_surf_summits(tab        = cvt(1:nx, 1:ny),     & ! in
                              nx         = nx,                  & ! in
                              ny         = ny,                  & ! in
                              valleys    = vall,                & ! out
                              peaks      = peak,                & ! out
                              saddles    = sadd,                & ! out
                              nb_summits = nb_extr(1:3))          ! out

      !.....................................................

      npits = nb_extr(1)

      do i = 1, npits

         tpits(i) = cvt( vall(i, 1), vall(i, 2) )     ! first values needed, so ascending sort is OK

      enddo

      call sort_array2( tab_inout = tpits(1:npits), &  ! inout
                        n         = npits )            ! in

      ! OUTPUT
      pits_curv = adim * sum( tpits(1:3) ) / 3.       ! mean of first 3 values, with the right dimension

      !.....................................................

      npeak = nb_extr(3)

      do i = 1, npeak

         tpeak(i) = -cvt( peak(i, 1), peak(i, 2) )    ! top values of peak needed, so reverse for ascending sort

      enddo

      call sort_array2( tab_inout = tpeak(1:npeak), &  ! inout
                        n         = npeak )            ! in

      ! OUTPUT
      peak_curv = - adim * sum( tpeak(1:3) ) / 3      ! mean of first 3 values, with the right dimension and the right sign


      deallocate( cvt, vall, peak, sadd )

   return
   endsubroutine peaks_and_pits_curvatures


   subroutine test_labelize_point()
   !================================================================================================
   !! Function to test the function "labelize_point" on a QU9 domain with a 2nd order polynomial
   !!        along x and y
   implicit none

      real(kind=R8) :: x, y, x0, y0
      real(kind=R8) :: aax, aay

      character(len=1) :: point_kind, label

      real(kind=R8), parameter :: xx0 = -0.1_R8, fxx0 = +2._R8
      real(kind=R8), parameter :: yy0 = +0.2_R8, fyy0 = +1._R8

      real(kind=R8), dimension(1:9) :: Ni, dNidx, dNidy, d2Nidx2, d2Nidxdy, d2Nidy2, tab

      real(kind=R8), dimension(1:9, 1:6) :: derivatives
      real(kind=R8), dimension(1:3, 1:3) :: height

      ! select what kind of point to check
      point_kind = 'S'

      selectcase( point_kind )
         case( 'V' )  ; aax = +1._R8 ; aay = +1._R8
         case( 'P' )  ; aax = -1._R8 ; aay = -1._R8
         case( 'S' )  ; aax = +1._R8 ; aay = -1._R8
         case default ; stop
      endselect

      ! "tab" is given row wise:
      !     Nodes order:
      !
      !     4---7---3
      !     |   |   |
      !     8---9-- 6
      !     |   |   |
      !     1---5---2
      !
      tab(1:9) = [ f(x = -1._R8, x0 = xx0, a = aax, fx0 = fxx0) * f(x = -1._R8, x0 = yy0, a = aay, fx0 = fyy0),   &  ! 1
                   f(x = +1._R8, x0 = xx0, a = aax, fx0 = fxx0) * f(x = -1._R8, x0 = yy0, a = aay, fx0 = fyy0),   &  ! 5
                   f(x = +1._R8, x0 = xx0, a = aax, fx0 = fxx0) * f(x = +1._R8, x0 = yy0, a = aay, fx0 = fyy0),   &  ! 2
                   f(x = -1._R8, x0 = xx0, a = aax, fx0 = fxx0) * f(x = +1._R8, x0 = yy0, a = aay, fx0 = fyy0),   &  ! 8
                   f(x =  0._R8, x0 = xx0, a = aax, fx0 = fxx0) * f(x = -1._R8, x0 = yy0, a = aay, fx0 = fyy0),   &  ! 9
                   f(x = +1._R8, x0 = xx0, a = aax, fx0 = fxx0) * f(x =  0._R8, x0 = yy0, a = aay, fx0 = fyy0),   &  ! 6
                   f(x =  0._R8, x0 = xx0, a = aax, fx0 = fxx0) * f(x = +1._R8, x0 = yy0, a = aay, fx0 = fyy0),   &  ! 4
                   f(x = -1._R8, x0 = xx0, a = aax, fx0 = fxx0) * f(x =  0._R8, x0 = yy0, a = aay, fx0 = fyy0),   &  ! 7
                   f(x =  0._R8, x0 = xx0, a = aax, fx0 = fxx0) * f(x =  0._R8, x0 = yy0, a = aay, fx0 = fyy0) ]     ! 3

      ! stack is done column wise for reshape
      height(1:3, 1:3) = reshape( [ tab(1),   &  !
                                    tab(8),   &  !
                                    tab(4),   &  !
                                    tab(5),   &  !
                                    tab(9),   &  !
                                    tab(7),   &  !
                                    tab(2),   &  !
                                    tab(6),   &  !
                                    tab(3) ], [3, 3] )

      ! random point to compare the analytic function to the quadratic approx
      x = -0.35_R8
      y = +0.52_R8

      call deriv_N(x = x, y = y, mat_d = derivatives(1:9, 1:6))

      Ni(1:9)       = derivatives(1:9, 1)
      dNidx(1:9)    = derivatives(1:9, 2)
      dNidy(1:9)    = derivatives(1:9, 3)
      d2Nidx2(1:9)  = derivatives(1:9, 4)
      d2Nidxdy(1:9) = derivatives(1:9, 5)
      d2Nidy2(1:9)  = derivatives(1:9, 6)

      write(*,*) 'numerical h      : ', sum(tab * Ni)        ,  ' ; theoretical h      : ',   f(x = x, x0 = xx0, a = aax, fx0 = fxx0) *   f(x = y, x0 = yy0, a = aay, fx0 = fyy0)
      write(*,*) 'numerical dhdx   : ', sum(tab * dNidx)     ,  ' ; theoretical dhdx   : ',  df(x = x, x0 = xx0, a = aax            ) *   f(x = y, x0 = yy0, a = aay, fx0 = fyy0)
      write(*,*) 'numerical dhdy   : ', sum(tab * dNidy)     ,  ' ; theoretical dhdy   : ',   f(x = x, x0 = xx0, a = aax, fx0 = fxx0) *  df(x = y, x0 = yy0, a = aay            )
      write(*,*) 'numerical d2hdx2 : ', sum(tab * d2Nidx2)   ,  ' ; theoretical d2hdx2 : ', d2f(                 a = aax            ) *   f(x = y, x0 = yy0, a = aay, fx0 = fyy0)
      write(*,*) 'numerical d2hdy2 : ', sum(tab * d2Nidy2)   ,  ' ; theoretical d2hdy2 : ',   f(x = x, x0 = xx0, a = aax, fx0 = fxx0) * d2f(                 a = aay            )
      write(*,*) 'numerical d2hdxdy: ', sum(tab * d2Nidxdy)  ,  ' ; theoretical d2hdxdy: ',  df(x = x, x0 = xx0, a = aax            ) *  df(x = y, x0 = yy0, a = aay            )

      !===========================================

      call labelize_point(height = height(1:3, 1:3), label = label, x = x0, y = y0)

      write(*,*) 'theoretical xx0: ', xx0, 'numerical xx0: ', x0
      write(*,*) 'theoretical yy0: ', yy0, 'numerical yy0: ', y0

      write(*,*) 'theoretical point: ', point_kind, ' numerical point: ', label

   contains

      real(kind=R8) function f(x, x0, a, fx0)
      implicit none
      real(kind=R8), intent(in) :: x
      real(kind=R8), intent(in) :: x0
      real(kind=R8), intent(in) :: a
      real(kind=R8), intent(in) :: fx0

         real(kind=R8) :: b, c

         ! 2nd order polynomial defined by:
         ! ->   a : +1 or -1 (curvature sign)
         ! ->  x0 : the extremum abscissa
         ! -> fx0 : the extremum value
         b = -2 * a * x0
         c = fx0 + a * x0**2

         f = a * x**2 + b * x + c

      return
      endfunction f

      real(kind=R8) function df(x, x0, a)
      implicit none
      real(kind=R8), intent(in) :: x
      real(kind=R8), intent(in) :: x0
      real(kind=R8), intent(in) :: a

         real(kind=R8) :: b

         b = -2 * a * x0

         df = 2 * a * x + b

      return
      endfunction df

      real(kind=R8) function d2f(a)
      implicit none
      real(kind=R8), intent(in) :: a

         d2f = 2 * a

      return
      endfunction d2f

   endsubroutine test_labelize_point


   subroutine test_label_surf_summits()
   !================================================================================================
   !! Function to test the capicity in detecting peaks, pits and saddles in a simple double sinus
   !!        surface.
   implicit none

      integer(kind=I4), allocatable, dimension(:,:) :: topo
      integer(kind=I4), allocatable, dimension(:,:) :: vall, peak, sadd

      real   (kind=R8), allocatable, dimension(:,:) :: heights

      integer(kind=I4), dimension(1:3) :: nb_null_derivatives

      integer(kind=I4), parameter :: nx = 750, ny = 1000

      integer(kind=I4) :: i, j, i1, j1

      type(SCALE_SURF) :: scal

      ! create a "digital surf" object
      call init_scal(scal = scal,          &      ! out; creates a surface type, containing ...
                       nx = nx,            &      !  in; ... the number of points along x ...
                       ny = ny,            &      !  in; ... the number of points along y ...
                       lx = nx * 1.e-6_R8, &      !  in; ... the length (default unit : m) ...
                       ly = ny * 1.e-6_R8, &      !  in; ... the width ...
                   unit_z = 'm'       )           !  in; ... and the unit along z.

      allocate( heights(1:nx, 1:ny) )
      allocate( topo   (1:nx, 1:ny) )

      do i = 1, nx
      do j = 1, ny

         heights(i, j) = sinsin(i, j, nx, ny)

      enddo
      enddo

      call write_surf( nom_fic = "out/test_sinus.sur",     &  !
                       tab_s   = heights(1:nx, 1:ny),      &  !
                       scal    = scal )                       !

      call label_surf_summits(tab        = heights(1:nx, 1:ny),      &  !
                              nx         = nx,                       &  !
                              ny         = ny,                       &  !
                              valleys    = vall,                     &  !
                              peaks      = peak,                     &  !
                              saddles    = sadd,                     &  !
                              nb_summits = nb_null_derivatives(1:3) )   !

      write(*,*) nb_null_derivatives(1:3)

      ! csv to compare computed and theoretical values
      open(10, file = 'out/extrema.csv')

         do i = 1, nb_null_derivatives(1)
            write(10,*) 'comput,', vall(i, 1), ',', vall(i, 2), ',', '1,VALLEY'
         enddo

         do i = 1, nb_null_derivatives(2)
            write(10,*) 'comput,', sadd(i, 1), ',', sadd(i, 2), ',', '2,SADDLE'
         enddo

         do i = 1, nb_null_derivatives(3)
            write(10,*) 'comput,', peak(i, 1), ',', peak(i, 2), ',', '3,PEAK'
         enddo

         do i = 0, 2
         i1 = int( nx * (0.5 + i) / 3. )

            do j = 0, 5
            j1 = int( ny * (0.5 + j) / 6.)

               if ( (j == 2 * (j/2)     .and. i == 2 * (i/2)    ) .or.     &  !
                    (j == 2 * (j/2) + 1 .and. i == 2 * (i/2) + 1) ) then      !

                  write(10,*) 'theory,', i1, ',', j1, ',', '3,PEAK'

               else

                  write(10,*) 'theory,', i1, ',', j1, ',', '1,VALLEY'

               endif

         enddo

      enddo

      do i = 0, 2
      i1 = int( nx * (1.0 + i) / 3. )

         do j = 0, 5
         j1 = int( ny * (1.0 + j) / 6.)

            if ( j1 /= ny .and. i1 /= nx ) write(10,*) 'theory,', i1, ',', j1, ',', '2,SADDLE'

         enddo

      enddo

      close(10)

      deallocate( heights, topo, vall, peak, sadd)

   contains

      real(kind=R8) function sinsin(i, j, nx, ny)
      implicit none
      integer(kind=I4), intent(in) :: i, j, nx, ny

         sinsin = sin( 3 * PI_R8 * i / nx + 1.e-8_R8 ) * sin( 6 * PI_R8 * j / ny + 2.e-8_R8 )

      return
      endfunction sinsin

   endsubroutine test_label_surf_summits


   subroutine test_peaks_and_pits_curvatures()
   !================================================================================================
   !! Function to test the function "peaks_and_pits_curvatures" on a real rough surface.
   !!        Outputs surface gradients and curvatures as 2D arrays or single values.
   implicit none

      real(kind=R8), allocatable, dimension(:,:) :: heights, bf_heights
      real(kind=R8), allocatable, dimension(:,:) :: gradx, grady, gradxx, gradxy, gradyy, cvt

      real(kind=R8)  :: dx, dy, dz, pcurv, pgrad, pmin, pmax, fft_cutoff

      integer(kind=I4) :: nx, ny, nx2, ny2, n_th

      type(SCALE_SURF) :: scal_surf

      call read_surf( nom_fic = "sur/test1.sur", &  ! IN
                        tab_s = heights,         &  ! OUT
                         scal = scal_surf )         ! OUT

      nx = scal_surf%xres
      ny = scal_surf%yres

      dx = scal_surf%dx * unit2IUf(scal_surf%dx_unit)
      dy = scal_surf%dy * unit2IUf(scal_surf%dy_unit)
      dz = 0

      fft_cutoff = dx / 5.e-6

      nx2 = 2 * ( nint(PAD_FFT * nx)/2 )
      ny2 = 2 * ( nint(PAD_FFT * ny)/2 )

      write(*,*) 'nx, ny = ', nx, ny
      write(*,*) 'dx, dy = ', dx, dy

      !=================================================================

      allocate( bf_heights(1:nx, 1:ny) )

      allocate( gradx(1:nx, 1:ny) )
      allocate( grady(1:nx, 1:ny) )

      allocate( gradxx(1:nx, 1:ny) )
      allocate( gradxy(1:nx, 1:ny) )
      allocate( gradyy(1:nx, 1:ny) )

      allocate( cvt(1:nx, 1:ny) )

      n_th = omp_get_max_threads()
      call fftw_plan_with_nthreads(nthreads = n_th)

      call init_fftw3(long = nx2,  & !
                      larg = ny2 )   !

      call fft_filter(tab       = heights(1:nx, 1:ny),      & ! in
                      long      = nx,                       & ! in
                      larg      = ny,                       & ! in
                      cutoff    = fft_cutoff,               & ! in
                      bf_tab    = bf_heights(1:nx, 1:ny),   & ! out
                      multi_fft = .FALSE.)                    ! in

      call end_fftw3()

      call write_surf( nom_fic = "out/test_bf_heights.sur", &  ! in
                       tab_s   = bf_heights(1:nx, 1:ny),    &  ! in
                       scal    = scal_surf )                   ! in

      call gradient(tab   = bf_heights(1:nx, 1:ny),   & ! in
                    nx    = nx,                       & ! in
                    ny    = ny,                       & ! in
                    dx    = dx,                       & ! in
                    dy    = dy,                       & ! in
                    gradx = gradx(1:nx, 1:ny),        & ! out
                    grady = grady(1:nx, 1:ny))          ! out

            call write_surf( nom_fic = "out/test_gradq.sur",  &  ! in
                             tab_s   = gradx(1:nx, 1:ny)**2 + &  ! in
                                       grady(1:nx, 1:ny)**2,  &  !
                             scal    = scal_surf )               ! in

            call write_surf( nom_fic = "out/test_gradx.sur",  &  ! in
                             tab_s   = gradx(1:nx, 1:ny),     &  ! in
                             scal    = scal_surf )               ! in

            call write_surf( nom_fic = "out/test_grady.sur",  &  ! in
                             tab_s   = grady(1:nx, 1:ny),     &  ! in
                             scal    = scal_surf )               ! in

      !=================================================================

      call gauss_curv(gradx  = gradx (1:nx, 1:ny),    & !
                      grady  = grady (1:nx, 1:ny),    & !
                      gradxx = gradxx(1:nx, 1:ny),    & !
                      gradxy = gradxy(1:nx, 1:ny),    & !
                      gradyy = gradyy(1:nx, 1:ny),    & !
                      nx     = nx,                    & !
                      ny     = ny,                    & !
                      dx     = dx,                    & !
                      dy     = dy)                      !

            call write_surf( nom_fic = "out/test_gradxx_1.sur",  &  ! in
                             tab_s   = gradxx(1:nx, 1:ny),       &  ! in
                             scal    = scal_surf )                  ! in

            call write_surf( nom_fic = "out/test_gradxy_1.sur",  &  ! in
                             tab_s   = gradxy(1:nx, 1:ny),       &  ! in
                             scal    = scal_surf )                  ! in

            call write_surf( nom_fic = "out/test_gradyy_1.sur",  &  ! in
                             tab_s   = gradyy(1:nx, 1:ny),       &  ! in
                             scal    = scal_surf )                  ! in

      call curv2(tab    = bf_heights(1:nx, 1:ny),    & ! in
                 nx     = nx,                        & ! in
                 ny     = ny,                        & ! in
                 dx     = dx,                        & ! in
                 dy     = dy,                        & ! in
                 gradxx = gradxx    (1:nx, 1:ny),    & ! out
                 gradyy = gradyy    (1:nx, 1:ny))      ! out

            call write_surf( nom_fic = "out/test_gradxx_2.sur",  &  ! in
                             tab_s   = gradxx(1:nx, 1:ny),       &  ! in
                             scal    = scal_surf )                  ! in

            call write_surf( nom_fic = "out/test_gradxy_2.sur",  &  ! in
                             tab_s   = gradxy(1:nx, 1:ny),       &  ! in
                             scal    = scal_surf )                  ! in

            call write_surf( nom_fic = "out/test_gradyy_2.sur",  &  ! in
                             tab_s   = gradyy(1:nx, 1:ny),       &  ! in
                             scal    = scal_surf )                  ! in

      call curvature(tab          = bf_heights(1:nx, 1:ny),   & ! in
                     nx           = nx,                       & ! in
                     ny           = ny,                       & ! in
                     dx           = dx,                       & ! in
                     dy           = dy,                       & ! in
                     S_param_grad = pgrad,                    & ! out
                     S_param_curv = pcurv,                    & ! out
                     gcurvt       = cvt(1:nx, 1:ny))            ! out

      write(*,*) "S_param_grad: ", pgrad

            call write_surf( nom_fic = "out/test_cvt.sur",       &  ! in
                             tab_s   = cvt(1:nx, 1:ny),          &  ! in
                             scal    = scal_surf )                  ! in

      !=================================================================

      call peaks_and_pits_curvatures( heights      = bf_heights(1:nx, 1:ny),   & ! in
                                      nx           = nx,                       & ! in
                                      ny           = ny,                       & ! in
                                      dx           = dx,                       & ! in
                                      dy           = dy,                       & ! in
                                      S_param_grad = pgrad,                    & ! out
                                      S_param_curv = pcurv,                    & ! out
                                      peak_curv    = pmax,                     & ! out
                                      pits_curv    = pmin )                      ! out

      write(*,*) "S_param_curv: ", pcurv
      write(*,*) "peak_curv: ",    pmin
      write(*,*) "pits_curv: ",    pmax

      !=================================================================

      deallocate( heights, bf_heights, gradx, grady, gradxy, gradxx, gradyy, cvt )

   return
   endsubroutine test_peaks_and_pits_curvatures

endmodule grad_curv
