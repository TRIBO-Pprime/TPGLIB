!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: april, 06 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Morphological operations.  Example of use.**
!<  </span>
program test_morpho
use data_arch,     only : I4, R8, PI_R8
use miscellaneous, only : get_unit
use surfile,       only : init_scal, read_surf, write_surf, SCALE_SURF
use sort_arrays,   only : sort_array2
use stat_mom,      only : calc_moments, moment_stat
use morpho,        only : calcul_normales, surf_area, count_cell, erode_dilate, def_masque
implicit none

real   (kind=R8), allocatable, dimension(:,:)  :: mask, array
real   (kind=R8), allocatable, dimension(:)    :: vec, vec_tmp
integer(kind=I4), allocatable, dimension(:,:)  :: imask

type(SCALE_SURF) :: scal_surf

integer(kind=I4) :: i, nnx, nny, nb_cells
real   (kind=R8) :: median_size, c1, c2, topo, mintab, maxtab

integer(kind=I4)  :: ni, nj, ii, jj, j, k, ua
real   (kind=R8)  :: angle, dx, dy, t, cp, sp, res_hori, area
type(moment_stat) :: mom

!----------------------------------------------------------------------------------
!---------- Returns the number of cells of a given mask ---------------------------
!----------------------------------------------------------------------------------

call read_surf( nom_fic = "sur/mask.sur",    &  !
                tab_s   = mask,              &  !
                scal    = scal_surf )           !

nnx = scal_surf%xres
nny = scal_surf%yres

allocate( imask(1:nnx, 1:nny) )

imask(1:nnx, 1:nny) = nint( mask(1:nnx, 1:nny) )

call count_cell(  msk      = imask(1:nnx, 1:nny),  &  !
                  long     = nnx,                  &  !
                  larg     = nny,                  &  !
                  nbr_cell = nb_cells,             &  !
                  med_cell = median_size )            !

call write_surf( nom_fic = "out/res_count_cell.sur",                 &  !
                 tab_s   = real( imask(1:nnx, 1:nny), kind = R8 ),   &  !
                 scal    = scal_surf )                                  !

write(*, '(a, T25, I3, a, E12.4, a)') 'initial mask:', nb_cells, ' cells, ', median_size, ' % surface'

!----------------------------------------------------------------------------------
!---- Returns the number of cells of a given mask after erosion -------------------
!----------------------------------------------------------------------------------

imask(1:nnx, 1:nny) = nint( mask(1:nnx, 1:nny) )

call erode_dilate( msk  = imask(1:nnx, 1:nny),  &  !
                   long = nnx,                  &  !
                   larg = nny,                  &  !
                   val  = 6,                    &  !
                   act  = 'erode' )                !

call count_cell(  msk      = imask(1:nnx, 1:nny),  &  !
                  long     = nnx,                  &  !
                  larg     = nny,                  &  !
                  nbr_cell = nb_cells,             &  !
                  med_cell = median_size )            !

call write_surf( nom_fic = "out/res_count_cell_erode.sur",           &  !
                 tab_s   = real( imask(1:nnx, 1:nny), kind = R8 ),   &  !
                 scal    = scal_surf )                                  !

write(*, '(a, T25, I3, a, E12.4, a)') 'mask after erode:', nb_cells, ' cells, ', median_size, ' % surface'

!----------------------------------------------------------------------------------
!----- Returns the number of cells of a given mask after erosion and dilation -----
!----------------------------------------------------------------------------------

where ( imask(1:nnx, 1:nny) >= 1 ) imask(1:nnx, 1:nny) = 1

call erode_dilate( msk  = imask(1:nnx, 1:nny),  &  !
                   long = nnx,                  &  !
                   larg = nny,                  &  !
                   val  = 6,                    &  !
                   act  = 'dilat' )                !

call count_cell(  msk      = imask(1:nnx, 1:nny),  &  !
                  long     = nnx,                  &  !
                  larg     = nny,                  &  !
                  nbr_cell = nb_cells,             &  !
                  med_cell = median_size )            !

call write_surf( nom_fic = "out/res_count_cell_erode_dilat.sur",     &  !
                 tab_s   = real( imask(1:nnx, 1:nny), kind = R8 ),   &  !
                 scal    = scal_surf )                                  !

write(*, '(a, T25, I3, a, E12.4, a)') 'mask erode + dilate:', nb_cells, ' cells, ', median_size, ' % surface'
write(*,*) '--------------------------------------------------------'

!----------------------------------------------------------------------------------
!---------- Returns the % age of heights above c2, once c1 heights are removed ----
!----------------------------------------------------------------------------------

allocate( array( 1:nnx, 1:nny ) )
allocate(   vec( 1:nnx  * nny ), vec_tmp( 1:nnx  * nny ) )

! thresholds
c1 = 0.15_R8
c2 = 0.85_R8

! surface made of increasing integers
vec( 1:nnx  * nny ) = [ ( i, i = 1, nnx * nny ) ]

mintab = minval( vec(1:nnx*nny) )
maxtab = maxval( vec(1:nnx*nny) )

! heights rescaled between 0 and 1
vec_tmp(1:nnx*nny) = ( vec(1:nnx*nny) - mintab )/( maxtab - mintab )

! shuffle heights
call random_number( vec(1:nnx*nny) )

call sort_array2(tab_inout =     vec(1:nnx*nny),   &  !
                      tab1 = vec_tmp(1:nnx*nny),   &  !
                         n = nnx*nny)                 !

array( 1:nnx, 1:nny ) = reshape( vec_tmp(1:nnx*nny), [nnx, nny] )

! mask heights above c2, when heights below c1 are ignored
call def_masque( msk = imask(1:nnx, 1:nny),  &  !
                 tab = array(1:nnx, 1:nny),  &  !
                long = nnx ,                 &  !
                larg = nny ,                 &  !
               crit1 = c1,                   &  !
               crit2 = c2,                   &  !
                 top = topo )                   !

call write_surf( nom_fic = "out/mask_after_threshold.sur",           &  !
                 tab_s   = real( imask(1:nnx, 1:nny), kind = R8 ),   &  !
                 scal    = scal_surf )                                  !

! topo is the fraction of heights that are involved
write(*,*) 'percentage of heights abobe 85%, without the first 15%: ', topo ! = 100. - 100 * ( c1 + c2 * (1. - c1) )
write(*,*) '--------------------------------------------------------'

deallocate( array, mask, imask, vec, vec_tmp )

!----------------------------------------------------------------------------------
!--------- Returns the % age of surface nearly horizontal -------------------------
!----------------------------------------------------------------------------------

nnx = 512 ; nny = 512

allocate( array( 1:nnx, 1:nny ) )
allocate(   vec( 1:nnx * nny )  )

! the domain is divided into 4x4 patches
ni = 512 / 4
nj = 512 / 4

! cone half angle
angle = 5._R8

! lag along x and y
dx    = 1.e-7_R8
dy    = 1.e-7_R8

k = 0
do i = 1, 4
   do j = 1, 4

      k = k + 1 ! subdomain number

      ! the subdomain is a plane of equation: z = -a.x -b.y -c
      ! the normal coordinates are (-a, -b, +1)/sqrt(a**2 + b**2 + 1)
      ! and the scalar product with \vec{k} is 1/sqrt(a**2 + b**2 + 1)
      ! which is also cos(\theta).
      ! Therefore a**2 + b**2 = tan(\theta)**2, so if a = cos(\phi).tan(\theta)
      ! b = sin(\phi).tan(\theta), the relationship is satisfied (whatever \phi).
      ! With \theta = k * 0.11 * angle, 9 subdomains are concerned.

      t  = tan( k * 0.11 * angle * PI_R8 / 180 )

      call random_number( cp )
      cp = 2 * (cp - 0.5)

      call random_number( sp )
      if (sp > 0.5_R8) then
         sp = + sqrt( 1._R8 - cp**2)
      else
         sp = - sqrt( 1._R8 - cp**2)
      endif

      do ii = (i - 1) * ni + 1, i * ni
      do jj = (j - 1) * nj + 1, j * nj

         array(ii, jj) = plane( a = cp * t,                 &  !
                                b = sp * t,                 &  !
                                c = 0._R8,                  &  !
                                x = ii * dx,                &  !
                                y = jj * dy )                  !

      enddo
      enddo

   enddo
enddo

call calc_moments( tab    = reshape( array(1:nnx, 1:nny), [nnx * nny] ),   &  !
                   mx     = mom,                                           &  !
                   nb_mom = 2 )                                               !

array(1:nnx, 1:nny) = ( array(1:nnx, 1:nny) - mom%mu ) / mom%si

call init_scal( scal   = scal_surf,    &  !
                nx     = nnx,          &  !
                ny     = nny,          &  !
                lx     = nnx * dx,     &  !
                ly     = nny * dy,     &  !
                unit_z = 'm ' )           !

call write_surf( nom_fic = "out/hori.sur",      &  !
                 tab_s   = array(1:nnx, 1:nny), &  !
                 scal    = scal_surf )             !

call calcul_normales( tab_in     = array(1:nnx, 1:nny),  &  !
                      long       = nnx,                  &  !
                      larg       = nny,                  &  !
                      scale_xyz  = [dx, dy, mom%si],     &  !
                      cone_angle = angle,                &  !
                      print_mask = .true.,               &  !
                      hori       = res_hori )               !

call surf_area( tab_in     = array(1:nnx, 1:nny),  &  !
                long       = nnx,                  &  !
                larg       = nny,                  &  !
                scale_xyz  = [dx, dy, mom%si],     &  !
                aire       = area )                   !

call get_unit(ua)
open( unit = ua, file = "out/mask_angle.txt")
   read( ua, * ) ( vec(i), i = 1, nnx * nny )
close( ua )

call write_surf( nom_fic = "out/mask_hori.sur",                   &  !
                 tab_s   = reshape( vec(1:nnx*nny), [nnx, nny] ), &  !
                 scal    = scal_surf )                               !

write(*, *) 'Percentage of nearly horizontal surface (+/- 5°): ', res_hori
write(*, *) 'Relativea area minus 1 and z scale:               ', area, mom%si

deallocate( array, vec )

contains

real(kind=R8) function plane(a, b, c, x, y)
implicit none
real(kind=R8), intent(in) :: a, b, c, x, y

   plane = a * x + b * y + c

return
endfunction plane

endprogram test_morpho
