!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: november, 01 2024
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<     **Firestone Abbott's curve. Example of use**
!<  </span>
program test_abbott
use data_arch,     only : I4, R8
use miscellaneous, only : get_unit
use surfile,       only : read_surf, write_surf, SCALE_SURF
use abbott,        only : abbott_param

implicit none

real(kind=R8), allocatable, dimension(:)   :: vec_heights
real(kind=R8), allocatable, dimension(:,:) :: heights

real(kind=R8), dimension(1:11) :: res

integer(kind=I4) :: nx, ny

type(SCALE_SURF) :: scal_surf

call read_surf( nom_fic = "sur/test.sur", &  ! IN
                  tab_s = heights,        &  ! OUT
                   scal = scal_surf )        ! OUT

nx = scal_surf%xres
ny = scal_surf%yres

allocate( vec_heights(1:nx * ny) )

vec_heights(1:nx * ny) = reshape( heights, [nx * ny] )

call abbott_param(   tab     = vec_heights(1:nx * ny),      &  !
                     lg      = nx * ny,                     &  !
                     nom     = 'out/test',                  &  !
                     curves  = [.true., .true., .true.],    &  !
                     results = res(1:11),                   &  !
                     omp     = .true. )                        !

write(*, *) 'smr1, iso 25178.................: ', res( 1)
write(*, *) 'smr2, iso 25178.................: ', res( 2)
write(*, *) 'spk , iso 25178.................: ', res( 3)
write(*, *) 'svk , iso 25178.................: ', res( 4)
write(*, *) 'off1, ordonnée de spk...........: ', res( 5)
write(*, *) 'off2, ordonnée de svk...........: ', res( 6)
write(*, *) 'sk  , iso 25178.................: ', res( 7)
write(*, *) 'core slope..................... : ', res( 8)
write(*, *) 'adjustment factor (tangent fit) : ', res( 9)
write(*, *) 'coeffa_tan        (tangent fit) : ', res(10)
write(*, *) 'coeffb_tan        (tangent fit) : ', res(11)

write(*, *)

write(*, *) 'NB: on a reduced Abbott curve, the tangent fit is not so good'   //    &  !
            ' because the beginning and the end of the data points are kept'  //    &  !
            ' making the curve too much sharp.'                                        !

deallocate( heights, vec_heights )

stop
endprogram test_abbott
