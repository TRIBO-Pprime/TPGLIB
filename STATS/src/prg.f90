!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: may, 03 2019
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Routines to calculate statistical moments. Example of use.**
!<  </span>
program test_stat
use data_arch, only : I4, R8, PI_R8
use stat_mom

implicit none

real(kind=R8),    allocatable, dimension(:) :: array
logical(kind=I4), allocatable, dimension(:) :: array_mask

real(kind=R8),    allocatable, dimension(:,:) :: array2D
logical(kind=I4), allocatable, dimension(:,:) :: array_mask2D

real(kind=R8), dimension(1:6) :: results

integer(kind=I4), parameter :: rn = 1e3
integer(kind=I4), parameter ::  n = rn**2

type(moment_stat) :: mom

integer(kind=I4) :: i
real(kind=R8)    :: x, median

real(kind=R8), allocatable, dimension(:) :: vx
real(kind=R8), parameter                 :: mu = 2.0_R8, sigma = 3.0_R8
integer(kind=I4), parameter              :: nn = 1e7
integer(kind=I4)                         :: ipow, iter

allocate(      array(1:n) )
allocate( array_mask(1:n) )

allocate(      array2D(1:rn, 1:rn) )
allocate( array_mask2D(1:rn, 1:rn) )

array_mask = .false.

do i = 0, n -1
   x = real(i, kind = R8) / (n - 1)
   array(i + 1) = sin( 4 * PI_R8 * x ) * exp( -5 * x ) + x**2

   if ( mod(i, 2) == 0 ) array_mask(i + 1) = .True.
enddo

array2D      = reshape(      array(1:n), [rn, rn] )
array_mask2D = reshape( array_mask(1:n), [rn, rn] )

!============================================================================== 1D WITH MASK

call calc_moments(tab = array, mask = array_mask, mx = mom, nb_mom = 4)
call calc_median(tab = array, mask = array_mask, md = median)

results = [ 0.4015711287531725,  0.08020113698321299,  0.28319805257666053,  0.07977645748158539,  2.1330678728336125,  0.42620050599180814 ] ! python
!         [ 0.40157112875317|66, 0.08020113698321|186, 0.2831980525766|5853, 0.0797764574815|4369, 2.1330678728336|414, 0.42620050599180814 ] ! fortran

write(*, *) '-----------------------------------------------------------------'
write(*, *) 'With mask :'
write(*, *) 'sum of abs error: ',   abs(mom%mu - results(1)) +    &  !
                                    abs(mom%va - results(2)) +    &  !
                                    abs(mom%si - results(3)) +    &  !
                                    abs(mom%Sk - results(4)) +    &  !
                                    abs(mom%Ku - results(5)) +    &  !
                                    abs(median - results(6))         !

write(*, *) 'detail: ', mom%mu, mom%va, mom%si, mom%Sk, mom%Ku, median

!============================================================================== 1D NO MASK

call calc_moments(tab = array, mx = mom, nb_mom = 4)
call calc_median(tab = array, md = median)

results = [ 0.4015716287568354,   0.08020123540993489,    0.28319822635379427,   0.07977715735614675,   2.1330700214753766,  0.42620146392518804 ]  ! python
!         [ 0.4015716287568|2495, 0.08020123540993|2986,  0.28319822635379|088,  0.079777157356|279285, 2.133070021475|4486, 0.426201463925188|10 ] ! fortran
write(*, *) '-----------------------------------------------------------------'
write(*, *) 'With NO mask :'
write(*, *) 'sum of abs error: ',   abs(mom%mu - results(1)) +    &  !
                                    abs(mom%va - results(2)) +    &  !
                                    abs(mom%si - results(3)) +    &  !
                                    abs(mom%Sk - results(4)) +    &  !
                                    abs(mom%Ku - results(5)) +    &  !
                                    abs(median - results(6))         !

write(*, *) 'detail: ', mom%mu, mom%va, mom%si, mom%Sk, mom%Ku, median

!============================================================================== 2D WITH MASK

call calc_moments(tab = array2D, mask = array_mask2D, mx = mom, nb_mom = 4)
call calc_median(tab = array, mask = array_mask, md = median)

results = [ 0.4015711287531725,  0.08020113698321299,  0.28319805257666053,  0.07977645748158539,  2.1330678728336125,  0.42620050599180814 ] ! python
!         [ 0.40157112875317|66, 0.08020113698321|186, 0.2831980525766|5853, 0.0797764574815|4369, 2.1330678728336|414, 0.42620050599180814 ] ! fortran

write(*, *) '-----------------------------------------------------------------'
write(*, *) 'With mask :'
write(*, *) 'sum of abs error: ',   abs(mom%mu - results(1)) +    &  !
                                    abs(mom%va - results(2)) +    &  !
                                    abs(mom%si - results(3)) +    &  !
                                    abs(mom%Sk - results(4)) +    &  !
                                    abs(mom%Ku - results(5)) +    &  !
                                    abs(median - results(6))         !

write(*, *) 'detail: ', mom%mu, mom%va, mom%si, mom%Sk, mom%Ku, median

!============================================================================== 2D NO MASK

call calc_moments(tab = array2D, mx = mom, nb_mom = 4)
call calc_median(tab = array, md = median)

results = [ 0.4015716287568354,   0.08020123540993489,    0.28319822635379427,   0.07977715735614675,   2.1330700214753766,  0.42620146392518804 ]  ! python
!         [ 0.4015716287568|2495, 0.08020123540993|2986,  0.28319822635379|088,  0.079777157356|279285, 2.133070021475|4486, 0.426201463925188|10 ] ! fortran
write(*, *) '-----------------------------------------------------------------'
write(*, *) 'With NO mask :'
write(*, *) 'sum of abs error: ',   abs(mom%mu - results(1)) +    &  !
                                    abs(mom%va - results(2)) +    &  !
                                    abs(mom%si - results(3)) +    &  !
                                    abs(mom%Sk - results(4)) +    &  !
                                    abs(mom%Ku - results(5)) +    &  !
                                    abs(median - results(6))         !

write(*, *) 'detail: ', mom%mu, mom%va, mom%si, mom%Sk, mom%Ku, median

















deallocate( array, array_mask, array2D, array_mask2D )

!============================================================================== TEST RND NORM

call random_seed()
allocate( vx(1:nn) )

write(*, *) '-----------------------------------------------------------------'
write(*, *) 'Test random normal'
write(*, "(*(a8))")                         "n", "mu", "sigma"
write(*,  "(i8,2f8.4)")                      nn, mu  ,  sigma
write(*,  "(/,'central moments',/,*(i10))") (ipow, ipow = 1, 4)

do iter = 1, 5

   vx = rnorm_vec(nn, mu, sigma)
   vx = vx - sum(vx)/nn

   write(*, "(*(f10.4))") ( sum( vx**ipow )/nn, ipow = 1, 4 )
enddo

write(*, *)            "theoretical"
write(*, "(*(f10.4))") 0.0_R8, sigma**2, 0.0_R8, 3*sigma**4

deallocate( vx )

stop
endprogram test_stat
