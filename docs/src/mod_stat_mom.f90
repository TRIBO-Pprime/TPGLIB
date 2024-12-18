!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: may, 03 2019
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Routines to calculate statistical moments, and some utilities**
!<  </span>

module stat_mom
use data_arch,    only : I4, R8, HIG_R8
use sort_arrays,  only : sort_array2

implicit none

private

   type moment_stat
   !! statistical moments
      real(kind=R8) :: mu  !! *mean*
      real(kind=R8) :: va  !! *variance*
      real(kind=R8) :: si  !! *standard deviation*
      real(kind=R8) :: Sk  !! *skewness*
      real(kind=R8) :: Ku  !! *kurtosis*
      real(kind=R8) :: Ss  !! *fifth moment*
      real(kind=R8) :: Kk  !! *sixth moment*
   endtype moment_stat

public :: moment_stat, calc_moments, calc_median, rnorm_vec, rnorm, scramble, random_normal

contains

   subroutine scramble(tab, lg)
   !================================================================================================
   !! scramble a vector of reals
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=i4), intent(in   )                  :: lg
   real(kind=r8)   , intent(inout), dimension(1:lg) :: tab

      real(kind=r8), dimension(1:lg) :: tmp
      integer(kind=i4) :: i

      call random_number( harvest = tmp(1:lg) )

      call sort_array2(tab_inout = tmp(1:lg),            &  !
                            tab1 = tab(1:lg), n = lg)       !

   return
   endsubroutine scramble


   function random_normal()
   implicit none
   !================================================================================================
   !< @note
   !<
   !< Adapted from the following fortran 77 code
   !<      algorithm 712, collected algorithms from acm.
   !<      this work published in transactions on mathematical software,
   !<      vol. 18, no. 4, december, 1992, pp. 434-435.
   !<
   !<  + The function random_normal() returns a normally distributed pseudo-random
   !<     number with zero mean and unit variance.
   !<  + The algorithm uses the ratio of uniforms method of a.j. kinderman
   !<    and j.f. monahan augmented with quadratic bounding curves.
   !<
   !<     Author:
   !<
   !<     + Alan Miller
   !<     + csiro division of mathematical & information sciences
   !<     + private bag 10, clayton south mdc
   !<     + clayton 3169, victoria, australia
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
      real(kind=r8) :: random_normal
      real(kind=r8) :: s, t, a, b, r1, r2, u, v, x, y, q

      s  =  0.449871_r8
      t  = -0.386595_r8
      a  =  0.19600_r8
      b  =  0.25472_r8
      r1 =  0.27597_r8
      r2 =  0.27846_r8

      !     generate p = (u,v) uniform in rectangle enclosing acceptance region

      do

        call random_number(u)
        call random_number(v)

        v = 1.7156 * (v - 0.5_r8)

      !     evaluate the quadratic form
        x = u - s
        y = abs(v) - t
        q = x**2 + y*(a*y - b*x)

      !     accept p if inside inner ellipse
        if (q < r1) exit
      !     reject p if outside outer ellipse
        if (q > r2) cycle
      !     reject p if outside acceptance region
        if (v**2 < -4.0*log(u)*u**2) exit

      enddo

      !     return ratio of p's coordinates as the normal deviate
      random_normal = v/u

   return
   endfunction random_normal


   subroutine calc_moments(tab, mask, mx, nb_mom)
   !================================================================================================
   !< @note Function to calculate the statistical moments of an array with mask, of shape dim. 1 or 2
   !<
   !< \begin{align*}
   !<     mu &= \frac{1}{n^2}\sum_{i,j=1}^{n}\eta_{i,j} \\
   !<     va &= \frac{1}{n^2}\sum_{i,j=1}^{n}(\eta_{i,j}-\mu)^2 \\
   !<     Sk &= \frac{1}{n^2}\sum_{i,j=1}^{n}\left(\frac{\eta_{i,j}-\mu}{\sigma}\right)^3 \\
   !<     Ku &= \frac{1}{n^2}\sum_{i,j=1}^{n}\left(\frac{\eta_{i,j}-\mu}{\sigma}\right)^4
   !< \end{align*}
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer (kind=I4), intent(in )                          :: nb_mom !! *number of desired moments*
   real    (kind=R8), intent(in ), dimension(..)           :: tab    !! *1D or 2D array*
   logical (kind=I4), intent(in ), dimension(..), optional :: mask   !! *1D or 2D mask*
   type(moment_stat), intent(out)                          :: mx     !! [[moment_stat]] *result*

      integer (kind=I4) :: size_tab

      real(kind=R8),    allocatable, dimension(:) :: tab_tmp
      logical(kind=I4), allocatable, dimension(:) :: msk_tmp

      select rank (tab)

         rank (1)

            if ( present( mask ) ) then

               select rank (mask)

                  rank (1)

                     call calc_moments_1D(tab, mask, mx, nb_mom)

                  rank default

                     stop "bad rank in mask 'calc_moments'"

               endselect

            else

               call calc_moments_1D(tab = tab, mx = mx, nb_mom = nb_mom)

            endif

         rank (2)

            size_tab = product( shape( tab ) )

            allocate( tab_tmp(1:size_tab) )

            tab_tmp = reshape(  tab, [size_tab] )

            if ( present( mask ) ) then

               allocate( msk_tmp(1:size_tab) )

               select rank (mask)

                  rank (2)

                     msk_tmp = reshape( mask, [size_tab] )

                     call calc_moments_1D(tab_tmp, msk_tmp, mx, nb_mom)

                     deallocate( msk_tmp )

                  rank default

                     stop "bad rank in mask 'calc_moments'"

               endselect

            else

               call calc_moments_1D(tab = tab_tmp, mx = mx, nb_mom = nb_mom)

            endif

            deallocate( tab_tmp )

         rank default

            stop "bad rank in 'calc_moments'"

      endselect

   return
   endsubroutine calc_moments


   subroutine calc_moments_1D(tab, mask, mx, nb_mom)
   !================================================================================================
   !! Function to calculate the statistical moments of a 1D array with mask, see [[calc_moments]]
   !------------------------------------------------------------------------------------------------
   implicit none
   integer (kind=I4), intent(in )                         :: nb_mom !! *number of desired moments*
   real    (kind=R8), intent(in ), dimension(:)           :: tab    !! *1D array*
   logical (kind=I4), intent(in ), dimension(:), optional :: mask   !! *1D mask*
   type(moment_stat), intent(out)                         :: mx     !! [[moment_stat]] *result*

      integer(kind=I4) :: lg, nz
      integer(kind=I4) :: i, ii
      real(kind=R8)    :: tmp

      real(kind=R8), allocatable, dimension(:) :: tab_tmp

      lg = size( tab )
      nz = lg

      if ( present(mask) ) then

         nz = count( mask )
         allocate( tab_tmp(1:nz) )

         ii = 0
         do i = 1, lg
            if ( mask(i) ) then
               ii = ii + 1
               tab_tmp(ii) = tab(i)
            endif
         enddo
         if (ii /= nz) stop 'error calc_moments'

      else

         tab_tmp = tab

      endif

      mx%mu = 0
      mx%si = 0
      mx%va = 0
      mx%Sk = 0
      mx%Ku = 0

      mx%Ss = 0
      mx%Kk = 0

      do i = 1, nz
         mx%mu = mx%mu +tab_tmp(i)/nz
      enddo
      if (nb_mom == 1) return

      do i = 1, nz
         mx%va = mx%va +( (tab_tmp(i) -mx%mu)**2 )/nz
      enddo
      mx%si = sqrt( mx%va )
      if (nb_mom==2) return         ! don't go further

      if (mx%si < 1.e-15_R8) then   ! if the standard deviation is too small, quit
         stop 'calc_moments, std too small'
      endif

      do i = 1, nz
         tmp = ( tab_tmp(i) -mx%mu )/mx%si
         mx%Sk = mx%Sk +(tmp**3)/nz
         mx%Ku = mx%Ku +(tmp**4)/nz
      enddo
      if (nb_mom == 4) return

      do i = 1, nz
         tmp = ( tab_tmp(i) -mx%mu )/mx%si
         mx%Ss = mx%Ss +(tmp**5)/nz
         mx%Kk = mx%Kk +(tmp**6)/nz
      enddo

      deallocate( tab_tmp )

   return
   endsubroutine calc_moments_1D


   subroutine calc_median(tab, mask, md)
   !================================================================================================
   !! Function to calculate the median value of a series.
   !!
   !! + Input array containing the values for which the median is to be calculated
   !! + Optional mask to include/exclude certain values from the array
   !! + Output: the calculated median value
   !------------------------------------------------------------------------------------------------
   implicit none
   real   (kind=R8), intent(in ), dimension(:)            :: tab  !! *series 1D array*
   logical(kind=I4), intent(in ), dimension(:), optional  :: mask !! *mask*
   real   (kind=R8), intent(out)                          :: md   !! *result: series median value*

      integer(kind=I4) :: lg, nz ! lg: size of the tab array; nz: number of elements to consider
      integer(kind=I4) :: i, ii  ! i: loop counter; ii: counter for tab_tmp

      real(kind=R8), allocatable, dimension(:) :: tab_tmp ! Temporary array to store filtered values

      md = 0._R8                    ! Initialize the median value to 0

      lg = size( tab )              ! Get the size of the input array
      nz = lg                       ! Initialize the number of elements to lg

      if ( present(mask) ) then     ! Check if a mask is provided

         nz = count( mask )         ! Count the number of true elements in the mask
         allocate( tab_tmp(1:nz) )  ! Allocate memory for the temporary array based on the number of elements to consider

         ii = 0                     ! Initialize the counter for tab_tmp
         do i = 1, lg               ! Loop through each element of the input array
            if ( mask(i) ) then     ! If the element is included in the mask
               ii = ii + 1          ! Increment the counter
               tab_tmp(ii) = tab(i) ! Copy the corresponding value into the temporary array
            endif
         enddo
         if (ii /= nz) stop 'error calc_median' ! Check if the number of copied elements matches nz; if not, stop the program

      else ! If no mask is provided

         tab_tmp = tab ! Copy the input array into tab_tmp

      endif

      if (nz == 1) then             ! If only one element is present
         md = tab_tmp(1)            ! The median is the single element
         return                     ! Exit the subroutine
      endif
      if (nz == 2) then                            ! If two elements are present
         md = 0.5_R8*(tab_tmp(1) + tab_tmp(2))     ! The median is the average of the two elements
         return                                    ! Exit the subroutine
      endif

      call sort_array2(tab_inout = tab_tmp(1:nz), n = nz)         ! Call a subroutine to sort the temporary array

      if ( mod(nz, 2) == 0 ) then                                 ! Check if the number of elements is even
         md = 0.5_R8 * ( tab_tmp( nz/2 ) + tab_tmp( nz/2 + 1) )   ! The median is the average of the two middle elements
      else                                                        ! If the number of elements is odd
         md = tab_tmp( (nz-1)/2 )                                 ! The median is the middle element
      endif

      deallocate( tab_tmp )                                       ! Free the allocated memory for tab_tmp

   return                                                         ! Exit the subroutine
   endsubroutine calc_median


   function rnorm_vec(n, mu, sigma) result(variates)
   !================================================================================================
   !! Vector of reals that follow a normal law
   !!
   !! [source](https://fortran-lang.discourse.group/t/normal-random-number-generator/3724/2)
   !!
   !! authors: Beliavsky, Miller
   !------------------------------------------------------------------------------------------------
   integer(kind=I4), intent(in )                 :: n           !! *vector size*
   real   (kind=R8), intent(in ), optional       :: mu          !! *distribution mean*
   real   (kind=R8), intent(in) , optional       :: sigma       !! *distribution std*
   real   (kind=R8), dimension(1:n)              :: variates    !! *output vector*

      integer(kind=I4) :: i

      do i = 1, n
         variates(i) = rnorm()
      enddo

      if ( present(sigma) ) variates = sigma*variates
      if ( present(mu)    ) variates = variates + mu

   return
   endfunction rnorm_vec


   function rnorm() result(fn_val)
   !================================================================================================
   !! Generate a random normal deviate using the polar method.
   !!
   !! *reference*: marsaglia,g. & bray,t.a. 'a convenient method for generating
   !!              normal variables', siam rev., vol.6, 260-264, 1964.
   !!
   !------------------------------------------------------------------------------------------------
   implicit none
   real(kind=R8) :: fn_val

      real(kind=R8)            :: u, sum
      real(kind=R8),    save   :: v, sln
      logical(kind=I4), save   :: second = .false.
      real(kind=R8), parameter :: one = 1.0_R8, vsmall = tiny( one )

      if (second) then
      ! if second, use the second random number generated on last call

         second = .false.
         fn_val = v*sln

      else
      ! first call; generate a pair of random normals

         second = .true.
         do
            call random_number( u )
            call random_number( v )
            u = scale( u, 1 ) - one
            v = scale( v, 1 ) - one
            sum = u*u + v*v + vsmall         ! vsmall added to prevent log(zero) / zero
            if( sum < one ) exit
         enddo

         sln = sqrt(-scale(log(sum),1)/sum)
         fn_val = u*sln

      endif

   return
   endfunction rnorm

endmodule stat_mom
