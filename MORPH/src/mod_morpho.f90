!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: april, 06 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Morphological operations**
!<  </span>
module morpho
use data_arch, only : I4, R8, PI_R8
use miscellaneous, only : get_unit
use stat_mom, only : calc_median
implicit none

private

public :: count_cell, erode_dilate, def_masque, topology, calcul_normales, surf_area

contains

   subroutine flood(masque, taille, nx, ny, niv)
   !================================================================================================
   !< @note
   !<
   !< Perform some kind of flood fill or connected component labeling on a grid (masque),
   !< starting from an initial '1' element found and spreading out to adjacent '1' elements,
   !< updating them to a specified value or zero if no value (niv) is specified.
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                        :: nx, ny          !! *Input dimensions*
   integer(kind=I4), intent(in   ), optional              :: niv             !! *Optional input level*
   integer(kind=I4), intent(inout), dimension(1:nx, 1:ny) :: masque          !! *Input/output matrix*
   integer(kind=I4), intent(  out)                        :: taille          !! *Output scalar*

   integer(kind=I4) :: i, j, p, q, k, ind, ival
   integer(kind=I4) :: nb_traite

   ! Arrays to hold position indices and a tracking matrix
   integer(kind=I4), dimension(nx*ny)  :: liste_x, liste_y
   integer(kind=I4), dimension(nx, ny) :: deja_fait

      ! Initialize arrays to zero
      liste_x = 0
      liste_y = 0
      deja_fait = 0

      ! Initial position
      i = 0
      j = 0

      ! Find the first occurrence of the element "1" in masque
o:    do i = 1, nx
         do j = 1, ny
            if ( masque(i, j) == 1 ) exit o
         enddo
      enddo o

      ! If no element is found, set output size to zero and exit subroutine
      if ( i == nx + 1 .and. j == ny + 1 ) then
         taille = 0
         return
      endif

      ! Initialize processing of found element
      ind = 1
      liste_x(ind) = i
      liste_y(ind) = j

      ! Mark element as processed
      deja_fait(i, j) = 1

      ! Start processing of elements
      nb_traite = 1

      ! Explore neighbors in 8 directions (up, down, left, right and combinations)
      do
         i = liste_x(nb_traite)
         j = liste_y(nb_traite)

         ! Break the loop if outside bounds
         if ( i*j == 0 ) exit

         ! Check south direction
         k = 1
         do

            p = i
            q = j - k

            if ( q < 1 ) exit
            if (    masque(p, q) == 0 ) exit
            if ( deja_fait(p, q) == 0 ) then
               ind = ind + 1
               liste_x(ind) = p
               liste_y(ind) = q
               deja_fait(p, q) = 1
            endif

            k = k + 1

         enddo

         ! Check south-west direction
         k = 1
         do

            p = i - k
            q = j - k

            if ( p < 1 ) exit
            if ( q < 1 ) exit
            if (    masque(p, q) == 0 ) exit
            if ( deja_fait(p, q) == 0 ) then
               ind = ind + 1
               liste_x(ind) = p
               liste_y(ind) = q
               deja_fait(p, q) = 1
            endif

            k = k + 1

         enddo

         ! Check south-east direction
         k = 1
         do

            p = i + k
            q = j - k

            if ( p > nx ) exit
            if ( q < 1  ) exit
            if (    masque(p, q) == 0 ) exit
            if ( deja_fait(p, q) == 0 ) then
               ind = ind + 1
               liste_x(ind) = p
               liste_y(ind) = q
               deja_fait(p, q) = 1
            endif

            k = k + 1

         enddo

         ! Check north direction
         k = 1
         do

            p = i
            q = j + k

            if ( q > nx  ) exit
            if (    masque(p, q) == 0 ) exit
            if ( deja_fait(p, q) == 0 ) then
               ind = ind + 1
               liste_x(ind) = p
               liste_y(ind) = q
               deja_fait(p, q) = 1
            endif

            k = k + 1

         enddo

         ! Check north-west direction
         k = 1
         do

            p = i - k
            q = j + k

            if ( p < 1  ) exit
            if ( q > ny ) exit
            if (    masque(p, q) == 0 ) exit
            if ( deja_fait(p, q) == 0 ) then
               ind = ind + 1
               liste_x(ind) = p
               liste_y(ind) = q
               deja_fait(p, q) = 1
            endif

            k = k + 1

         enddo

         ! Check north-east direction
         k = 1
         do

            p = i + k
            q = j + k

            if ( p > nx  ) exit
            if ( q > ny ) exit
            if (    masque(p, q) == 0 ) exit
            if ( deja_fait(p, q) == 0 ) then
               ind = ind + 1
               liste_x(ind) = p
               liste_y(ind) = q
               deja_fait(p, q) = 1
            endif

            k = k + 1

         enddo

         ! Check west direction
         k = 1
         do

            p = i - k
            q = j

            if ( p < 1  ) exit
            if (    masque(p, q) == 0 ) exit
            if ( deja_fait(p, q) == 0 ) then
               ind = ind + 1
               liste_x(ind) = p
               liste_y(ind) = q
               deja_fait(p, q) = 1
            endif

            k = k + 1

         enddo

         ! Check east direction
         k = 1
         do

            p = i + k
            q = j

            if ( p > nx  ) exit
            if (    masque(p, q) == 0 ) exit
            if ( deja_fait(p, q) == 0 ) then
               ind = ind + 1
               liste_x(ind) = p
               liste_y(ind) = q
               deja_fait(p, q) = 1
            endif

            k = k + 1

         enddo

         ! Increment processed count
         nb_traite = nb_traite + 1

      enddo

      ! Set output size to number of processed elements
      taille = ind

      ! Update masque values based on presence of niv and processed nodes
      ival = 0
      if ( present( niv ) ) ival = niv
      where( deja_fait == 1) masque = ival ! Fill the connected region in 'masque'

   return
   endsubroutine flood


   subroutine count_cell(msk, long, larg, nbr_cell, med_cell)
   !================================================================================================
   !! Calculate the number of cells in a mask, as well as the cell median size
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                            :: long         !! *2D array length*
   integer(kind=I4), intent(in   )                            :: larg         !! *2D array height*
   integer(kind=I4), intent(out  )                            :: nbr_cell     !! *number of cells*
   real   (kind=R8), intent(out  ), optional                  :: med_cell     !! *median size of the cells*
   integer(kind=I4), intent(inout), dimension(1:long, 1:larg) :: msk          !! *mask*

      integer(kind=I4) :: i, j, k, nb

      real(kind=R8), dimension(1:long*larg) :: cell

      cell(1:long*larg) = 0

      k = 1
      do

         call flood( masque = msk(1:long, 1:larg),       &  ! INOUT
                     taille = nb,                        &  ! OUT
                     nx     = long,                      &  ! IN
                     ny     = larg,                      &  ! IN
                     niv    = k + 1)                        ! IN

         if ( nb == 0 ) exit

         cell(k) = nb

         k = k + 1

      enddo

      nbr_cell = k - 1
      med_cell = 0._R8

      if (.not.present(med_cell)) return

      if (nbr_cell > 0) then

         call calc_median(  tab = cell(1:nbr_cell),   &  ! IN
                             md = med_cell )             ! OUT

         med_cell = 100*med_cell/(long*larg)

      else

         med_cell = 0

      endif

   return
   endsubroutine count_cell


   subroutine erode_dilate(msk, long, larg, val, act)
   !================================================================================================
   !! Perform erosion or dilation on a binary mask depending on the value of act.
   !! The operations utilize a defined kernel to affect neighboring pixels based on the specified val
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                             :: long  !! *2D mask length*
   integer(kind=I4), intent(in   )                             :: larg  !! *2D mask height*
   integer(kind=I4), intent(in   )                             :: val   !! *size of the structuring element for the erosion/dilation operation*
   character(len=5), intent(in   )                             :: act   !! *action to be performed, either "erode" or another operation, presumably "dilate".*
   integer(kind=I4), intent(inout), dimension(1:long,1:larg)   :: msk   !! *2D mask*

      integer(kind=I4) :: i, j, fin, ival
      integer(kind=I4), dimension(-val:val, -val:val) :: ker

      ! If the action specified is "erode", the mask is inverted (0s become 1s and 1s become 0s). This prepares the mask for the erosion operation.
      if (act =="erode") msk(1:long, 1:larg) = 1 - msk(1:long, 1:larg)
      fin = +2

      ! The forall construct fills the kernel with the value fin (2) for all points (i, j) within a circular area defined by the radius val.
      ! The center of the kernel (ker(0, 0)) is explicitly set to 0.
      ker(-val:val, -val:val) = 0
      forall (i = -val:+val, j = -val:+val, i**2 +j**2 <= val**2) ker(i, j) = fin
      ker(0, 0) = 0

      do j = 1 +val, larg - val
      do i = 1 +val, long - val

         ! If the pixel value at position (i, j) in the mask is even, the loop skips any further processing for that pixel and continues to the next iteration.
         if (mod(msk(i, j), 2) == 0) cycle
         ! If the pixel value is odd, the corresponding area in the mask defined by the kernel is updated by adding the kernel values to it.
         ! This implements the dilation operation.
         msk(i - val:i + val, j - val:j + val) = msk(i - val:i + val, j - val:j + val) + ker(-val:val, -val:val)

      enddo
      enddo

      ! This line sets all elements of the mask that are greater than or equal to 1 to 1, effectively binarizing the mask after dilation.
      where(msk(1:long,1:larg)>=1) msk = 1

      ! If the action specified is "erode", the mask is inverted again, restoring it to its original form after dilation.
      if (act=="erode") msk(1:long,1:larg) = 1 -msk(1:long,1:larg)

   return
   endsubroutine erode_dilate


   subroutine topology(tab, long, larg, res)
   !================================================================================================
   !< @note
   !<
   !< The function performs the following operations on a surface:
   !<
   !< + mask the heights heights that are above 85 % of the (100% - 15 %) core, to avoid pits and peaks
   !< + erode then dilate (opening) the mask
   !< + count cells and the median size
   !<
   !< Reproduce the preceding steps with thresholds 15% and 95%.
   !<
   !< The results are put in the vector 'res'
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                             :: long  !! *2D array length*
   integer(kind=I4), intent(in )                             :: larg  !! *2D array height*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg)  :: tab   !! *heights 2D array*
   real   (kind=R8), intent(out), dimension(1:6)             :: res   !! *results*

      real   (kind=R8), dimension(1:long, 1:larg) :: tab_tmp1
      integer(kind=I4), dimension(1:long, 1:larg) :: msk

      real   (kind=R8) :: mintab, maxtab, top01, top02, med_cell01, med_cell02
      integer(kind=I4) :: nbr_cell01, nbr_cell02

      mintab = minval( tab(1:long, 1:larg) )
      maxtab = maxval( tab(1:long, 1:larg) )

      tab_tmp1(1:long, 1:larg) = ( tab(1:long, 1:larg) - mintab )/( maxtab - mintab )

      call def_masque(msk = msk,          &  ! OUT
                      tab = tab_tmp1,     &  ! IN
                     long = long,         &  ! IN
                     larg = larg,         &  ! IN
                    crit1 = 0.15_R8,      &  ! IN
                    crit2 = 0.85_R8,      &  ! IN
                      top = top01)           ! OUT
      !......................................!
      call erode_dilate(msk = msk,        &  ! INOUT
                       long = long,       &  ! IN
                       larg = larg,       &  ! IN
                        val = 5,          &  ! IN
                        act = "erode")       ! IN

      call erode_dilate(msk = msk,        &  ! INOUT
                       long = long,       &  ! IN
                       larg = larg,       &  ! IN
                        val = 5,          &  ! IN
                        act = "dilat")       ! IN

      call count_cell(msk = msk,          &  ! INOUT
                     long = long,         &  ! IN
                     larg = larg,         &  ! IN
                 nbr_cell = nbr_cell01,   &  ! IN
                 med_cell = med_cell01)      ! IN
      !......................................!
      call def_masque(msk = msk,          &  ! OUT
                      tab = tab_tmp1,     &  ! IN
                     long = long,         &  ! IN
                     larg = larg,         &  ! IN
                    crit1 = 0.15_R8,      &  ! IN
                    crit2 = 0.95_R8,      &  ! IN
                      top = top02)           ! OUT

      call erode_dilate(msk = msk,        &  ! INOUT
                       long = long,       &  ! IN
                       larg = larg,       &  ! IN
                        val = 5,          &  ! IN
                        act = "erode")       ! IN

      call erode_dilate(msk = msk,        &  ! INOUT
                       long = long,       &  ! IN
                       larg = larg,       &  ! IN
                        val = 5,          &  ! IN
                        act = "dilat")       ! IN

      call count_cell(msk = msk,          &  ! INOUT
                     long = long,         &  ! IN
                     larg = larg,         &  ! IN
                 nbr_cell = nbr_cell02,   &  ! OUT
                 med_cell = med_cell02)      ! OUT

      res(1:6) = [real(nbr_cell01, kind=R8), med_cell01, top01, &  !
                  real(nbr_cell02, kind=R8), med_cell02, top02]

   return
   endsubroutine topology


   subroutine def_masque(msk, tab, long, larg, crit1, crit2, top)
   !================================================================================================
   !! Height mask without deepest pits and highest peaks
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long     !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg     !! *2D array height*
   real   (kind=R8), intent(in )                            :: crit1    !! *%age for deepest pits to remove*
   real   (kind=R8), intent(in )                            :: crit2    !! *%age for highest peaks to remove*
   real   (kind=R8), intent(out)                            :: top      !! *%age of surface masked*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab      !! *heights 2D array*
   integer(kind=I4), intent(out), dimension(1:long, 1:larg) :: msk      !! *mask*

      msk(1:long, 1:larg) = 0
      ! height masked : height above c1 + c2 * (1 - c1)
      ! ex: c1 = 15 %, c2 = 85 % -> heights that are above 85 % of the (1 - 15 %) core are masked
      where( tab > crit1 + crit2 - crit1 * crit2 ) msk = 1
      top = 100 * real( sum( msk(1:long, 1:larg) ), kind=R8 ) / ( long * larg )

   return
   endsubroutine def_masque


   subroutine make_mask(x0, y0, a, b, shap, msk, long, larg)
   !================================================================================================
   !! Mask a region within a given shape (for the moment an ellipsis)
   !------------------------------------------------------------------------------------------------
   implicit none
   integer  (kind=I4), intent(in   )                              :: long  !! *2D array length*
   integer  (kind=I4), intent(in   )                              :: larg  !! *2D array height*
   integer  (kind=I4), intent(in   )                              :: x0    !! *mask shape center 1st coordinate*
   integer  (kind=I4), intent(in   )                              :: y0    !! *mask shape center 2nd coordinate*
   integer  (kind=I4), intent(in   )                              :: a     !! *ellipsis semi-length*
   integer  (kind=I4), intent(in   )                              :: b     !! *ellipsis semi-height*
   character(len = 8), intent(in   )                              :: shap  !! *kind of mask shape*
   integer  (kind=I4), intent(inout), dimension(1:long, 1:larg)   :: msk   !! *mask*

      integer(kind=I4) :: i, j

      select case (shap)

         case('ellipsis')
            forall( i = 1:long, j = 1:larg, (real(i - x0, kind = R8) / a)**2 + (real(j - y0, kind = R8) / b)**2 < 1. ) msk(i, j) = 1

         case default
            stop 'make_mask, bad choice'

      endselect

   return
   endsubroutine make_mask


   subroutine make_composite_mask(msk, n_cells, locus, a, b, width, height, shap, long, larg)
   !================================================================================================
   !< @note
   !<
   !< The subroutine generates a composite mask based on specified parameters.
   !< The mask is filled with a number of shapes (in this case, ellipses) placed randomly within a defined area, while ensuring that the shapes do not overlap.
   !<
   !< + If the `locus` is set to 'center', the subroutine immediately centers a single shape in the mask.
   !< + If not centered, it attempts to randomly place the specified number of shapes (`n_cells`) within the mask.
   !<   The algorithm limits attempts to place shapes to avoid overlaps.
   !< + The routine uses random numbers to determine the position of each shape, ensuring that they fit within the given dimensions and do not exceed the boundaries of the mask.
   !< + A separate temporary mask (`msk_quad`) is used to track where shapes have already been placed, preventing overlaps.
   !< + The process continues until the desired number of shapes is successfully placed or until a maximum number of attempts is reached (to avoid infinite loops).
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer  (kind=I4), intent(in)                              :: n_cells   !! *Number of cells to create in the composite mask*
   integer  (kind=I4), intent(in)                              :: long      !! *Length (dimensions) of the mask*
   integer  (kind=I4), intent(in)                              :: larg      !! *Width (dimensions) of the mask*
   integer  (kind=I4), intent(in)                              :: width     !! *Width of each shape to be drawn in the mask*
   integer  (kind=I4), intent(in)                              :: height    !! *Height of each shape to be drawn in the mask*
   integer  (kind=I4), intent(in)                              :: a         !! *ellipsis first parameter*
   integer  (kind=I4), intent(in)                              :: b         !! *ellipsis second parameter*
   character(len = 8), intent(in)                              :: shap      !! *Type of shape to be drawn (here, 'ellipsis' is used)*
   character(len = 6), intent(in)                              :: locus     !! *Position of the shape (e.g., 'center' for centering)*
   integer  (kind=I4), intent(out), dimension(1:long, 1:larg)  :: msk       !! *Output mask that will be filled with shapes*

      ! Declaration of local variables
      integer(kind=I4) :: icells, x0, y0, iter, itry
      real   (kind=R8) :: x, y

      ! Temporary masks to manage the shapes
      integer(kind=I4), dimension(1:long, 1:larg) :: msk_quad, msk_shap

      ! Initialize the output mask to zero
      msk(1:long, 1:larg) = 0

      ! If the position is centered, calculate the central coordinates
      if (locus == 'center') then

         x0 = long / 2 + 1       ! x-coordinate of the center
         y0 = larg / 2 + 1       ! y-coordinate of the center

         ! Call the subroutine to create a centered ellipse mask
         call make_mask(x0 = x0,          &  !
                        y0 = y0,          &  !
                         a = a,           &  !
                         b = b,           &  !
                      shap = 'ellipsis',  &  !
                       msk = msk,         &  !
                      long = long,        &  !
                      larg = larg)           !
         return

      endif

      ! Initialize the attempt counter for placing cells
      itry = 0
l1:   do

      itry = itry + 1               ! Increment the attempt counter
      if (itry > 100) exit l1       ! Limit the number of tries to 100

      msk_quad(1:long, 1:larg) = 0  ! Reset the quadrilateral mask
      msk(1:long, 1:larg) = 0       ! Reset the final mask

      icells = 0                    ! Cell placement counter
      iter = 0                      ! Iteration counter for placement

l2:   do

         iter = iter + 1            ! Increment the iteration counter
         if (iter > 200) cycle l1   ! Limit the number of iterations to 200

         ! Display the current state of the iteration and the number of cells
         write(*,*) iter, icells

         call random_number(x)                           ! Generate a random number for x
         x0 = width / 2 + 1 + x * (long - width - 1)     ! Calculate the x-coordinate of the center of the shape
         call random_number(y)                           ! Generate a random number for y
         y0 = height / 2 + 1 + y * (larg - height - 1)   ! Calculate the y-coordinate of the center of the shape

         msk_shap(1:long, 1:larg) = 0                    ! Reset the temporary shape mask

         ! Define the shape mask over a square area centered on (x0, y0)
         msk_shap(x0 - width  / 2:x0 + width / 2,        &  !
                  y0 - height / 2:y0 + height / 2) = 1      !

         ! Check if the shape does not overlap with an already placed shape
         if (sum(msk_quad(1:long, 1:larg) * msk_shap(1:long, 1:larg)) > 0) cycle l2

         ! Call the subroutine to create an ellipse mask at coordinates (x0, y0)
         call make_mask(x0 = x0,             &  !
                        y0 = y0,             &  !
                         a = a,              &  !
                         b = b,              &  !
                      shap = 'ellipsis',     &  !
                       msk = msk,            &  !
                      long = long,           &  !
                      larg = larg)              !

         msk_quad(x0 - width  / 2:x0 + width / 2,        &  !
                  y0 - height / 2:y0 + height / 2) = 1      ! Mark the area of the placed shape

         icells = icells + 1                                ! Increment the count of placed cells

         if (icells == n_cells) exit l1                     ! Exit if the desired number of cells has been reached

      enddo l2

      enddo l1

   return
   endsubroutine make_composite_mask


   subroutine calcul_normales(tab_in, long, larg, scale_xyz, cone_angle, hori, print_mask)
   !================================================================================================
   !! Function that returns the fraction of surface nearly horizontal (normal less than 5 degrees
   !!        from a vertical line)
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                           :: long         !! *surface array length*
   integer(kind=I4), intent(in )                           :: larg         !! *surface array width*
   real   (kind=R8), intent(in )                           :: cone_angle   !! *cone angle*
   real   (kind=R8), intent(out)                           :: hori         !! *fraction of facets nearly horizontal*
   logical(kind=I4), intent(in ), optional                 :: print_mask   !! *mask output ?*
   real   (kind=R8), intent(in ), dimension(1:long,1:larg) :: tab_in       !! *surface array*
   real   (kind=R8), intent(in ), dimension(1:3)           :: scale_xyz    !! *lag along x, y and scale z*

      integer(kind=I4) :: i, j, ino, nb_no, ua
      real   (kind=R8) :: z1, z2, z3, z4
      real   (kind=R8) :: n1x, n1y, n1z, n2x, n2y, n2z, nx, ny, nz, norme
      real   (kind=R8) :: phi, angle_moy, ech_x, ech_y, ech_z, ech_r, d13, d24, d42

      character(len=16) :: str

      real(kind=R8), allocatable, dimension(:)     :: angle
      real(kind=R8), allocatable, dimension(:,:,:) :: vec

      nb_no = long*larg !  nombre de noeuds auquels on associera une normale

      allocate( angle(1:nb_no) )                ! recevra les angles à la verticale pour chaque noeud
      allocate(   vec(1:long, 1:larg, 1:3) )    ! recevra pour chaque noeud la contribution des facettes triangulaires connectées

      ech_x = scale_xyz(1) !SCALE_IMG%dx * unit2IUf(SCALE_IMG%dx_unit)
      ech_y = scale_xyz(2) !SCALE_IMG%dy * unit2IUf(SCALE_IMG%dy_unit)
      ech_z = scale_xyz(3)

      ech_r = ech_x/ech_y

      d13 = (2./3.)*sqrt(  ech_x**2    +  ech_y**2    )/2
      d42 = (2./3.)*sqrt( (ech_x/2)**2 +  ech_y**2    )
      d24 = (2./3.)*sqrt(  ech_x**2    + (ech_y/2)**2 )

      ! Raisonnement sur chaque carré du domaine
      vec(1:long, 1:larg, 1:3) = 0
      do j = 1, larg -1
      do i = 1, long -1

         z1 = tab_in(i   , j   ) * ech_z
         z4 = tab_in(i   , j +1) * ech_z
         z2 = tab_in(i +1, j   ) * ech_z
         z3 = tab_in(i +1, j +1) * ech_z

         ! Triangle 1
         n1x = (z1 -z2)/ech_x
         n1y = (z1 -z4)/ech_y
         n1z = 1

         norme = sqrt( n1x**2 +n1y**2 +n1z**2)
         n1x = n1x/norme
         n1y = n1y/norme
         n1z = n1z/norme


         ! Triangle 2
         n2x = (z4 -z3)/ech_x
         n2y = (z2 -z3)/ech_y
         n2z = 1

         norme = sqrt( n2x**2 +n2y**2 +n2z**2)
         n2x = n2x/norme
         n2y = n2y/norme
         n2z = n2z/norme

         ! triangle 1
         vec(i +0, j +0, 1:3) = vec(i +0, j +0, 1:3) +[n1x, n1y, n1z]/d13 ! node 1
         vec(i +0, j +1, 1:3) = vec(i +0, j +1, 1:3) +[n1x, n1y, n1z]/d24 ! node 2
         vec(i +1, j +0, 1:3) = vec(i +1, j +0, 1:3) +[n1x, n1y, n1z]/d42 ! node 4

         ! triangle 2
         vec(i +1, j +1, 1:3) = vec(i +1, j +1, 1:3) +[n2x, n2y, n2z]/d13 ! node 3
         vec(i +0, j +1, 1:3) = vec(i +0, j +1, 1:3) +[n2x, n2y, n2z]/d42 ! node 2
         vec(i +1, j +0, 1:3) = vec(i +1, j +0, 1:3) +[n2x, n2y, n2z]/d24 ! node 4


      enddo
      enddo

      ino = 0
      do j = 1, larg
      do i = 1, long

         nx = vec(i, j, 1)
         ny = vec(i, j, 2)
         nz = vec(i, j, 3)

         norme = sqrt( nx**2 +ny**2 +nz**2)

         nx = nx/norme
         ny = ny/norme
         nz = nz/norme

         phi = acos( nz )*180./PI_R8
         ino = ino +1
         angle(ino) = phi

      enddo
      enddo

      ! angle moyen de la surface
      angle_moy = 0!sum( angle(1:nb_no) )/nb_no

      ! tout ce qui dépasse l'angle moyen (à cone_angle° près) est masqué à 0
      where( abs(angle - angle_moy) < cone_angle )
         angle = 1.
      elsewhere
         angle = 0.
      endwhere
      ! pourcentage de facettes quasi horizontales
      hori = 100*sum( angle(1:nb_no) )/nb_no

      if ( present(print_mask) ) then

         str = repeat( ' ', len(str) )
         write( str, '(I12.12)' ) ino
         call get_unit(ua)
         open( unit = ua, file = "out/mask_angle.txt")
            write( ua, '('//trim(str)//'I2.1)' ) ( int( angle(i) ), i = 1, long * larg )
         close( ua )

      endif

      deallocate( angle, vec )

   return
   endsubroutine calcul_normales


   subroutine surf_area(tab_in, long, larg, scale_xyz, aire)
   !================================================================================================
   !! Function that returns the relative area of a surface minus 1.
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long        !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg        !! *2D array width*
   real   (kind=R8), intent(in ), dimension(1:3)            :: scale_xyz   !! *scale along x, y, z*
   real   (kind=R8), intent(out)                            :: aire        !! *computed area*
   real   (kind=R8), intent(in ), dimension(1:long,1:larg)  :: tab_in      !! *surface array*

      integer(kind=I4) :: i, j
      real   (kind=R8) :: z1, z2, z3, z4, hx, hy, hz

      hx = scale_xyz(1)
      hy = scale_xyz(2)
      hz = scale_xyz(3)

      ! Raisonnement sur chaque carré du domaine
      aire = 0.
      do j = 1, larg -1
      do i = 1, long -1

         z1 = tab_in(i   , j   )*hz
         z2 = tab_in(i   , j +1)*hz
         z3 = tab_in(i +1, j +1)*hz
         z4 = tab_in(i +1, j   )*hz

         aire = aire +0.5_R8*( sqrt( 1._R8 + ( (z1-z2)/hx )**2 + ( (z1-z4)/hy )**2 ) +    &  !
                               sqrt( 1._R8 + ( (z3-z2)/hy )**2 + ( (z3-z4)/hx )**2 ) )       !

      enddo
      enddo
      aire = aire/( (long -1)*(larg -1) ) - 1._R8

   return
   endsubroutine surf_area

endmodule morpho
