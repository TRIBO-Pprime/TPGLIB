!< author: Arthur Francisco
!<  version: 1.1.0
!<  date: april, 07 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Returns the Firestone Abbott's curve as well as some ISO 25178 parameters**
!<  </span>
module abbott
use data_arch,       only : I4, R4, R8, UN, EPS_R8, PI_R8
use miscellaneous,   only : get_unit
use stat_mom,        only : calc_moments, moment_stat
use sort_arrays,     only : sort_array2
use least_squares,   only : moindres_carres_lineaire
use gnufor,          only : run_gnuplot
use pikaia_oop,      only : pikaia_class

!$ use omp_lib
implicit none

contains

   subroutine abbott_param(tab, lg, nom, curves, results, omp)
   !================================================================================================
   !< @note Function that returns the Abbott's curve in a svg file as well as smrk1, smrk2, spk, svk, sk
   !
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                              :: lg       !! *surface total number of points*
   character(len=*), intent(in   )                              :: nom      !! *output generic name*
   logical(kind=I4), intent(in   ), dimension(1:3)              :: curves   !! *if true, generates a svg drawing*
                                                                            !! 1: histogram 2: Abbott 3: tangent fit
   real   (kind=R8), intent(inout), dimension(1:lg )            :: tab      !! *surface in a 1D vector*
   real   (kind=R8), intent(  out), dimension(1:11 )            :: results  !! *surface parameters output*
                                                                            !!
                                                                            !!  +  1 **smrk1**, iso 25178
                                                                            !!  +  2 **smrk2**, iso 25178
                                                                            !!  +  3 **spk**  , iso 25178
                                                                            !!  +  4 **svk**  , iso 25178
                                                                            !!  +  5 **off1** , ordonnée de spk
                                                                            !!  +  6 **off2** , ordonnée de svk
                                                                            !!  +  7 **sk**   , iso 25178
                                                                            !!  +  8 **core slope**
                                                                            !!  +  9 **adjustment factor** (tangent fit)
                                                                            !!  + 10 **coeffa_tan**        (tangent fit)
                                                                            !!  + 11 **coeffb_tan**        (tangent fit)
   logical(kind=I4), intent(in   )                              :: omp      !! *if true, openmp used*

      integer(kind=I4), parameter :: nb_points = 4 * 4096   ! nb points kept for the curve

      integer(kind=I4) :: i, ua, icat, ncat, nb_paquets, inc, reste
      integer(kind=I4) :: len_reg, beg_reg, end_reg
      integer(kind=I4) :: smrk1, smrk2, ios
      integer(kind=I4) :: status                            ! PIKAIA: status

      real(kind=R8) :: hmin, hmax, delt, seuil, reduction
      real(kind=R8) :: tmp, off1, off2, spk, svk, sk
      real(kind=R8) :: f                                    ! PIKAIA: best cost

      type(moment_stat) :: mx

      real(kind=R8), dimension(1:2) :: vec_reg

      real(kind=R8), dimension(1:4) :: xx                   ! PIKAIA: chromosom
      real(kind=R8), dimension(1:4) :: xl, xu               ! PIKAIA: lower and upper bonds of xx

      integer(kind=I4), allocatable, dimension(:) :: tab_abbott

      real(kind=R8), allocatable, dimension(:)   :: tab_moy
      real(kind=R8), allocatable, dimension(:,:) :: jac_reg

      type(pikaia_class) :: p                               ! PIKAIA: class instanciation


      ! center heights
      call calc_moments(   tab = tab(1:lg),  &  !
                            mx = mx,         &  !
                        nb_mom = 2 )            !

      tab(1:lg) = tab(1:lg) - mx%mu

      ! sort in decreasing order
      tab(1:lg) = -tab(1:lg)
      call sort_array2( tab_inout = tab(1:lg),    &  !
                        n         = lg )             !
      tab(1:lg) = -tab(1:lg)

      ! reduction factor to stay between 1 and 100 because there are
      ! too many points in array tab
      reduction = nb_points / 100.

      ! array of points for the curve
      allocate( tab_moy(1:nb_points) )

      ! nb heights of tab representing each point of the curve
      nb_paquets = lg/( nb_points - 1 )
      ! ... and remaining heights
      reste      = lg - nb_paquets*( nb_points - 1 )

      ! array of points; the first point and the last one are
      ! the first and last height to keep the extrema.
      inc = 1
      do i = 1, nb_points - 1

         tab_moy(i) = tab(inc)
         inc = inc + nb_paquets
         if ( i == nb_points/2 ) inc = inc + reste - 1 ! the remaing heights are added to
                                                       ! the middle points, where it makes less difference
      enddo
      tab_moy(nb_points) = tab(lg)

      if ( curves(1) ) then ! print histogram ?

         ncat = 100                                      !  nombre de "classes" pour l'histogramme
         allocate( tab_abbott(1:ncat) ) ; tab_abbott(1:ncat) = 0

         hmin = tab(lg) - 100*EPS_R8                     !  pour être sûr d'attraper le min
         hmax = tab( 1) + 100*EPS_R8                     !  idem pour le max

         delt = ( hmax - hmin ) / ncat                   ! categories width (heights)

         seuil = hmax - delt

         i = 0
         do icat = 1, ncat

            do
               i = i + 1 ; if ( i > lg ) exit
               if ( tab(i) < seuil ) exit
               tab_abbott(icat) = tab_abbott(icat) + 1
            enddo

            seuil = seuil - delt

         enddo

         call get_unit( ua )
         open( unit = ua, file = trim(nom)//'_histo.txt')

            do icat = 1, ncat

               write(ua, *) real( hmax + delt/2 - icat * delt, kind = R4 ), 100 * real( tab_abbott(icat), kind = R4 ) / lg

            enddo

         close(ua)

         deallocate( tab_abbott )

      endif

      ! nb points for the line regression! between 30% and 60% of the curve
      len_reg = int( 0.4 * nb_points )
      beg_reg = int( 0.3 * nb_points )
      end_reg = beg_reg + len_reg - 1

      allocate( jac_reg(1:len_reg, 1:2) )

      ! jacobian regression : vec_reg(1) * Xi + vec_reg(2) * 1
      jac_reg( 1:len_reg, 1 ) = [ ( beg_reg + i, i = 0, len_reg - 1 ) ]
      jac_reg( 1:len_reg, 2 ) = 1._R8

      call moindres_carres_lineaire(   nb_var = 2,                            &  ! number of parameters to be determined
                                       nb_pts = len_reg,                      &  ! number of points for function evaluation
                                       hij    = tab_moy( beg_reg:end_reg ),   &  ! vector of evaluation points (1:nb_pts)
                                       beta   = vec_reg( 1:2 ),               &  ! parameters vector (1:nb_var)
                                       Jf     = jac_reg( 1:len_reg, 1:2 ) )      ! Jacobian (1:nb_pts, 1:nb_var)

      ! f(smrk1) coordinates
      off1 = vec_reg(2)                               ! (0._R8, off1)
      do i = 1, nb_points
         if ( tab_moy(i) < off1 ) exit
      enddo
      smrk1 = i                                       ! (UN*smrk1, off1)

      tmp  = 0.
      do i = 1, smrk1
         tmp = tmp +( tab_moy(i) - off1 )
      enddo
      spk = 2*tmp/smrk1                               ! areas equivalency
      spk = abs(spk)

      ! f(smrk2) coordinates
      off2  = vec_reg(1) * nb_points + vec_reg(2)     ! (nb_points, off2)
      do i = nb_points, 1, -1
         if ( tab_moy(i) > off2 ) exit
      enddo
      smrk2 = i                                       ! (UN*smrk2, off2)

      tmp = 0.
      do i = smrk2, nb_points
         tmp = tmp +( tab_moy(i) - off2 )
      enddo
      svk = 2*tmp/( nb_points - smrk2 + 1 )           ! areas equivalency

      sk = tab_moy(smrk1) - tab_moy(smrk2)

      if ( curves(2) ) then ! print Abbott ?
         call get_unit( ua )

         open(unit = ua, file = trim(nom)//'.dat')

            do i = 1, nb_points
               write(ua, *) 100 * real( i - 1, kind = R4 )/( nb_points - 1),        &  !
                                  real( tab_moy(i), kind = R4),                     &  !
                                  real( vec_reg(1) * i + vec_reg(2), kind = R4 )       !
            enddo

         close(ua)

         open(unit = ua, file = trim(nom)//'.gpl', status = 'replace', iostat = ios )

            write( ua, '(a)' ) 'set terminal svg dashed size 350,262 font "Verdana, 10"'
            write( ua, '(a)' ) 'set output "'//trim(nom)//'.svg"'
            write( ua, '(a)' ) 'set title "Abbott"'
            write( ua, '(a)' ) 'set xlabel "%age"'
            write( ua, '(a)' ) 'set ylabel "h"'
            write( ua, '(a)' ) "set style line 1 lc rgb 'dark-green' lt 1 lw 1"
            write( ua, '(a)' ) "set style line 2 lc rgb 'dark-green' lt 5 lw 1"
            ! tracé du segment horizontal de A1 (cf norme 25178)
            write( ua, * ) 'set arrow 1 from ',        0._R8/reduction, ',', off1, ' to ',     UN*smrk1/reduction, ',',      off1, ' nohead front ls 1'
            ! tracé du segment incliné de A1
            write( ua, * ) 'set arrow 2 from ',     UN*smrk1/reduction, ',', off1, ' to ',        0._R8/reduction, ',', spk +off1, ' nohead front ls 1'
            ! tracé de la droite haute horizontale
            write( ua, * ) 'set arrow 5 from ',        0._R8/reduction, ',', off1, ' to ', UN*nb_points/reduction, ',',      off1, ' nohead front ls 2'
            ! tracé du segment horizontal de A2 (cf norme 25178)
            write( ua, * ) 'set arrow 3 from ',     UN*smrk2/reduction, ',', off2, ' to ', UN*nb_points/reduction, ',',      off2, ' nohead front ls 1'
            ! tracé du segment incliné de A2
            write( ua, * ) 'set arrow 4 from ',     UN*smrk2/reduction, ',', off2, ' to ', UN*nb_points/reduction, ',', svk +off2, ' nohead front ls 1'
            ! tracé de la droite basse horizontale
            write( ua, * ) 'set arrow 6 from ', UN*nb_points/reduction, ',', off2, ' to ',        0._R8/reduction, ',',      off2, ' nohead front ls 2'
            ! tracé du segment vertical correspondant à sk
            write( ua, * ) 'set arrow 7 from ', 50., ',', off2, ' to ', 50., ',', off2 +sk, ' nohead front ls 1'
            write( ua, '(a)' ) 'plot "' // trim(nom)//'.dat" ' // 'using 1:2  title "Abbott-Firestone curve" with lines,\'
            write( ua, '(a)' )     ' "' // trim(nom)//'.dat" ' // 'using 1:3  notitle                        with lines ls 1'

         close(unit = ua)

         call run_gnuplot (trim(nom)//'.gpl')
      endif

      deallocate( jac_reg )

      !-----------------------------------------------------------

      call sort_array2( tab_inout = tab_moy(1:nb_points), &  !
                        n         = nb_points )              !

      xx(1:4) = 0.0_R8 ! vector of parameters
      xl(1:4) = 0.0_R8 ! lower bound
      xu(1:4) = 1.0_R8 ! upper bound

      !initialize the class:
      call p%init(           n = 4,                   &  ! IN           ; the parameter space dimension, i.e., the number of
                                                         !                adjustable parameters (size of the x vector).
                            xl = xl,                  &  ! IN, DIM(n)   ; vector of lower bounds for x
                            xu = xu,                  &  ! IN, DIM(n)   ; vector of upper bounds for x
                             f = cost,                &  !              ; user-supplied scalar function of n variables, which
                                                         !                must have the pikaia_func procedure interface.
                        status = status,              &  ! OUT          ; status output flag (0 if there were no errors)
                            np = 100,                 &  ! IN, OPT      ; number of individuals in a population (default is 100)
                          ngen = 1000,                &  ! IN, OPT      ; maximum number of iterations
                            nd = 9,                   &  ! IN           ; number of significant digits (i.e., number of genes)
                                                         !                retained in chromosomal encoding
                        pcross = 0.85_R8,             &  ! IN, OPT      ; crossover probability; must be <= 1.0 (default is 0.85)
                                                         !                If crossover takes place, either one or two splicing points are used,
                                                         !                with equal probabilities
                        pmutmn = 0.0005_R8,           &  ! IN, OPT      ; minimum mutation rate; must be >= 0.0 (default is 0.0005)
                        pmutmx = 0.25_R8,             &  ! IN, OPT      ; maximum mutation rate; must be <= 1.0 (default is 0.25)
                          pmut = 0.005_R8,            &  ! IN, OPT      ; initial mutation rate; should be small (default is 0.005)
                                                         !                (Note: the mutation rate is the probability that any one gene
                                                         !                 locus will mutate in any one generation.)
                          imut = 2,                   &  ! IN, OPT      ; mutation mode; 1/2/3/4/5 (default is 2).
                                                         !                1=one-point mutation, fixed rate.
                                                         !                2=one-point, adjustable rate based on fitness.
                                                         !                3=one-point, adjustable rate based on distance.
                                                         !                4=one-point+creep, fixed rate.
                                                         !                5=one-point+creep, adjustable rate based on fitness.
                                                         !                6=one-point+creep, adjustable rate based on distance.
                          fdif = 1._R8,               &  ! IN, OPT      ; relative fitness differential; range from 0 (none) to 1 (maximum).
                                                         !                (default is 1.0)
                          irep = 3,                   &  ! IN, OPT      ; reproduction plan; 1/2/3=Full generational replacement/
                                                         !                                         Steady-state-replace-random/
                                                         !                                         Steady- state-replace-worst (default is 3)
                        ielite = 0,                   &  ! IN, OPT      ; elitism flag; 0/1=off/on (default is 0)
                                                         !                (Applies only to reproduction plans 1 and 2)
                          ivrb = 0,                   &  ! IN, OPT      ; printed output 0/1/2=None/Minimal/Verbose
               convergence_tol = 1.0e-6_R8,           &  ! IN, OPT      ; convergence tolerance; must be > 0.0 (default is 0.0001)
            convergence_window = 400,                 &  ! IN, OPT      ; convergence window; must be >= 0 This is the number of consecutive
                                                         !                solutions within the tolerance for convergence to be declared (default is 20)
            initial_guess_frac = 0.1_R8,              &  ! IN, OPT      ; raction of the initial population to set equal to the initial guess.
                                                         !                Range from 0 (none) to 1.0 (all). (default is 0.1 or 10%).
                         iseed = 999)                    ! IN, OPT      ; random seed value; must be > 0 (default is 999)

      !Now call pikaia:
      call p%solve(      x = xx(1:4),     &  ! INOUT, DIM(*) ;
                         f = f,           &  !   OUT         ;
                    status = status,      &  !   OUT         ;
                       omp = omp )           ! IN

      xx(3:4) = xx(3:4) * mx%si

      if ( curves(3) ) then ! print tangent fit ?

         call get_unit(ua)

         open( unit = ua, file = trim(nom)//'_tan.datt')

            do i = 1, nb_points

               write(ua, *) 100 * real( i - 1, kind = R4) / (nb_points - 1),                                      &  !
                                  real( tab_moy(i), kind = R4 ),                                                  &  !
                                  real( tg( real( i - 1, kind = R8) / ( nb_points - 1 ), xx(1:4) ), kind = R4 )      !

            enddo

         close(ua)

         open( unit = ua, file = trim(nom)//'_tan.gplt', status = 'replace', iostat = ios )

            write( ua, '(a)' ) 'set terminal svg size 350,262 font "Verdana, 10"'
            write( ua, '(a)' ) 'set output "'//trim(nom)//'_tan.svgt"'
            write( ua, '(a)' ) 'set title "Abbott-Firestone curve"'
            write( ua, '(a)' ) 'set xlabel "%age"'
            write( ua, '(a)' ) 'set ylabel "h"'
            write( ua, '(a)' ) 'plot "' // trim(nom)//'_tan.datt" ' // 'using 1:2  title "Abbott-Firestone"   with lines,\'
            write( ua, '(a)' )     ' "' // trim(nom)//'_tan.datt" ' // 'using 1:3  title "tangent regression" with lines'

         close(unit = ua)

         call run_gnuplot (trim(nom)//'_tan.gplt')

      endif

      results(1:11) = [ real(smrk1, kind=R8)*100/nb_points,  & !  1 smrk1, iso 25178
                        real(smrk2, kind=R8)*100/nb_points,  & !  2 smrk2, iso 25178
                        spk,                                 & !  3 spk , iso 25178
                        abs(svk),                            & !  4 svk , iso 25178
                        spk + off1,                          & !  5 off1, ordonnée de spk
                        abs(svk + off2),                     & !  6 off2, ordonnée de svk
                        sk,                                  & !  7 sk  , iso 25178
                        vec_reg(1) * nb_points / 100.,       & !  8 core slope
                        f,                                   & !  9 adjustment factor (tangent fit)
                        xx(1),                               & ! 10 coeffa_tan        (tangent fit)
                        xx(2) ]                                ! 11 coeffb_tan        (tangent fit)

      deallocate( tab_moy )

      contains

      subroutine cost(me, x, f)
      implicit none
      class(pikaia_class), intent(inout)               :: me
      real(kind=R8)      , intent(in   ), dimension(:) :: x
      real(kind=R8)      , intent(  out)               :: f

         f = 1./(1. + diff_abb_tan(chrom = x(1:4)))

      return
      endsubroutine cost

      real(kind=R8) function diff_abb_tan(chrom)
      implicit none
      real(kind=R8), intent(in), dimension(1:4) :: chrom

         integer(kind=I4) :: i

         ! the height array is standardize to avoid too small values
         diff_abb_tan = sum( [ ( ( tg( real(i - 1, kind = R8) / ( nb_points - 1), chrom ) - tab_moy(i) / mx%si ) ** 2., i = 1, nb_points, 8 ) ] )

      return
      endfunction diff_abb_tan

      real(kind=R8) function tg(xi, chrom)
      implicit none
      real(kind=R8), intent(in)                 :: xi
      real(kind=R8), intent(in), dimension(1:4) :: chrom

         tg = chrom(3) * tan( (PI_R8/2)*( - (1._R8 - chrom(1)) + xi * ( 2._R8 - ( chrom(1) + chrom(2) ) ) ) ) + chrom(4)

      return
      endfunction tg


   endsubroutine abbott_param


endmodule abbott
