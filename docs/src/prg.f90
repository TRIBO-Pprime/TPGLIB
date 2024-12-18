!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: may, 03 2019
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Asfc. Example of use**
!<  </span>
program test_asfc
use data_arch,     only : I4, R8
use miscellaneous, only : get_unit
use surfile,       only : read_surf, SCALE_SURF
use asfc,          only : calcul_asfc_hermite, indice_fractal
use files,         only : list_files, clean_scratch
implicit none

   type(SCALE_SURF) :: scal_surf                           !! *object [[SCALE_SURF]]*

   real(kind=R8), dimension(:,:), allocatable :: tab_surf  !! *height array*
   real(kind=R8), dimension(1:2)              :: res_asfc  !! *result: asfc, adjustment factor*
   real(kind=R8), dimension(1:3)              :: ind_frac  !! *result: indice fractal*

   character(len = 512), allocatable, dimension(:) :: list_sur
   character(len = 512), allocatable, dimension(:) :: list_sur1
   character(len = 512), allocatable, dimension(:) :: list_sur2

   integer(kind = I4) :: i_g, n_g, n1_g, n2_g, nx, ny

   call clean_scratch()

   call list_files(dir = "sur", list = list_sur1, ext = "sur")
   call list_files(dir = "sur", list = list_sur2, ext = "SUR")

   n1_g = ubound( list_sur1, 1 )
   n2_g = ubound( list_sur2, 1 )

   n_g  = n1_g + n2_g

   allocate( list_sur(1:n_g) )

   list_sur(       1:n1_g) = list_sur1(1:n1_g)
   list_sur(n1_g + 1:n_g ) = list_sur2(1:n2_g)

   do i_g = 1, n_g

      write(*,*) '==============================================='
      write(*,*) trim( list_sur(i_g) )

      call read_surf(nom_fic = trim( list_sur(i_g) ), &  ! IN
                          mu =  0._R8,                &  ! IN , OPT
                       tab_s = tab_surf,              &  ! OUT
                        scal = scal_surf)                ! OUT

      nx = scal_surf%xres
      ny = scal_surf%yres

      call calcul_asfc_hermite(tab_in = tab_surf,     &  !
                                 scal = scal_surf,    &  !
                             asfc_res = res_asfc,     &  !
                                  omp = .true.)          !

      call indice_fractal( tab_in = tab_surf(1:nx, 1:ny),   &  !
                           long   = nx,                     &  !
                           larg   = ny,                     &  !
                           indf   = ind_frac(1:3) )            !

      write(*,*) 'Asfc2 (asfc2 + correlation):             ', res_asfc(1:2)
      write(*,*) 'Box counting (frac. ind. + correlation): ', ind_frac(1), ind_frac(3)

   enddo

   deallocate( list_sur, list_sur1, list_sur2 )

endprogram test_asfc
