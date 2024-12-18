!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: may, 3 2024
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Surface gradients and curvatures. Example of use.**
!<  </span>
program test_grad_curv
use grad_curv, only : test_labelize_point, test_label_surf_summits, test_peaks_and_pits_curvatures
implicit none

call test_label_surf_summits()

call test_labelize_point()

call test_peaks_and_pits_curvatures()

stop

!~  expected output:

!~            9          10           9
!~  numerical h      :    1.8512999999999997       ; theoretical h      :    1.8512999999999999
!~  numerical dhdx   :  -0.44879999999999987       ; theoretical dhdx   :  -0.44879999999999992
!~  numerical dhdy   :   -1.3200000000000001       ; theoretical dhdy   :   -1.3200000000000001
!~  numerical d2hdx2 :    1.7952000000000008       ; theoretical d2hdx2 :    1.7951999999999999
!~  numerical d2hdy2 :   -4.1249999999999991       ; theoretical d2hdy2 :   -4.1250000000000000
!~  numerical d2hdxdy:   0.31999999999999984       ; theoretical d2hdxdy:   0.31999999999999995
!~  theoretical xx0:  -0.10000000000000001      numerical xx0:   -9.9999999744180723E-002
!~  theoretical yy0:   0.20000000000000001      numerical yy0:   0.19999999780599384
!~  theoretical point: S numerical point: S
!~  nx, ny =         1551        1555
!~  dx, dy =    1.2899998437499999E-007   1.2899998437499999E-007
!~  S_param_grad:    9.1818702650203546E-002
!~  S_param_curv:    32734.346139266941
!~  peak_curv:   -15965527969.430788
!~  pits_curv:    66825282090.434196

endprogram test_grad_curv
