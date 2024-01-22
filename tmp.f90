   subroutine jacobian_qed_ccsd_bilinear_s2(wf, rho_vovo, ct_vovo, d_oo, d_vv, cs_vo, cs_vovo, cs, s_vovo, d_ov, v_vovo, ct_vo, t_vovo, cv_vovo)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(out) :: rho_vovo
!
      real(dp), intent(in) :: cs
      real(dp), dimension(wf%n_o,wf%n_o), intent(in) :: d_oo
      real(dp), dimension(wf%n_o,wf%n_v), intent(in) :: d_ov
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: cs_vo, ct_vo
      real(dp), dimension(wf%n_v,wf%n_v), intent(in) :: d_vv
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: ct_vovo, cs_vovo, s_vovo, v_vovo, t_vovo, cv_vovo
!
      real(dp), dimension(:,:), allocatable :: X1_ov, X2_vo, X3_ov, X4_vo, X5_oo, X8_oo, X9_vo, d_ov_21, X11_oo, X12_vo, X14_oo, X15_oo, X16_vo, X17_oo, X18_vo, X19_oo, X22_oo, X24_oo
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_vovo_1243, rho_vovo_2134, cs_vovo_1243, s_vovo_1243, rho_vovo_1342, s_vovo_1342, X6_ovoo, X7_ovoo, X10_vooo, t_vovo_1243, X10_vooo_2341, t_vovo_1342, X13_vooo, X13_vooo_2341, X20_ovoo, ct_vovo_1243, X21_vooo, X21_vooo_2341, X23_ovoo
!

!
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -two, &
                 ct_vovo, &
                 wf%n_v**2*wf%n_o, &
                 d_oo, &
                 wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v**2*wf%n_o)
!
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 two, &
                 d_vv, &
                 wf%n_v, &
                 ct_vovo, &
                 wf%n_v, &
                 one, &
                 rho_vovo, &
                 wf%n_v)
!
      call mem%alloc(X1_ov, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_v, &
                 wf%n_v, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 d_vv, &
                 wf%n_v, &
                 zero, &
                 X1_ov, &
                 wf%n_o)
!
      call mem%alloc(rho_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call zero_array(rho_vovo_1243, wf%n_v**2*wf%n_o**2)
!
      call dger(wf%n_v*wf%n_o, &
                wf%n_v*wf%n_o, &
                one, &
                cs_vo, 1, &
                X1_ov, 1, &
                rho_vovo_1243, 1)
!
      call add_1243_to_1234(one, rho_vovo_1243, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1243)
      call mem%dealloc(X1_ov)
      call mem%alloc(X2_vo, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 d_oo, &
                 wf%n_o, &
                 zero, &
                 X2_vo, &
                 wf%n_v)
!
!
      call dger(wf%n_v*wf%n_o, &
                wf%n_v*wf%n_o, &
                -one, &
                cs_vo, 1, &
                X2_vo, 1, &
                rho_vovo, 1)
!
      call mem%dealloc(X2_vo)
      call mem%alloc(X3_ov, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_v, &
                 wf%n_o, &
                 -one, &
                 d_oo, &
                 wf%n_o, &
                 cs_vo, &
                 wf%n_v, &
                 zero, &
                 X3_ov, &
                 wf%n_o)
!
      call mem%alloc(rho_vovo_2134, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call zero_array(rho_vovo_2134, wf%n_v**2*wf%n_o**2)
!
      call dger(wf%n_v*wf%n_o, &
                wf%n_v*wf%n_o, &
                one, &
                X3_ov, 1, &
                wf%s1, 1, &
                rho_vovo_2134, 1)
!
      call add_2134_to_1234(one, rho_vovo_2134, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_2134)
      call mem%dealloc(X3_ov)
      call mem%alloc(X4_vo, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 d_vv, &
                 wf%n_v, &
                 cs_vo, &
                 wf%n_v, &
                 zero, &
                 X4_vo, &
                 wf%n_v)
!
!
      call dger(wf%n_v*wf%n_o, &
                wf%n_v*wf%n_o, &
                one, &
                X4_vo, 1, &
                wf%s1, 1, &
                rho_vovo, 1)
!
      call mem%dealloc(X4_vo)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 two*-wf%s0, &
                 cs_vovo, &
                 wf%n_v**2*wf%n_o, &
                 d_oo, &
                 wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%alloc(cs_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(cs_vovo, cs_vovo_1243, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_v, &
                 two*wf%s0, &
                 cs_vovo_1243, &
                 wf%n_v**2*wf%n_o, &
                 d_vv, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1243, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(cs_vovo_1243)
      call add_1243_to_1234(one, rho_vovo_1243, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1243)
      call mem%alloc(s_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(s_vovo, s_vovo_1243, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 cs, &
                 d_vv, &
                 wf%n_v, &
                 s_vovo_1243, &
                 wf%n_v**2*wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(s_vovo_1243)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%alloc(s_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(s_vovo, s_vovo_1342, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -cs, &
                 s_vovo_1342, &
                 wf%n_v*wf%n_o**2, &
                 d_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(s_vovo_1342)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%alloc(X5_oo, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 d_ov, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X5_oo, &
                 wf%n_o)
!
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -four, &
                 cs_vovo, &
                 wf%n_v**2*wf%n_o, &
                 X5_oo, &
                 wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X5_oo)
      call mem%alloc(X6_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(cs_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(cs_vovo, cs_vovo_1243, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_o, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 one, &
                 d_ov, &
                 wf%n_o, &
                 cs_vovo_1243, &
                 wf%n_v**2*wf%n_o, &
                 zero, &
                 X6_ovoo, &
                 wf%n_o)
!
      call mem%dealloc(cs_vovo_1243)
      call mem%alloc(rho_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_o, &
                 -four, &
                 X6_ovoo, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1243, &
                 wf%n_v*wf%n_o**2)
!
      call add_1243_to_1234(one, rho_vovo_1243, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1243)
      call mem%dealloc(X6_ovoo)
      call mem%alloc(X7_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(s_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(s_vovo, s_vovo_1243, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_o, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 one, &
                 d_ov, &
                 wf%n_o, &
                 s_vovo_1243, &
                 wf%n_v**2*wf%n_o, &
                 zero, &
                 X7_ovoo, &
                 wf%n_o)
!
      call mem%dealloc(s_vovo_1243)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 -two, &
                 cs_vo, &
                 wf%n_v, &
                 X7_ovoo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X7_ovoo)
      call mem%alloc(X8_oo, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 d_ov, &
                 wf%n_o, &
                 cs_vo, &
                 wf%n_v, &
                 zero, &
                 X8_oo, &
                 wf%n_o)
!
      call mem%alloc(s_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(s_vovo, s_vovo_1342, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -two, &
                 s_vovo_1342, &
                 wf%n_v*wf%n_o**2, &
                 X8_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(s_vovo_1342)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X8_oo)
      call mem%alloc(X9_vo, wf%n_v, wf%n_o)
      call mem%alloc(d_ov_21, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, d_ov_21, wf%n_o, wf%n_v)
!
      call dgemv('N', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 v_vovo, &
                 wf%n_v*wf%n_o, &
                 d_ov_21, 1, &
                 zero, &
                 X9_vo, 1)
!
      call mem%dealloc(d_ov_21)
!
      call dger(wf%n_v*wf%n_o, &
                wf%n_v*wf%n_o, &
                one, &
                cs_vo, 1, &
                X9_vo, 1, &
                rho_vovo, 1)
!
      call mem%dealloc(X9_vo)
      call mem%alloc(X10_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(t_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, t_vovo_1243, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 t_vovo_1243, &
                 wf%n_v**2*wf%n_o, &
                 d_ov, &
                 wf%n_o, &
                 zero, &
                 X10_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(t_vovo_1243)
      call mem%alloc(X10_vooo_2341, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2341(X10_vooo, X10_vooo_2341, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 -one, &
                 ct_vo, &
                 wf%n_v, &
                 X10_vooo_2341, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X10_vooo_2341)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X10_vooo)
      call mem%alloc(X11_oo, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 d_ov, &
                 wf%n_o, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X11_oo, &
                 wf%n_o)
!
      call mem%alloc(t_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, t_vovo_1342, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 t_vovo_1342, &
                 wf%n_v*wf%n_o**2, &
                 X11_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(t_vovo_1342)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X11_oo)
      call mem%alloc(X12_vo, wf%n_v, wf%n_o)
      call mem%alloc(d_ov_21, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, d_ov_21, wf%n_o, wf%n_v)
!
      call dgemv('N', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 cv_vovo, &
                 wf%n_v*wf%n_o, &
                 d_ov_21, 1, &
                 zero, &
                 X12_vo, 1)
!
      call mem%dealloc(d_ov_21)
!
      call dger(wf%n_v*wf%n_o, &
                wf%n_v*wf%n_o, &
                two, &
                X12_vo, 1, &
                wf%s1, 1, &
                rho_vovo, 1)
!
      call mem%dealloc(X12_vo)
      call mem%alloc(X13_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(t_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, t_vovo_1243, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 t_vovo_1243, &
                 wf%n_v**2*wf%n_o, &
                 d_ov, &
                 wf%n_o, &
                 zero, &
                 X13_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(t_vovo_1243)
      call mem%alloc(X13_vooo_2341, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2341(X13_vooo, X13_vooo_2341, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 -wf%s0, &
                 cs_vo, &
                 wf%n_v, &
                 X13_vooo_2341, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X13_vooo_2341)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X13_vooo)
      call mem%alloc(X14_oo, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 d_ov, &
                 wf%n_o, &
                 cs_vo, &
                 wf%n_v, &
                 zero, &
                 X14_oo, &
                 wf%n_o)
!
      call mem%alloc(t_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, t_vovo_1342, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 wf%s0, &
                 t_vovo_1342, &
                 wf%n_v*wf%n_o**2, &
                 X14_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(t_vovo_1342)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X14_oo)
      call mem%alloc(X15_oo, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 d_ov, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X15_oo, &
                 wf%n_o)
!
      call mem%alloc(X16_vo, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o, &
                 wf%n_o, &
                 -one, &
                 ct_vo, &
                 wf%n_v, &
                 X15_oo, &
                 wf%n_o, &
                 zero, &
                 X16_vo, &
                 wf%n_v)
!
      call mem%dealloc(X15_oo)
!
      call dger(wf%n_v*wf%n_o, &
                wf%n_v*wf%n_o, &
                one, &
                X16_vo, 1, &
                wf%s1, 1, &
                rho_vovo, 1)
!
      call mem%dealloc(X16_vo)
      call mem%alloc(X17_oo, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 d_ov, &
                 wf%n_o, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X17_oo, &
                 wf%n_o)
!
      call mem%alloc(X18_vo, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 X17_oo, &
                 wf%n_o, &
                 zero, &
                 X18_vo, &
                 wf%n_v)
!
      call mem%dealloc(X17_oo)
!
      call dger(wf%n_v*wf%n_o, &
                wf%n_v*wf%n_o, &
                one, &
                X18_vo, 1, &
                wf%s1, 1, &
                rho_vovo, 1)
!
      call mem%dealloc(X18_vo)
      call mem%alloc(X19_oo, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 d_ov, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X19_oo, &
                 wf%n_o)
!
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 two*-wf%s0, &
                 ct_vovo, &
                 wf%n_v**2*wf%n_o, &
                 X19_oo, &
                 wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X19_oo)
      call mem%alloc(X20_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(ct_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(ct_vovo, ct_vovo_1243, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_o, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 one, &
                 d_ov, &
                 wf%n_o, &
                 ct_vovo_1243, &
                 wf%n_v**2*wf%n_o, &
                 zero, &
                 X20_ovoo, &
                 wf%n_o)
!
      call mem%dealloc(ct_vovo_1243)
      call mem%alloc(rho_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_o, &
                 two*-wf%s0, &
                 X20_ovoo, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1243, &
                 wf%n_v*wf%n_o**2)
!
      call add_1243_to_1234(one, rho_vovo_1243, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1243)
      call mem%dealloc(X20_ovoo)
      call mem%alloc(X21_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(s_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(s_vovo, s_vovo_1243, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 s_vovo_1243, &
                 wf%n_v**2*wf%n_o, &
                 d_ov, &
                 wf%n_o, &
                 zero, &
                 X21_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(s_vovo_1243)
      call mem%alloc(X21_vooo_2341, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2341(X21_vooo, X21_vooo_2341, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 -wf%s0, &
                 ct_vo, &
                 wf%n_v, &
                 X21_vooo_2341, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X21_vooo_2341)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X21_vooo)
      call mem%alloc(X22_oo, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 d_ov, &
                 wf%n_o, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X22_oo, &
                 wf%n_o)
!
      call mem%alloc(s_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(s_vovo, s_vovo_1342, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 wf%s0, &
                 s_vovo_1342, &
                 wf%n_v*wf%n_o**2, &
                 X22_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(s_vovo_1342)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X22_oo)
      call mem%alloc(X23_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(t_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, t_vovo_1243, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_o, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 one, &
                 d_ov, &
                 wf%n_o, &
                 t_vovo_1243, &
                 wf%n_v**2*wf%n_o, &
                 zero, &
                 X23_ovoo, &
                 wf%n_o)
!
      call mem%dealloc(t_vovo_1243)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 -cs, &
                 wf%s1, &
                 wf%n_v, &
                 X23_ovoo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X23_ovoo)
      call mem%alloc(X24_oo, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 d_ov, &
                 wf%n_o, &
                 zero, &
                 X24_oo, &
                 wf%n_o)
!
      call mem%alloc(t_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, t_vovo_1342, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -cs, &
                 t_vovo_1342, &
                 wf%n_v*wf%n_o**2, &
                 X24_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(t_vovo_1342)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X24_oo)
!
   end subroutine jacobian_qed_ccsd_bilinear_s2
