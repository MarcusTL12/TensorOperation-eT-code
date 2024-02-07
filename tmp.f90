   subroutine jacobian_qed_ccsd_bilinear_t2_qed_ccsd(wf, rho_vovo, cs, cs_vo, cs_vovo, ct_vo, ct_vovo, cu_vovo, d_oo, d_ov, d_vo, d_vv, s_vovo, t_vovo)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout) :: rho_vovo
!
      real(dp), intent(in) :: cs
      real(dp), dimension(wf%n_o,wf%n_o), intent(in) :: d_oo
      real(dp), dimension(wf%n_o,wf%n_v), intent(in) :: d_ov
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: cs_vo, ct_vo, d_vo
      real(dp), dimension(wf%n_v,wf%n_v), intent(in) :: d_vv
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: cs_vovo, ct_vovo, cu_vovo, s_vovo, t_vovo
!
      real(dp), dimension(:,:), allocatable :: X3, X4, X11, X14, X17, X18, X19, X20, X23, X26, X29, X30, X31, X34
      real(dp), dimension(:,:,:,:), allocatable :: X1, X2, X5, X6, X7, X8, X9, X10, X12, X13, X15, X16, X21, X22, X24, X25, X27, X28, X32, X33, X35, X36
!
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         cs_vo, 1, &
         d_vo, 1, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         -one, &
         cs_vovo, &
         wf%n_v**2*wf%n_o, &
         d_oo, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%alloc(X1, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(cs_vovo, X1, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X2, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         wf%n_v, &
         one, &
         X1, &
         wf%n_v*wf%n_o**2, &
         d_vv, &
         wf%n_v, &
         zero, &
         X2, &
         wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X1)
      call add_1243_to_1234(one, X2, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X2)
      call mem%alloc(X3, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_o, &
         wf%n_o, &
         -one, &
         ct_vo, &
         wf%n_v, &
         d_oo, &
         wf%n_o, &
         zero, &
         X3, &
         wf%n_v)
!
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X3, 1, &
         wf%s1, 1, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X3)
      call mem%alloc(X4, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_o, &
         wf%n_v, &
         one, &
         d_vv, &
         wf%n_v, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X4, &
         wf%n_v)
!
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X4, 1, &
         wf%s1, 1, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X4)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         -wf%s0, &
         ct_vovo, &
         wf%n_v**2*wf%n_o, &
         d_oo, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%alloc(X5, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(ct_vovo, X5, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X6, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         wf%n_v, &
         wf%s0, &
         X5, &
         wf%n_v*wf%n_o**2, &
         d_vv, &
         wf%n_v, &
         zero, &
         X6, &
         wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X5)
      call add_1243_to_1234(one, X6, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X6)
      call mem%alloc(X7, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, X7, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X8, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         cs, &
         d_vv, &
         wf%n_v, &
         X7, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X8, &
         wf%n_v)
!
      call mem%dealloc(X7)
      call add_1342_to_1234(one, X8, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X8)
      call mem%alloc(X9, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X9, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X10, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         -cs, &
         X9, &
         wf%n_v**2*wf%n_o, &
         d_oo, &
         wf%n_o, &
         zero, &
         X10, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X9)
      call add_1342_to_1234(one, X10, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X10)
      call mem%alloc(X11, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v, &
         wf%n_o, &
         -one, &
         cs_vo, &
         wf%n_v, &
         d_ov, &
         wf%n_o, &
         zero, &
         X11, &
         wf%n_v)
!
      call mem%alloc(X12, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, X12, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X13, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X11, &
         wf%n_v, &
         X12, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X13, &
         wf%n_v)
!
      call mem%dealloc(X11)
      call mem%dealloc(X12)
      call add_1342_to_1234(one, X13, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X13)
      call mem%alloc(X14, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         cs_vo, &
         wf%n_v, &
         d_ov, &
         wf%n_o, &
         zero, &
         X14, &
         wf%n_o)
!
      call mem%alloc(X15, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X15, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X16, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         X15, &
         wf%n_v**2*wf%n_o, &
         X14, &
         wf%n_o, &
         zero, &
         X16, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X14)
      call mem%dealloc(X15)
      call add_1342_to_1234(one, X16, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X16)
      call mem%alloc(X17, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, X17, wf%n_o, wf%n_v)
      call mem%alloc(X18, wf%n_v, wf%n_o)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         wf%u_aibj, &
         wf%n_v*wf%n_o, &
         X17, 1, &
         zero, &
         X18, 1)
!
      call mem%dealloc(X17)
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         cs_vo, 1, &
         X18, 1, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X18)
      call mem%alloc(X19, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         d_ov, &
         wf%n_o, &
         wf%s1, &
         wf%n_v, &
         zero, &
         X19, &
         wf%n_o)
!
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         ct_vovo, &
         wf%n_v**2*wf%n_o, &
         X19, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X19)
      call mem%alloc(X20, wf%n_v, wf%n_v)
!
      call dgemm('T', 'T', &
         wf%n_v, &
         wf%n_v, &
         wf%n_o, &
         -one, &
         d_ov, &
         wf%n_o, &
         wf%s1, &
         wf%n_v, &
         zero, &
         X20, &
         wf%n_v)
!
      call mem%alloc(X21, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(ct_vovo, X21, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X22, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         wf%n_v, &
         one, &
         X21, &
         wf%n_v*wf%n_o**2, &
         X20, &
         wf%n_v, &
         zero, &
         X22, &
         wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X21)
      call mem%dealloc(X20)
      call add_1243_to_1234(one, X22, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X22)
      call mem%alloc(X23, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v, &
         wf%n_o, &
         -one, &
         ct_vo, &
         wf%n_v, &
         d_ov, &
         wf%n_o, &
         zero, &
         X23, &
         wf%n_v)
!
      call mem%alloc(X24, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(s_vovo, X24, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X25, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X23, &
         wf%n_v, &
         X24, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X25, &
         wf%n_v)
!
      call mem%dealloc(X23)
      call mem%dealloc(X24)
      call add_1342_to_1234(one, X25, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X25)
      call mem%alloc(X26, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         ct_vo, &
         wf%n_v, &
         d_ov, &
         wf%n_o, &
         zero, &
         X26, &
         wf%n_o)
!
      call mem%alloc(X27, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(s_vovo, X27, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X28, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         X27, &
         wf%n_v**2*wf%n_o, &
         X26, &
         wf%n_o, &
         zero, &
         X28, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X26)
      call mem%dealloc(X27)
      call add_1342_to_1234(one, X28, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X28)
      call mem%alloc(X29, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, X29, wf%n_o, wf%n_v)
      call mem%alloc(X30, wf%n_v, wf%n_o)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         cu_vovo, &
         wf%n_v*wf%n_o, &
         X29, 1, &
         zero, &
         X30, 1)
!
      call mem%dealloc(X29)
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X30, 1, &
         wf%s1, 1, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X30)
      call mem%alloc(X31, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v, &
         wf%n_o, &
         -one, &
         ct_vo, &
         wf%n_v, &
         d_ov, &
         wf%n_o, &
         zero, &
         X31, &
         wf%n_v)
!
      call mem%alloc(X32, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, X32, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X33, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         wf%s0, &
         X31, &
         wf%n_v, &
         X32, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X33, &
         wf%n_v)
!
      call mem%dealloc(X31)
      call mem%dealloc(X32)
      call add_1342_to_1234(one, X33, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X33)
      call mem%alloc(X34, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         ct_vo, &
         wf%n_v, &
         d_ov, &
         wf%n_o, &
         zero, &
         X34, &
         wf%n_o)
!
      call mem%alloc(X35, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X35, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X36, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         wf%s0, &
         X35, &
         wf%n_v**2*wf%n_o, &
         X34, &
         wf%n_o, &
         zero, &
         X36, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X34)
      call mem%dealloc(X35)
      call add_1342_to_1234(one, X36, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X36)
!
   end subroutine jacobian_qed_ccsd_bilinear_t2_qed_ccsd
