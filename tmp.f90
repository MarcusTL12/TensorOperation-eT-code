   subroutine jacobian_qed_ccsd_electronic_s2_sym_qed_ccsd(wf, rho_vovo, cs_vovo, ct_vovo, g_oooo, g_ovov, s_vovo, t_vovo)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout) :: rho_vovo
!
      real(dp), dimension(wf%n_o,wf%n_o,wf%n_o,wf%n_o), intent(in) :: g_oooo
      real(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: g_ovov
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: cs_vovo, ct_vovo, s_vovo, t_vovo
!
      real(dp), dimension(:,:,:,:), allocatable :: X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, X16, X17, X18, X19, X20, X21, X22, X23
!
      call mem%alloc(X1, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(cs_vovo, X1, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X2, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1324(g_oooo, X2, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X3, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2, &
         wf%n_o**2, &
         wf%n_o**2, &
         one, &
         X1, &
         wf%n_v**2, &
         X2, &
         wf%n_o**2, &
         zero, &
         X3, &
         wf%n_v**2)
!
      call mem%dealloc(X1)
      call mem%dealloc(X2)
      call add_1324_to_1234(one, X3, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X3)
      call mem%alloc(X4, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_ovov, X4, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X5, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(t_vovo, X5, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X6, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o**2, &
         wf%n_o**2, &
         wf%n_v**2, &
         one, &
         X4, &
         wf%n_o**2, &
         X5, &
         wf%n_v**2, &
         zero, &
         X6, &
         wf%n_o**2)
!
      call mem%dealloc(X4)
      call mem%dealloc(X5)
      call mem%alloc(X7, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(cs_vovo, X7, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X8, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2, &
         wf%n_o**2, &
         wf%n_o**2, &
         one, &
         X7, &
         wf%n_v**2, &
         X6, &
         wf%n_o**2, &
         zero, &
         X8, &
         wf%n_v**2)
!
      call mem%dealloc(X7)
      call mem%dealloc(X6)
      call add_1324_to_1234(one, X8, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X8)
      call mem%alloc(X9, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(cs_vovo, X9, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X10, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_ovov, X10, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X11, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o**2, &
         wf%n_o**2, &
         wf%n_v**2, &
         one, &
         X9, &
         wf%n_v**2, &
         X10, &
         wf%n_o**2, &
         zero, &
         X11, &
         wf%n_o**2)
!
      call mem%dealloc(X9)
      call mem%dealloc(X10)
      call mem%alloc(X12, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(t_vovo, X12, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X13, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v**2, &
         wf%n_o**2, &
         wf%n_o**2, &
         one, &
         X12, &
         wf%n_v**2, &
         X11, &
         wf%n_o**2, &
         zero, &
         X13, &
         wf%n_v**2)
!
      call mem%dealloc(X11)
      call mem%dealloc(X12)
      call add_1324_to_1234(one, X13, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X13)
      call mem%alloc(X14, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_ovov, X14, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X15, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(s_vovo, X15, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X16, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o**2, &
         wf%n_o**2, &
         wf%n_v**2, &
         one, &
         X14, &
         wf%n_o**2, &
         X15, &
         wf%n_v**2, &
         zero, &
         X16, &
         wf%n_o**2)
!
      call mem%dealloc(X14)
      call mem%dealloc(X15)
      call mem%alloc(X17, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(ct_vovo, X17, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X18, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2, &
         wf%n_o**2, &
         wf%n_o**2, &
         one, &
         X17, &
         wf%n_v**2, &
         X16, &
         wf%n_o**2, &
         zero, &
         X18, &
         wf%n_v**2)
!
      call mem%dealloc(X17)
      call mem%dealloc(X16)
      call add_1324_to_1234(one, X18, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X18)
      call mem%alloc(X19, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(ct_vovo, X19, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X20, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_ovov, X20, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X21, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o**2, &
         wf%n_o**2, &
         wf%n_v**2, &
         one, &
         X19, &
         wf%n_v**2, &
         X20, &
         wf%n_o**2, &
         zero, &
         X21, &
         wf%n_o**2)
!
      call mem%dealloc(X19)
      call mem%dealloc(X20)
      call mem%alloc(X22, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(s_vovo, X22, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X23, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v**2, &
         wf%n_o**2, &
         wf%n_o**2, &
         one, &
         X22, &
         wf%n_v**2, &
         X21, &
         wf%n_o**2, &
         zero, &
         X23, &
         wf%n_v**2)
!
      call mem%dealloc(X21)
      call mem%dealloc(X22)
      call add_1324_to_1234(one, X23, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X23)
!
   end subroutine jacobian_qed_ccsd_electronic_s2_sym_qed_ccsd
