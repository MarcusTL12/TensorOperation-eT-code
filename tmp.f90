   subroutine rho_test_qed_ccsd(wf, rho_vovo, cs_vo, ct_vo, d_ov, t_vovo)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout) :: rho_vovo
!
      real(dp), dimension(wf%n_o,wf%n_v), intent(in) :: d_ov
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: cs_vo, ct_vo
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_vovo
!
      real(dp), dimension(:,:), allocatable :: X1, X2, X3, X4
!
      call mem%alloc(X1, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, X1, wf%n_o, wf%n_v)
      call mem%alloc(X2, wf%n_v, wf%n_o)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         wf%u_aibj, &
         wf%n_v*wf%n_o, &
         X1, 1, &
         zero, &
         X2, 1)
!
      call mem%dealloc(X1)
!
      call dger(wf%n_v*wf%n_o,
         wf%n_v*wf%n_o, &
         one, &
         cs_vo, 1, &
         X2, 1, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X2)
      call mem%alloc(X3, wf%n_v, wf%n_v)
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
         X3, &
         wf%n_v)
!
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         wf%s0, &
         X3, &
         wf%n_v, &
         t_vovo, &
         wf%n_v, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X3)
      call mem%alloc(X4, wf%n_o, wf%n_o)
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
         X4, &
         wf%n_o)
!
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         wf%s0, &
         t_vovo, &
         wf%n_v**2*wf%n_o, &
         X4, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X4)
!
   end subroutine rho_test_qed_ccsd

