   subroutine rho_test_qed_ccsd(wf, rho, L_ovov, cs_vo, cs_vovo, ct_vo)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), intent(inout) :: rho
!
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: cs_vo, ct_vo
      real(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: L_ovov
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: cs_vovo
!
      real(dp), dimension(:,:), allocatable :: X1, X3, X4, X5
      real(dp), dimension(:,:,:,:), allocatable :: X2
!
      call mem%alloc(X1, wf%n_v, wf%n_o)
      call sort_to_21(wf%fock_ia, X1, wf%n_o, wf%n_v)
      rho = rho + two * ddot(wf%n_v*wf%n_o, X1, 1, cs_vo, 1)
      call mem%dealloc(X1)
      call mem%alloc(X2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_2143(L_ovov, X2, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      rho = rho + two * ddot(wf%n_v**2*wf%n_o**2, X2, 1, cs_vovo, 1)
      call mem%dealloc(X2)
      call mem%alloc(X3, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, X3, wf%n_v, wf%n_o)
      call mem%alloc(X4, wf%n_o, wf%n_v)
!
      call dgemv('T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         two, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X3, 1, &
         zero, &
         X4, 1)
!
      call mem%dealloc(X3)
      call mem%alloc(X5, wf%n_v, wf%n_o)
      call sort_to_21(X4, X5, wf%n_o, wf%n_v)
      call mem%dealloc(X4)
      rho = rho + ddot(wf%n_v*wf%n_o, X5, 1, wf%s1, 1)
      call mem%dealloc(X5)
!
   end subroutine rho_test_qed_ccsd

