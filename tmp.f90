   subroutine correlation_qed_ccsd(wf, E, L_J_ov, d_oo, d_ov, s_vo, u_vovo, γ₁)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), intent(inout) :: E
!
      real(dp), intent(in) :: γ₁
      real(dp), dimension(wf%n_o,wf%n_o), intent(in) :: d_oo
      real(dp), dimension(wf%n_o,wf%n_v), intent(in) :: d_ov
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: s_vo
      real(dp), dimension(wf%eri_t1%n_J,wf%n_o,wf%n_v), intent(in) :: L_J_ov
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: u_vovo
!
      real(dp) :: X1
      real(dp), dimension(:,:), allocatable :: X2
      real(dp), dimension(:,:,:), allocatable :: X3, X4, X5
!
      integer :: i1
!
      real(dp), external :: ddot
!
      X1 = zero
!
      do i1 = 1, wf%n_o
         X1 = X1 + d_oo(i1,i1)
      end do
!
      E = E + two*X1 * γ₁
      call mem%alloc(X2, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, X2, wf%n_o, wf%n_v)
      E = E + two * ddot(wf%n_v*wf%n_o, X2, 1, s_vo, 1)
      call mem%dealloc(X2)
      call mem%alloc(X3, wf%eri_t1%n_J, wf%n_v, wf%n_o)
      call sort_to_132(L_J_ov, X3, wf%eri_t1%n_J, wf%n_o, wf%n_v)
      call mem%alloc(X4, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%eri_t1%n_J, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X3, &
         wf%eri_t1%n_J, &
         u_vovo, &
         wf%n_v*wf%n_o, &
         zero, &
         X4, &
         wf%eri_t1%n_J)
!
      call mem%dealloc(X3)
      call mem%alloc(X5, wf%eri_t1%n_J, wf%n_o, wf%n_v)
      call sort_to_132(X4, X5, wf%eri_t1%n_J, wf%n_v, wf%n_o)
      call mem%dealloc(X4)
      E = E + ddot(wf%n_v*wf%eri_t1%n_J*wf%n_o, X5, 1, L_J_ov, 1)
      call mem%dealloc(X5)
!
   end subroutine correlation_qed_ccsd
