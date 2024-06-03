   subroutine omega_1_ai_qed_ccsd_2(wf, omega_vo, d_oo, d_ov, d_vo, d_vv, v₂_vovo)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd_2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o), intent(inout) :: omega_vo
!
      real(dp), dimension(wf%n_o,wf%n_o), intent(in) :: d_oo
      real(dp), dimension(wf%n_o,wf%n_v), intent(in) :: d_ov
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: d_vo
      real(dp), dimension(wf%n_v,wf%n_v), intent(in) :: d_vv
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: v₂_vovo
!
      real(dp) :: X1
      real(dp), dimension(:,:), allocatable :: X2, X3
!
      integer :: i1
!
      call daxpy(wf%n_v*wf%n_o, two*wf%s0_1, d_vo, 1, omega_vo, 1)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_o, &
         wf%n_v, &
         two, &
         d_vv, &
         wf%n_v, &
         wf%s1_2, &
         wf%n_v, &
         one, &
         omega_vo, &
         wf%n_v)
!
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_o, &
         wf%n_o, &
         -two, &
         wf%s1_2, &
         wf%n_v, &
         d_oo, &
         wf%n_o, &
         one, &
         omega_vo, &
         wf%n_v)
!
      X1 = zero
!
      do i1 = 1, wf%n_o
         X1 = X1 + d_oo(i1,i1)
      end do
!
      call daxpy(wf%n_v*wf%n_o, four*X1, wf%s1_2, 1, omega_vo, 1)
      call mem%alloc(X2, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, X2, wf%n_o, wf%n_v)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         two, &
         v₂_vovo, &
         wf%n_v*wf%n_o, &
         X2, 1, &
         one, &
         omega_vo, 1)
!
      call mem%dealloc(X2)
      call mem%alloc(X3, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, X3, wf%n_o, wf%n_v)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         two*wf%s0_1, &
         wf%u_aibj, &
         wf%n_v*wf%n_o, &
         X3, 1, &
         one, &
         omega_vo, 1)
!
      call mem%dealloc(X3)
!
   end subroutine omega_1_ai_qed_ccsd_2
