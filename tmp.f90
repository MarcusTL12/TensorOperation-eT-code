   subroutine test_routine_ccsd(wf, X, A_vovo, B_oo)
!!
!! Generated function
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_v), intent(inout) :: X
!
      real(dp), dimension(wf%n_o,wf%n_o), intent(in) :: B_oo
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: A_vovo
!
      real(dp), dimension(:,:,:,:), allocatable :: X1
!
      call mem%alloc(X1, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(A_vovo, X1, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemv('N', &
         wf%n_v**2, &
         wf%n_o**2, &
         one, &
         X1, &
         wf%n_v**2, &
         B_oo, 1, &
         one, &
         X, 1)
!
      call mem%dealloc(X1)
!
   end subroutine test_routine_ccsd
