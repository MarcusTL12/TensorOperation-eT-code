   subroutine jacobian_qed_ccsd_bilinear_t1(wf, rho_vo, d_vo, d_oo, cs_vo, d_vv, ct_vo, d_ov, cv_vovo, s_vo, cu_vovo, u_vovo)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o), intent(out) :: rho_vo
!
      real(dp), dimension(wf%n_o,wf%n_o), intent(in) :: d_oo
      real(dp), dimension(wf%n_o,wf%n_v), intent(in) :: d_ov
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: d_vo, cs_vo, ct_vo, s_vo
      real(dp), dimension(wf%n_v,wf%n_v), intent(in) :: d_vv
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: cv_vovo, cu_vovo, u_vovo
!
      real(dp) :: X1, X4
!
      real(dp), dimension(:,:), allocatable :: d_ov_21, X2_oo, X3_oo, ct_vo_21
!
      integer :: i1
      real(dp), external :: ddot
!
      call daxpy(wf%n_v*wf%n_o, cγ, d_vo, 1, rho_vo, 1)
      X1 = zero
!
!$omp parallel do schedule(static) private(i1)
      do i1 = 1, wf%n_o
         X1 = X1 + d_oo(i1,i1)
      end do
!$omp end parallel do
!
      call daxpy(wf%n_v*wf%n_o, two*X1, cs_vo, 1, rho_vo, 1)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o, &
                 wf%n_o, &
                 -one, &
                 cs_vo, &
                 wf%n_v, &
                 d_oo, &
                 wf%n_o, &
                 one, &
                 rho_vo, &
                 wf%n_v)
!
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
                 one, &
                 rho_vo, &
                 wf%n_v)
!
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o, &
                 wf%n_o, &
                 -γ₁, &
                 ct_vo, &
                 wf%n_v, &
                 d_oo, &
                 wf%n_o, &
                 one, &
                 rho_vo, &
                 wf%n_v)
!
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o, &
                 wf%n_v, &
                 γ₁, &
                 d_vv, &
                 wf%n_v, &
                 ct_vo, &
                 wf%n_v, &
                 one, &
                 rho_vo, &
                 wf%n_v)
!
      call mem%alloc(d_ov_21, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, d_ov_21, wf%n_o, wf%n_v)
!
      call dgemv('N', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two, &
                 cv_vovo, &
                 wf%n_v*wf%n_o, &
                 d_ov_21, 1, &
                 one, &
                 rho_vo, 1)
!
      call mem%dealloc(d_ov_21)
      call mem%alloc(X2_oo, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 s_vo, &
                 wf%n_v, &
                 d_ov, &
                 wf%n_o, &
                 zero, &
                 X2_oo, &
                 wf%n_o)
!
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_o, &
                 wf%n_o, &
                 -one, &
                 ct_vo, &
                 wf%n_v, &
                 X2_oo, &
                 wf%n_o, &
                 one, &
                 rho_vo, &
                 wf%n_v)
!
      call mem%dealloc(X2_oo)
      call mem%alloc(X3_oo, wf%n_o, wf%n_o)
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
                 X3_oo, &
                 wf%n_o)
!
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 s_vo, &
                 wf%n_v, &
                 X3_oo, &
                 wf%n_o, &
                 one, &
                 rho_vo, &
                 wf%n_v)
!
      call mem%dealloc(X3_oo)
      call mem%alloc(ct_vo_21, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, ct_vo_21, wf%n_v, wf%n_o)
      X4 = ddot(wf%n_v*wf%n_o, d_ov, 1, ct_vo_21, 1)
      call daxpy(wf%n_v*wf%n_o, two*X4, s_vo, 1, rho_vo, 1)
      call mem%alloc(d_ov_21, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, d_ov_21, wf%n_o, wf%n_v)
!
      call dgemv('N', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two*γ₁, &
                 cu_vovo, &
                 wf%n_v*wf%n_o, &
                 d_ov_21, 1, &
                 one, &
                 rho_vo, 1)
!
      call mem%dealloc(d_ov_21)
      call mem%alloc(d_ov_21, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, d_ov_21, wf%n_o, wf%n_v)
!
      call dgemv('N', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 cγ, &
                 u_vovo, &
                 wf%n_v*wf%n_o, &
                 d_ov_21, 1, &
                 one, &
                 rho_vo, 1)
!
      call mem%dealloc(d_ov_21)
!
   end subroutine jacobian_qed_ccsd_bilinear_t1
