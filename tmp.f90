   subroutine jacobian_transpose_qed_ccsd_bilinear_t2_qed_ccsd(wf, sigma_vovo, bs_vo, bs_vovo, bt_vo, bt_vovo, d_oo, d_ov, d_vv)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout) :: sigma_vovo
!
      real(dp), dimension(wf%n_o,wf%n_o), intent(in) :: d_oo
      real(dp), dimension(wf%n_o,wf%n_v), intent(in) :: d_ov
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: bs_vo, bt_vo
      real(dp), dimension(wf%n_v,wf%n_v), intent(in) :: d_vv
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: bs_vovo, bt_vovo
!
      real(dp), dimension(:,:), allocatable :: X1, X3, X5, X6, X7, X8, X9, X11, X12
      real(dp), dimension(:,:,:,:), allocatable :: X2, X4, X10
!
      call mem%alloc(X1, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, X1, wf%n_o, wf%n_v)
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         two, &
         bs_vo, 1, &
         X1, 1, &
         sigma_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X1)
      call mem%alloc(X2, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         bs_vo, 1, &
         d_ov, 1, &
         X2, &
         wf%n_v*wf%n_o)
!
      call add_1423_to_1234(one, X2, sigma_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X2)
      call mem%alloc(X3, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, X3, wf%n_o, wf%n_v)
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         two*wf%s0, &
         bt_vo, 1, &
         X3, 1, &
         sigma_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X3)
      call mem%alloc(X4, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -wf%s0, &
         bt_vo, 1, &
         d_ov, 1, &
         X4, &
         wf%n_v*wf%n_o)
!
      call add_1423_to_1234(one, X4, sigma_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X4)
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         -one, &
         bs_vovo, &
         wf%n_v**2*wf%n_o, &
         d_oo, &
         wf%n_o, &
         one, &
         sigma_vovo, &
         wf%n_v**2*wf%n_o)
!
!
      call dgemm('T', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         d_vv, &
         wf%n_v, &
         bs_vovo, &
         wf%n_v, &
         one, &
         sigma_vovo, &
         wf%n_v)
!
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         -wf%s0, &
         bt_vovo, &
         wf%n_v**2*wf%n_o, &
         d_oo, &
         wf%n_o, &
         one, &
         sigma_vovo, &
         wf%n_v**2*wf%n_o)
!
!
      call dgemm('T', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         wf%s0, &
         d_vv, &
         wf%n_v, &
         bt_vovo, &
         wf%n_v, &
         one, &
         sigma_vovo, &
         wf%n_v)
!
      call mem%alloc(X5, wf%n_o, wf%n_o)
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
         X5, &
         wf%n_o)
!
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         bt_vovo, &
         wf%n_v**2*wf%n_o, &
         X5, &
         wf%n_o, &
         one, &
         sigma_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X5)
      call mem%alloc(X6, wf%n_v, wf%n_v)
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
         X6, &
         wf%n_v)
!
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X6, &
         wf%n_v, &
         bt_vovo, &
         wf%n_v, &
         one, &
         sigma_vovo, &
         wf%n_v)
!
      call mem%dealloc(X6)
      call mem%alloc(X7, wf%n_v, wf%n_o)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         two, &
         bt_vovo, &
         wf%n_v*wf%n_o, &
         wf%s1, 1, &
         zero, &
         X7, 1)
!
      call mem%alloc(X8, wf%n_v, wf%n_o)
      call sort_to_21(d_ov, X8, wf%n_o, wf%n_v)
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X7, 1, &
         X8, 1, &
         sigma_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X7)
      call mem%dealloc(X8)
      call mem%alloc(X9, wf%n_v, wf%n_o)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         bt_vovo, &
         wf%n_v*wf%n_o, &
         wf%s1, 1, &
         zero, &
         X9, 1)
!
      call mem%alloc(X10, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dger(wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X9, 1, &
         d_ov, 1, &
         X10, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X9)
      call add_1423_to_1234(one, X10, sigma_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X10)
      call mem%alloc(X11, wf%n_o, wf%n_o)
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
         X11, &
         wf%n_o)
!
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         wf%s0, &
         bs_vovo, &
         wf%n_v**2*wf%n_o, &
         X11, &
         wf%n_o, &
         one, &
         sigma_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X11)
      call mem%alloc(X12, wf%n_v, wf%n_v)
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
         X12, &
         wf%n_v)
!
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         wf%s0, &
         X12, &
         wf%n_v, &
         bs_vovo, &
         wf%n_v, &
         one, &
         sigma_vovo, &
         wf%n_v)
!
      call mem%dealloc(X12)
!
   end subroutine jacobian_transpose_qed_ccsd_bilinear_t2_qed_ccsd
