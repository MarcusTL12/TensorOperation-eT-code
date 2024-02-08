   subroutine omega_singles_X3N5_scc2(wf, omega_cor_ai, L_ovov, Ra_vo, Rb_vo, zeta)
!!
!! Generated function
!!
      implicit none
!
      class(scc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o), intent(inout) :: omega_cor_ai
!
      real(dp), intent(in) :: zeta
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: Ra_vo, Rb_vo
      real(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: L_ovov
!
      real(dp) :: X4, X17
      real(dp), dimension(:,:), allocatable :: X1, X2, X3, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, X16
!
      real(dp), external :: ddot
!
      call mem%alloc(X1, wf%n_o, wf%n_v)
      call sort_to_21(Ra_vo, X1, wf%n_v, wf%n_o)
      call mem%alloc(X2, wf%n_o, wf%n_v)
!
      call dgemv('T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         four, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X1, 1, &
         zero, &
         X2, 1)
!
      call mem%dealloc(X1)
      call mem%alloc(X3, wf%n_v, wf%n_o)
      call sort_to_21(X2, X3, wf%n_o, wf%n_v)
      call mem%dealloc(X2)
      X4 = ddot(wf%n_v*wf%n_o, X3, 1, Rb_vo, 1)
      call mem%dealloc(X3)
      call daxpy(wf%n_v*wf%n_o, zeta, Ra_vo, 1, omega_cor_ai, 1)
      call mem%alloc(X5, wf%n_o, wf%n_v)
      call sort_to_21(Rb_vo, X5, wf%n_v, wf%n_o)
      call mem%alloc(X6, wf%n_o, wf%n_v)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -two, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X5, 1, &
         zero, &
         X6, 1)
!
      call mem%dealloc(X5)
      call mem%alloc(X7, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         one, &
         X6, &
         wf%n_o, &
         Ra_vo, &
         wf%n_v, &
         zero, &
         X7, &
         wf%n_o)
!
      call mem%dealloc(X6)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_o, &
         wf%n_o, &
         zeta, &
         Ra_vo, &
         wf%n_v, &
         X7, &
         wf%n_o, &
         one, &
         omega_cor_ai, &
         wf%n_v)
!
      call mem%dealloc(X7)
      call mem%alloc(X8, wf%n_o, wf%n_v)
      call sort_to_21(Ra_vo, X8, wf%n_v, wf%n_o)
      call mem%alloc(X9, wf%n_o, wf%n_v)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -two, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X8, 1, &
         zero, &
         X9, 1)
!
      call mem%dealloc(X8)
      call mem%alloc(X10, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         one, &
         X9, &
         wf%n_o, &
         Rb_vo, &
         wf%n_v, &
         zero, &
         X10, &
         wf%n_o)
!
      call mem%dealloc(X9)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_o, &
         wf%n_o, &
         zeta, &
         Ra_vo, &
         wf%n_v, &
         X10, &
         wf%n_o, &
         one, &
         omega_cor_ai, &
         wf%n_v)
!
      call mem%dealloc(X10)
      call mem%alloc(X11, wf%n_o, wf%n_v)
      call sort_to_21(Ra_vo, X11, wf%n_v, wf%n_o)
      call mem%alloc(X12, wf%n_o, wf%n_v)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -two, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X11, 1, &
         zero, &
         X12, 1)
!
      call mem%dealloc(X11)
      call mem%alloc(X13, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         one, &
         X12, &
         wf%n_o, &
         Ra_vo, &
         wf%n_v, &
         zero, &
         X13, &
         wf%n_o)
!
      call mem%dealloc(X12)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_o, &
         wf%n_o, &
         zeta, &
         Rb_vo, &
         wf%n_v, &
         X13, &
         wf%n_o, &
         one, &
         omega_cor_ai, &
         wf%n_v)
!
      call mem%dealloc(X13)
      call mem%alloc(X14, wf%n_o, wf%n_v)
      call sort_to_21(Ra_vo, X14, wf%n_v, wf%n_o)
      call mem%alloc(X15, wf%n_o, wf%n_v)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         two, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X14, 1, &
         zero, &
         X15, 1)
!
      call mem%dealloc(X14)
      call mem%alloc(X16, wf%n_v, wf%n_o)
      call sort_to_21(X15, X16, wf%n_o, wf%n_v)
      call mem%dealloc(X15)
      X17 = ddot(wf%n_v*wf%n_o, X16, 1, Ra_vo, 1)
      call mem%dealloc(X16)
      call daxpy(wf%n_v*wf%n_o, zeta, Rb_vo, 1, omega_cor_ai, 1)
!
   end subroutine omega_singles_X3N5_scc2
