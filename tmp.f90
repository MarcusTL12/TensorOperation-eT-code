   subroutine jacobian_qed_ccsd_electronic_s2_qed_ccsd(wf, rho_vovo, F_oo, F_ov, F_vv, L_J_vv, L_ooov, L_ovov, L_voov, L_vvov, cs_vo, cs_vovo, ct_vo, ct_vovo, g_oooo, g_ooov, g_ovoo, g_ovov, g_vooo, g_voov, g_vvoo, g_vvov, g_vvvo, s_vovo, t_vovo, v_vovo)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout) :: rho_vovo
!
      real(dp), dimension(wf%n_o,wf%n_o), intent(in) :: F_oo
      real(dp), dimension(wf%n_o,wf%n_v), intent(in) :: F_ov
      real(dp), dimension(wf%n_v,wf%n_o), intent(in) :: cs_vo, ct_vo
      real(dp), dimension(wf%n_v,wf%n_v), intent(in) :: F_vv
      real(dp), dimension(wf%eri_t1%n_J,wf%n_v,wf%n_v), intent(in) :: L_J_vv
      real(dp), dimension(wf%n_o,wf%n_o,wf%n_o,wf%n_o), intent(in) :: g_oooo
      real(dp), dimension(wf%n_o,wf%n_o,wf%n_o,wf%n_v), intent(in) :: L_ooov, g_ooov
      real(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_o), intent(in) :: g_ovoo
      real(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: L_ovov, g_ovov
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: g_vooo
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_o,wf%n_v), intent(in) :: L_voov, g_voov
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: cs_vovo, ct_vovo, s_vovo, t_vovo, v_vovo
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in) :: g_vvoo
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_v), intent(in) :: L_vvov, g_vvov
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o), intent(in) :: g_vvvo
!
      real(dp), dimension(:,:), allocatable :: X10, X11, X14, X28, X29, X30, X31, X32, X33, X37, X38, X39, X40, X43, X44, X115, X117, X122, X124, X137, X139, X152, X154, X158, X159, X161, X162, X163, X164, X165, X167, X168, X169
      real(dp), dimension(:,:,:), allocatable :: X22, X23
      real(dp), dimension(:,:,:,:), allocatable :: X1, X2, X3, X4, X5, X6, X7, X8, X9, X12, X13, X15, X16, X17, X18, X19, X20, X21, X24, X25, X26, X27, X34, X35, X36, X41, X42, X45, X46, X47, X48, X49, X50, X51, X52, X53, X54, X55, X56, X57, X58, X59, X60, X61, X62, X63, X64, X65, X66, X67, X68, X69, X70, X71, X72, X73, X74, X75, X76, X77, X78, X79, X80, X81, X82, X83, X84, X85, X86, X87, X88, X89, X90, X91, X92, X93, X94, X95, X96, X97, X98, X99, X100, X101, X102, X103, X104, X105, X106, X107, X108, X109, X110, X111, X112, X113, X114, X116, X118, X119, X120, X121, X123, X125, X126, X127, X128, X129, X130, X131, X132, X133, X134, X135, X136, X138, X140, X141, X142, X143, X144, X145, X146, X147, X148, X149, X150, X151, X153, X155, X156, X157, X160, X166, X170, X171, X172, X173, X174, X175, X176, X177, X178, X179, X180, X181, X182, X183, X184, X185, X186, X187, X188, X189, X190, X191, X192, X193, X194, X195, X196, X197, X198
!
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         F_vv, &
         wf%n_v, &
         cs_vovo, &
         wf%n_v, &
         one, &
         rho_vovo, &
         wf%n_v)
!
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         -one, &
         cs_vovo, &
         wf%n_v**2*wf%n_o, &
         F_oo, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%alloc(X1, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_4123(g_vooo, X1, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         -one, &
         cs_vo, &
         wf%n_v, &
         X1, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X1)
      call mem%alloc(X2, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1342(g_vvvo, X2, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(X3, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         one, &
         X2, &
         wf%n_v**2*wf%n_o, &
         cs_vo, &
         wf%n_v, &
         zero, &
         X3, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X2)
      call add_1342_to_1234(one, X3, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X3)
      call mem%alloc(X4, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(L_voov, X4, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X4, &
         wf%n_v*wf%n_o, &
         cs_vovo, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X4)
      call mem%alloc(X5, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(cs_vovo, X5, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         X5, &
         wf%n_v*wf%n_o, &
         g_voov, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X5)
      call mem%alloc(X6, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(cs_vovo, X6, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X7, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(g_vvoo, X7, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X8, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         X6, &
         wf%n_v*wf%n_o, &
         X7, &
         wf%n_v*wf%n_o, &
         zero, &
         X8, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X6)
      call mem%dealloc(X7)
      call add_1432_to_1234(one, X8, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X8)
      call mem%alloc(X9, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         -one, &
         F_ov, &
         wf%n_o, &
         t_vovo, &
         wf%n_v, &
         zero, &
         X9, &
         wf%n_o)
!
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         cs_vo, &
         wf%n_v, &
         X9, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X9)
      call mem%alloc(X10, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         F_ov, &
         wf%n_o, &
         cs_vo, &
         wf%n_v, &
         zero, &
         X10, &
         wf%n_o)
!
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         t_vovo, &
         wf%n_v**2*wf%n_o, &
         X10, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X10)
      call mem%alloc(X11, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         wf%s1, &
         wf%n_v, &
         F_ov, &
         wf%n_o, &
         zero, &
         X11, &
         wf%n_o)
!
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         ct_vovo, &
         wf%n_v**2*wf%n_o, &
         X11, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X11)
      call mem%alloc(X12, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         ct_vovo, &
         wf%n_v, &
         F_ov, &
         wf%n_o, &
         zero, &
         X12, &
         wf%n_v*wf%n_o**2)
!
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X12, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X12)
      call mem%alloc(X13, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         s_vovo, &
         wf%n_v, &
         F_ov, &
         wf%n_o, &
         zero, &
         X13, &
         wf%n_v*wf%n_o**2)
!
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X13, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X13)
      call mem%alloc(X14, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         F_ov, &
         wf%n_o, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X14, &
         wf%n_o)
!
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         s_vovo, &
         wf%n_v**2*wf%n_o, &
         X14, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X14)
      call mem%alloc(X15, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         -one, &
         wf%s1, &
         wf%n_v, &
         g_voov, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X15, &
         wf%n_o)
!
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X15, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X15)
      call mem%alloc(X16, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1342(g_vvoo, X16, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X17, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         X16, &
         wf%n_v*wf%n_o**2, &
         wf%s1, &
         wf%n_v, &
         zero, &
         X17, &
         wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X16)
      call mem%alloc(X18, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_3142(X17, X18, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X17)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X18, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X18)
      call mem%alloc(X19, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1342(g_oooo, X19, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X20, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
         wf%n_o**3, &
         wf%n_v, &
         wf%n_o, &
         one, &
         X19, &
         wf%n_o, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X20, &
         wf%n_o**3)
!
      call mem%dealloc(X19)
      call mem%alloc(X21, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(X20, X21, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X20)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X21, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X21)
      call mem%alloc(X22, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%eri_t1%n_J, &
         wf%n_o, &
         wf%n_v, &
         one, &
         L_J_vv, &
         wf%n_v*wf%eri_t1%n_J, &
         wf%s1, &
         wf%n_v, &
         zero, &
         X22, &
         wf%n_v*wf%eri_t1%n_J)
!
      call mem%alloc(X23, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%eri_t1%n_J, &
         wf%n_o, &
         wf%n_v, &
         one, &
         L_J_vv, &
         wf%n_v*wf%eri_t1%n_J, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X23, &
         wf%n_v*wf%eri_t1%n_J)
!
!
      call dgemm('T', 'N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%eri_t1%n_J, &
         one, &
         X23, &
         wf%eri_t1%n_J, &
         X22, &
         wf%eri_t1%n_J, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X22)
      call mem%dealloc(X23)
      call mem%alloc(X24, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1342(g_vvoo, X24, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X25, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         X24, &
         wf%n_v*wf%n_o**2, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X25, &
         wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X24)
      call mem%alloc(X26, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_3142(X25, X26, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X25)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X26, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X26)
      call mem%alloc(X27, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         -one, &
         ct_vo, &
         wf%n_v, &
         g_voov, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X27, &
         wf%n_o)
!
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X27, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X27)
      call mem%alloc(X28, wf%n_o, wf%n_v)
      call sort_to_21(cs_vo, X28, wf%n_v, wf%n_o)
      call mem%alloc(X29, wf%n_v, wf%n_v)
!
      call dgemv('N', &
         wf%n_v**2, &
         wf%n_v*wf%n_o, &
         one, &
         L_vvov, &
         wf%n_v**2, &
         X28, 1, &
         zero, &
         X29, 1)
!
      call mem%dealloc(X28)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X29, &
         wf%n_v, &
         t_vovo, &
         wf%n_v, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X29)
      call mem%alloc(X30, wf%n_o, wf%n_v)
      call sort_to_21(cs_vo, X30, wf%n_v, wf%n_o)
      call mem%alloc(X31, wf%n_o, wf%n_o)
!
      call dgemv('N', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         -one, &
         L_ooov, &
         wf%n_o**2, &
         X30, 1, &
         zero, &
         X31, 1)
!
      call mem%dealloc(X30)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         t_vovo, &
         wf%n_v**2*wf%n_o, &
         X31, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X31)
      call mem%alloc(X32, wf%n_o, wf%n_v)
      call sort_to_21(wf%s1, X32, wf%n_v, wf%n_o)
      call mem%alloc(X33, wf%n_v, wf%n_v)
!
      call dgemv('N', &
         wf%n_v**2, &
         wf%n_v*wf%n_o, &
         one, &
         L_vvov, &
         wf%n_v**2, &
         X32, 1, &
         zero, &
         X33, 1)
!
      call mem%dealloc(X32)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X33, &
         wf%n_v, &
         ct_vovo, &
         wf%n_v, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X33)
      call mem%alloc(X34, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(L_vvov, X34, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X35, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         one, &
         X34, &
         wf%n_v**2*wf%n_o, &
         wf%s1, &
         wf%n_v, &
         zero, &
         X35, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X34)
      call mem%alloc(X36, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X35, X36, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X35)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X36, &
         wf%n_v*wf%n_o, &
         ct_vovo, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X36)
      call mem%alloc(X37, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, X37, wf%n_v, wf%n_o)
      call mem%alloc(X38, wf%n_v, wf%n_v)
!
      call dgemv('N', &
         wf%n_v**2, &
         wf%n_v*wf%n_o, &
         one, &
         L_vvov, &
         wf%n_v**2, &
         X37, 1, &
         zero, &
         X38, 1)
!
      call mem%dealloc(X37)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X38, &
         wf%n_v, &
         s_vovo, &
         wf%n_v, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X38)
      call mem%alloc(X39, wf%n_o, wf%n_v)
      call sort_to_21(wf%s1, X39, wf%n_v, wf%n_o)
      call mem%alloc(X40, wf%n_o, wf%n_o)
!
      call dgemv('N', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         -one, &
         L_ooov, &
         wf%n_o**2, &
         X39, 1, &
         zero, &
         X40, 1)
!
      call mem%dealloc(X39)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         ct_vovo, &
         wf%n_v**2*wf%n_o, &
         X40, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X40)
      call mem%alloc(X41, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(L_ooov, X41, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(X42, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         X41, &
         wf%n_o**2, &
         ct_vovo, &
         wf%n_v*wf%n_o, &
         zero, &
         X42, &
         wf%n_o**2)
!
      call mem%dealloc(X41)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X42, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X42)
      call mem%alloc(X43, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, X43, wf%n_v, wf%n_o)
      call mem%alloc(X44, wf%n_o, wf%n_o)
!
      call dgemv('N', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         -one, &
         L_ooov, &
         wf%n_o**2, &
         X43, 1, &
         zero, &
         X44, 1)
!
      call mem%dealloc(X43)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         s_vovo, &
         wf%n_v**2*wf%n_o, &
         X44, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X44)
      call mem%alloc(X45, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_vvov, X45, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X46, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X46, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X47, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v**2, &
         -one, &
         X46, &
         wf%n_v**2, &
         X45, &
         wf%n_v*wf%n_o, &
         zero, &
         X47, &
         wf%n_o**2)
!
      call mem%dealloc(X45)
      call mem%dealloc(X46)
      call mem%alloc(X48, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(X47, X48, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X47)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         cs_vo, &
         wf%n_v, &
         X48, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X48)
      call mem%alloc(X49, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovoo, X49, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X50, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X49, &
         wf%n_o**2, &
         t_vovo, &
         wf%n_v*wf%n_o, &
         zero, &
         X50, &
         wf%n_o**2)
!
      call mem%dealloc(X49)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         cs_vo, &
         wf%n_v, &
         X50, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X50)
      call mem%alloc(X51, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(t_vovo, X51, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X52, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovoo, X52, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X53, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X52, &
         wf%n_o**2, &
         X51, &
         wf%n_v*wf%n_o, &
         zero, &
         X53, &
         wf%n_o**2)
!
      call mem%dealloc(X51)
      call mem%dealloc(X52)
      call mem%alloc(X54, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X53, X54, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X53)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         cs_vo, &
         wf%n_v, &
         X54, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X54)
      call mem%alloc(X55, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         g_vvov, &
         wf%n_v**2*wf%n_o, &
         cs_vo, &
         wf%n_v, &
         zero, &
         X55, &
         wf%n_v**2*wf%n_o)
!
      call mem%alloc(X56, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(X55, X56, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(X55)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X56, &
         wf%n_v*wf%n_o, &
         t_vovo, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X56)
      call mem%alloc(X57, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v**2*wf%n_o, &
         wf%n_v, &
         -one, &
         cs_vo, &
         wf%n_v, &
         g_vvov, &
         wf%n_v**2*wf%n_o, &
         zero, &
         X57, &
         wf%n_o)
!
      call mem%alloc(X58, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(t_vovo, X58, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X59, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X58, &
         wf%n_v*wf%n_o, &
         X57, &
         wf%n_v*wf%n_o, &
         zero, &
         X59, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X57)
      call mem%dealloc(X58)
      call add_1423_to_1234(one, X59, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X59)
      call mem%alloc(X60, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o**3, &
         wf%n_v, &
         one, &
         cs_vo, &
         wf%n_v, &
         g_ooov, &
         wf%n_o**3, &
         zero, &
         X60, &
         wf%n_o)
!
      call mem%alloc(X61, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1342(X60, X61, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X60)
      call mem%alloc(X62, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(t_vovo, X62, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X63, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v**2, &
         wf%n_o**2, &
         wf%n_o**2, &
         one, &
         X62, &
         wf%n_v**2, &
         X61, &
         wf%n_o**2, &
         zero, &
         X63, &
         wf%n_v**2)
!
      call mem%dealloc(X61)
      call mem%dealloc(X62)
      call add_1324_to_1234(one, X63, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X63)
      call mem%alloc(X64, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(g_ooov, X64, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(X65, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         X64, &
         wf%n_o**2, &
         wf%u_aibj, &
         wf%n_v*wf%n_o, &
         zero, &
         X65, &
         wf%n_o**2)
!
      call mem%dealloc(X64)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         cs_vo, &
         wf%n_v, &
         X65, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X65)
      call mem%alloc(X66, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_vvov, X66, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X67, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         one, &
         X66, &
         wf%n_v**2*wf%n_o, &
         cs_vo, &
         wf%n_v, &
         zero, &
         X67, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X66)
      call mem%alloc(X68, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X67, X68, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X67)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X68, &
         wf%n_v*wf%n_o, &
         wf%u_aibj, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X68)
      call mem%alloc(X69, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(s_vovo, X69, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X70, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_vvov, X70, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X71, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v**2, &
         -one, &
         X69, &
         wf%n_v**2, &
         X70, &
         wf%n_v*wf%n_o, &
         zero, &
         X71, &
         wf%n_o**2)
!
      call mem%dealloc(X69)
      call mem%dealloc(X70)
      call mem%alloc(X72, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(X71, X72, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X71)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X72, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X72)
      call mem%alloc(X73, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovoo, X73, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X74, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X73, &
         wf%n_o**2, &
         s_vovo, &
         wf%n_v*wf%n_o, &
         zero, &
         X74, &
         wf%n_o**2)
!
      call mem%dealloc(X73)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X74, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X74)
      call mem%alloc(X75, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1432(g_ovoo, X75, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X76, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(s_vovo, X76, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X77, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X75, &
         wf%n_o**2, &
         X76, &
         wf%n_v*wf%n_o, &
         zero, &
         X77, &
         wf%n_o**2)
!
      call mem%dealloc(X75)
      call mem%dealloc(X76)
      call mem%alloc(X78, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X77, X78, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X77)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X78, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X78)
      call mem%alloc(X79, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o**3, &
         wf%n_o, &
         wf%n_v, &
         one, &
         g_ooov, &
         wf%n_o**3, &
         wf%s1, &
         wf%n_v, &
         zero, &
         X79, &
         wf%n_o**3)
!
      call mem%alloc(X80, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(ct_vovo, X80, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X81, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1324(X79, X81, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X79)
      call mem%alloc(X82, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2, &
         wf%n_o**2, &
         wf%n_o**2, &
         one, &
         X80, &
         wf%n_v**2, &
         X81, &
         wf%n_o**2, &
         zero, &
         X82, &
         wf%n_v**2)
!
      call mem%dealloc(X80)
      call mem%dealloc(X81)
      call add_1324_to_1234(one, X82, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X82)
      call mem%alloc(X83, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_vvov, X83, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X84, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v**2*wf%n_o, &
         wf%n_v, &
         -one, &
         wf%s1, &
         wf%n_v, &
         X83, &
         wf%n_v**2*wf%n_o, &
         zero, &
         X84, &
         wf%n_o)
!
      call mem%dealloc(X83)
      call mem%alloc(X85, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X85, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X86, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X85, &
         wf%n_v*wf%n_o, &
         X84, &
         wf%n_v*wf%n_o, &
         zero, &
         X86, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X85)
      call mem%dealloc(X84)
      call add_1243_to_1234(one, X86, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X86)
      call mem%alloc(X87, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X87, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X88, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         g_ovoo, &
         wf%n_v*wf%n_o, &
         X87, &
         wf%n_v*wf%n_o, &
         zero, &
         X88, &
         wf%n_o**2)
!
      call mem%dealloc(X87)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X88, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X88)
      call mem%alloc(X89, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v**2*wf%n_o, &
         wf%n_v, &
         -one, &
         wf%s1, &
         wf%n_v, &
         g_vvov, &
         wf%n_v**2*wf%n_o, &
         zero, &
         X89, &
         wf%n_o)
!
      call mem%alloc(X90, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(ct_vovo, X90, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X91, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X90, &
         wf%n_v*wf%n_o, &
         X89, &
         wf%n_v*wf%n_o, &
         zero, &
         X91, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X90)
      call mem%dealloc(X89)
      call add_1423_to_1234(one, X91, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X91)
      call mem%alloc(X92, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X92, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X93, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1432(g_ooov, X93, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(X94, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X93, &
         wf%n_v*wf%n_o, &
         X92, &
         wf%n_v*wf%n_o, &
         zero, &
         X94, &
         wf%n_o**2)
!
      call mem%dealloc(X92)
      call mem%dealloc(X93)
      call mem%alloc(X95, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X94, X95, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X94)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X95, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X95)
      call mem%alloc(X96, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         g_vvov, &
         wf%n_v**2*wf%n_o, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X96, &
         wf%n_v**2*wf%n_o)
!
      call mem%alloc(X97, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(X96, X97, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(X96)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X97, &
         wf%n_v*wf%n_o, &
         s_vovo, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X97)
      call mem%alloc(X98, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v**2*wf%n_o, &
         wf%n_v, &
         -one, &
         ct_vo, &
         wf%n_v, &
         g_vvov, &
         wf%n_v**2*wf%n_o, &
         zero, &
         X98, &
         wf%n_o)
!
      call mem%alloc(X99, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(s_vovo, X99, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X100, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X99, &
         wf%n_v*wf%n_o, &
         X98, &
         wf%n_v*wf%n_o, &
         zero, &
         X100, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X98)
      call mem%dealloc(X99)
      call add_1423_to_1234(one, X100, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X100)
      call mem%alloc(X101, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o**3, &
         wf%n_v, &
         one, &
         ct_vo, &
         wf%n_v, &
         g_ooov, &
         wf%n_o**3, &
         zero, &
         X101, &
         wf%n_o)
!
      call mem%alloc(X102, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1342(X101, X102, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X101)
      call mem%alloc(X103, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(s_vovo, X103, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X104, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v**2, &
         wf%n_o**2, &
         wf%n_o**2, &
         one, &
         X103, &
         wf%n_v**2, &
         X102, &
         wf%n_o**2, &
         zero, &
         X104, &
         wf%n_v**2)
!
      call mem%dealloc(X102)
      call mem%dealloc(X103)
      call add_1324_to_1234(one, X104, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X104)
      call mem%alloc(X105, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_vvov, X105, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X106, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(ct_vovo, X106, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X107, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v**2, &
         -one, &
         X106, &
         wf%n_v**2, &
         X105, &
         wf%n_v*wf%n_o, &
         zero, &
         X107, &
         wf%n_o**2)
!
      call mem%dealloc(X105)
      call mem%dealloc(X106)
      call mem%alloc(X108, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(X107, X108, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X107)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X108, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X108)
      call mem%alloc(X109, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(g_ooov, X109, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(X110, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         X109, &
         wf%n_o**2, &
         v_vovo, &
         wf%n_v*wf%n_o, &
         zero, &
         X110, &
         wf%n_o**2)
!
      call mem%dealloc(X109)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X110, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X110)
      call mem%alloc(X111, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_vvov, X111, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X112, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         one, &
         X111, &
         wf%n_v**2*wf%n_o, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X112, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X111)
      call mem%alloc(X113, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X112, X113, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X112)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X113, &
         wf%n_v*wf%n_o, &
         v_vovo, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X113)
      call mem%alloc(X114, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1432(L_ovov, X114, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X115, wf%n_v, wf%n_v)
!
      call dgemm('T', 'T', &
         wf%n_v, &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         -one, &
         X114, &
         wf%n_v*wf%n_o**2, &
         cs_vovo, &
         wf%n_v, &
         zero, &
         X115, &
         wf%n_v)
!
      call mem%dealloc(X114)
!
      call dgemm('T', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X115, &
         wf%n_v, &
         t_vovo, &
         wf%n_v, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X115)
      call mem%alloc(X116, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(cs_vovo, X116, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X117, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v**2*wf%n_o, &
         -one, &
         X116, &
         wf%n_v**2*wf%n_o, &
         L_ovov, &
         wf%n_o, &
         zero, &
         X117, &
         wf%n_o)
!
      call mem%dealloc(X116)
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         t_vovo, &
         wf%n_v**2*wf%n_o, &
         X117, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X117)
      call mem%alloc(X118, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(cs_vovo, X118, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X119, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X118, &
         wf%n_v*wf%n_o, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         zero, &
         X119, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X118)
      call mem%alloc(X120, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(X119, X120, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X119)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X120, &
         wf%n_v*wf%n_o, &
         wf%u_aibj, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X120)
      call mem%alloc(X121, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1432(L_ovov, X121, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X122, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         -one, &
         ct_vovo, &
         wf%n_v, &
         X121, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X122, &
         wf%n_v)
!
      call mem%dealloc(X121)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X122, &
         wf%n_v, &
         s_vovo, &
         wf%n_v, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X122)
      call mem%alloc(X123, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(ct_vovo, X123, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X124, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v**2*wf%n_o, &
         -one, &
         X123, &
         wf%n_v**2*wf%n_o, &
         L_ovov, &
         wf%n_o, &
         zero, &
         X124, &
         wf%n_o)
!
      call mem%dealloc(X123)
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         s_vovo, &
         wf%n_v**2*wf%n_o, &
         X124, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X124)
      call mem%alloc(X125, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1243(L_ovov, X125, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X126, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         v_vovo, &
         wf%n_v*wf%n_o, &
         X125, &
         wf%n_v*wf%n_o, &
         zero, &
         X126, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X125)
      call mem%alloc(X127, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(X126, X127, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X126)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         ct_vovo, &
         wf%n_v*wf%n_o, &
         X127, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X127)
      call mem%alloc(X128, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovov, X128, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X129, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(cs_vovo, X129, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X130, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X128, &
         wf%n_v*wf%n_o, &
         X129, &
         wf%n_v*wf%n_o, &
         zero, &
         X130, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X128)
      call mem%dealloc(X129)
!
      call dgemm('T', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X130, &
         wf%n_v*wf%n_o, &
         t_vovo, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X130)
      call mem%alloc(X131, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1432(g_ovov, X131, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X132, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(t_vovo, X132, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X133, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X131, &
         wf%n_v*wf%n_o, &
         X132, &
         wf%n_v*wf%n_o, &
         zero, &
         X133, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X131)
      call mem%dealloc(X132)
      call mem%alloc(X134, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(cs_vovo, X134, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X135, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X134, &
         wf%n_v*wf%n_o, &
         X133, &
         wf%n_v*wf%n_o, &
         zero, &
         X135, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X134)
      call mem%dealloc(X133)
      call add_1432_to_1234(one, X135, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X135)
      call mem%alloc(X136, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(wf%u_aibj, X136, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X137, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v**2*wf%n_o, &
         -one, &
         X136, &
         wf%n_v**2*wf%n_o, &
         g_ovov, &
         wf%n_o, &
         zero, &
         X137, &
         wf%n_o)
!
      call mem%dealloc(X136)
!
      call dgemm('N', 'T', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         cs_vovo, &
         wf%n_v**2*wf%n_o, &
         X137, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X137)
      call mem%alloc(X138, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1432(g_ovov, X138, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X139, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         -one, &
         wf%u_aibj, &
         wf%n_v, &
         X138, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X139, &
         wf%n_v)
!
      call mem%dealloc(X138)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X139, &
         wf%n_v, &
         cs_vovo, &
         wf%n_v, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X139)
      call mem%alloc(X140, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1243(g_ovov, X140, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X141, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         X140, &
         wf%n_v*wf%n_o, &
         wf%u_aibj, &
         wf%n_v*wf%n_o, &
         zero, &
         X141, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X140)
      call mem%alloc(X142, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(cs_vovo, X142, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X142, &
         wf%n_v*wf%n_o, &
         X141, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X142)
      call mem%dealloc(X141)
      call mem%alloc(X143, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovov, X143, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X144, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X143, &
         wf%n_v*wf%n_o, &
         s_vovo, &
         wf%n_v*wf%n_o, &
         zero, &
         X144, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X143)
      call mem%alloc(X145, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X145, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X145, &
         wf%n_v*wf%n_o, &
         X144, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X145)
      call mem%dealloc(X144)
      call mem%alloc(X146, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1432(g_ovov, X146, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X147, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X147, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X148, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X146, &
         wf%n_v*wf%n_o, &
         X147, &
         wf%n_v*wf%n_o, &
         zero, &
         X148, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X146)
      call mem%dealloc(X147)
      call mem%alloc(X149, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(s_vovo, X149, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X150, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X148, &
         wf%n_v*wf%n_o, &
         X149, &
         wf%n_v*wf%n_o, &
         zero, &
         X150, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X148)
      call mem%dealloc(X149)
      call add_1432_to_1234(one, X150, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X150)
      call mem%alloc(X151, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(v_vovo, X151, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X152, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v**2*wf%n_o, &
         -one, &
         g_ovov, &
         wf%n_o, &
         X151, &
         wf%n_v**2*wf%n_o, &
         zero, &
         X152, &
         wf%n_o)
!
      call mem%dealloc(X151)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         ct_vovo, &
         wf%n_v**2*wf%n_o, &
         X152, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X152)
      call mem%alloc(X153, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1432(g_ovov, X153, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X154, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         -one, &
         v_vovo, &
         wf%n_v, &
         X153, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X154, &
         wf%n_v)
!
      call mem%dealloc(X153)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X154, &
         wf%n_v, &
         ct_vovo, &
         wf%n_v, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X154)
      call mem%alloc(X155, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X155, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X156, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         X155, &
         wf%n_v*wf%n_o, &
         g_ovov, &
         wf%n_v*wf%n_o, &
         zero, &
         X156, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X155)
      call mem%alloc(X157, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(v_vovo, X157, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X156, &
         wf%n_v*wf%n_o, &
         X157, &
         wf%n_v*wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X157)
      call mem%dealloc(X156)
      call mem%alloc(X158, wf%n_o, wf%n_v)
      call sort_to_21(wf%s1, X158, wf%n_v, wf%n_o)
      call mem%alloc(X159, wf%n_o, wf%n_v)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X158, 1, &
         zero, &
         X159, 1)
!
      call mem%dealloc(X158)
      call mem%alloc(X160, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         wf%n_v, &
         one, &
         t_vovo, &
         wf%n_v, &
         X159, &
         wf%n_o, &
         zero, &
         X160, &
         wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X159)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X160, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X160)
      call mem%alloc(X161, wf%n_o, wf%n_v)
      call sort_to_21(wf%s1, X161, wf%n_v, wf%n_o)
      call mem%alloc(X162, wf%n_o, wf%n_v)
!
      call dgemv('N', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X161, 1, &
         zero, &
         X162, 1)
!
      call mem%dealloc(X161)
      call mem%alloc(X163, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         one, &
         X162, &
         wf%n_o, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X163, &
         wf%n_o)
!
      call mem%dealloc(X162)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         t_vovo, &
         wf%n_v**2*wf%n_o, &
         X163, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X163)
      call mem%alloc(X164, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, X164, wf%n_v, wf%n_o)
      call mem%alloc(X165, wf%n_o, wf%n_v)
!
      call dgemv('T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X164, 1, &
         zero, &
         X165, 1)
!
      call mem%dealloc(X164)
      call mem%alloc(X166, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         X165, &
         wf%n_o, &
         t_vovo, &
         wf%n_v, &
         zero, &
         X166, &
         wf%n_o)
!
      call mem%dealloc(X165)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X166, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X166)
      call mem%alloc(X167, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, X167, wf%n_v, wf%n_o)
      call mem%alloc(X168, wf%n_o, wf%n_v)
!
      call dgemv('T', &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         -one, &
         L_ovov, &
         wf%n_v*wf%n_o, &
         X167, 1, &
         zero, &
         X168, 1)
!
      call mem%dealloc(X167)
      call mem%alloc(X169, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o, &
         wf%n_o, &
         wf%n_v, &
         one, &
         X168, &
         wf%n_o, &
         wf%s1, &
         wf%n_v, &
         zero, &
         X169, &
         wf%n_o)
!
      call mem%dealloc(X168)
!
      call dgemm('N', 'N', &
         wf%n_v**2*wf%n_o, &
         wf%n_o, &
         wf%n_o, &
         one, &
         t_vovo, &
         wf%n_v**2*wf%n_o, &
         X169, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X169)
      call mem%alloc(X170, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X170, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X171, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_ovov, X171, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X172, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_o**2, &
         wf%n_o**2, &
         wf%n_v**2, &
         one, &
         X171, &
         wf%n_o**2, &
         X170, &
         wf%n_v**2, &
         zero, &
         X172, &
         wf%n_o**2)
!
      call mem%dealloc(X170)
      call mem%dealloc(X171)
      call mem%alloc(X173, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
         wf%n_o**3, &
         wf%n_v, &
         wf%n_o, &
         one, &
         X172, &
         wf%n_o, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X173, &
         wf%n_o**3)
!
      call mem%dealloc(X172)
      call mem%alloc(X174, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(X173, X174, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X173)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X174, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X174)
      call mem%alloc(X175, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         wf%n_v, &
         one, &
         g_ovov, &
         wf%n_v*wf%n_o**2, &
         wf%s1, &
         wf%n_v, &
         zero, &
         X175, &
         wf%n_v*wf%n_o**2)
!
      call mem%alloc(X176, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(X175, X176, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(X175)
      call mem%alloc(X177, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X176, &
         wf%n_o**2, &
         t_vovo, &
         wf%n_v*wf%n_o, &
         zero, &
         X177, &
         wf%n_o**2)
!
      call mem%dealloc(X176)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X177, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X177)
      call mem%alloc(X178, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         wf%s1, &
         wf%n_v, &
         g_ovov, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X178, &
         wf%n_o)
!
      call mem%alloc(X179, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(t_vovo, X179, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X180, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         one, &
         X179, &
         wf%n_v*wf%n_o, &
         X178, &
         wf%n_o**2, &
         zero, &
         X180, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X179)
      call mem%dealloc(X178)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X180, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X180)
      call mem%alloc(X181, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         ct_vo, &
         wf%n_v, &
         g_ovov, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X181, &
         wf%n_o)
!
      call mem%alloc(X182, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X181, &
         wf%n_o**2, &
         t_vovo, &
         wf%n_v*wf%n_o, &
         zero, &
         X182, &
         wf%n_o**2)
!
      call mem%dealloc(X181)
      call mem%alloc(X183, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(X182, X183, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X182)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X183, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X183)
      call mem%alloc(X184, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         ct_vo, &
         wf%n_v, &
         g_ovov, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X184, &
         wf%n_o)
!
      call mem%alloc(X185, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(t_vovo, X185, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X186, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_v*wf%n_o, &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         one, &
         X185, &
         wf%n_v*wf%n_o, &
         X184, &
         wf%n_o**2, &
         zero, &
         X186, &
         wf%n_v*wf%n_o)
!
      call mem%dealloc(X184)
      call mem%dealloc(X185)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X186, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X186)
      call mem%alloc(X187, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         one, &
         ct_vo, &
         wf%n_v, &
         g_ovov, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X187, &
         wf%n_o)
!
      call mem%alloc(X188, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(X187, X188, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X187)
      call mem%alloc(X189, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_o**3, &
         wf%n_v, &
         one, &
         wf%s1, &
         wf%n_v, &
         X188, &
         wf%n_o**3, &
         zero, &
         X189, &
         wf%n_o)
!
      call mem%dealloc(X188)
      call mem%alloc(X190, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(t_vovo, X190, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X191, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_v**2, &
         wf%n_o**2, &
         wf%n_o**2, &
         one, &
         X190, &
         wf%n_v**2, &
         X189, &
         wf%n_o**2, &
         zero, &
         X191, &
         wf%n_v**2)
!
      call mem%dealloc(X189)
      call mem%dealloc(X190)
      call add_1342_to_1234(one, X191, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X191)
      call mem%alloc(X192, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o, &
         wf%n_v*wf%n_o**2, &
         wf%n_v, &
         -one, &
         wf%s1, &
         wf%n_v, &
         g_ovov, &
         wf%n_v*wf%n_o**2, &
         zero, &
         X192, &
         wf%n_o)
!
      call mem%alloc(X193, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X192, X193, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X192)
      call mem%alloc(X194, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X193, &
         wf%n_o**2, &
         wf%u_aibj, &
         wf%n_v*wf%n_o, &
         zero, &
         X194, &
         wf%n_o**2)
!
      call mem%dealloc(X193)
      call mem%alloc(X195, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(X194, X195, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X194)
!
      call dgemm('N', 'T', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         ct_vo, &
         wf%n_v, &
         X195, &
         wf%n_v*wf%n_o**2, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X195)
      call mem%alloc(X196, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         wf%n_v, &
         -one, &
         g_ovov, &
         wf%n_v*wf%n_o**2, &
         ct_vo, &
         wf%n_v, &
         zero, &
         X196, &
         wf%n_v*wf%n_o**2)
!
      call mem%alloc(X197, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2134(X196, X197, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(X196)
      call mem%alloc(X198, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
         wf%n_o**2, &
         wf%n_v*wf%n_o, &
         wf%n_v*wf%n_o, &
         one, &
         X197, &
         wf%n_v*wf%n_o, &
         wf%u_aibj, &
         wf%n_v*wf%n_o, &
         zero, &
         X198, &
         wf%n_o**2)
!
      call mem%dealloc(X197)
!
      call dgemm('N', 'N', &
         wf%n_v, &
         wf%n_v*wf%n_o**2, &
         wf%n_o, &
         one, &
         wf%s1, &
         wf%n_v, &
         X198, &
         wf%n_o, &
         one, &
         rho_vovo, &
         wf%n_v)
!
      call mem%dealloc(X198)
!
   end subroutine jacobian_qed_ccsd_electronic_s2_qed_ccsd
