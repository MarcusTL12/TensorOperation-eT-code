   subroutine jacobian_qed_ccsd_electronic_s2_qed_ccsd(wf, rho_vovo, F_oo, F_ov, F_vv, L_ooov, L_ovov, L_voov, L_vvov, cs_vo, cs_vovo, ct_vo, ct_vovo, g_oooo, g_ooov, g_ovoo, g_ovov, g_vooo, g_voov, g_vvoo, g_vvov, g_vvvo, g_vvvv, s_vovo, t_vovo, v_vovo)
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
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_v), intent(in) :: g_vvvv
!
      real(dp), dimension(:,:), allocatable :: X12_oo, X14_oo, X21_oo, X38_vv, X39_ov, X41_oo, X42_ov, X44_vv, X45_ov, X50_vv, X51_ov, X53_oo, X54_ov, X59_oo, X60_ov, X132_vv, X135_oo, X142_vv, X145_oo, X160_oo, X162_vv, X177_oo, X179_vv, X186_ov, X187_ov, X191_ov, X192_ov, X193_oo, X195_ov, X196_ov, X200_ov, X201_ov, X202_oo
!
      real(dp), dimension(:,:,:,:), allocatable :: X1_voov, rho_vovo_1342, X2_vvoo, X3_vooo, X4_vvov, X5_vovo, X6_voov, X7_vovo, X8_vovo, rho_vovo_1432, X9_vooo, X10_voov, X11_ooov, X13_vvoo, X15_vooo, X16_voov, X17_ooov, rho_vovo_1243, X18_vooo, X19_voov, X20_ooov, X22_vvoo, X23_ovoo, X24_vooo, X25_oovo, X26_voov, X25_oovo_2341, X27_ovoo, X28_ovoo, X29_oooo, X28_ovoo_2341, X30_ooov, rho_vovo_1324, X31_vovv, X31_vovv_2341, X32_vvvo, rho_vovo_4132, X33_ovoo, X34_voov, X35_oovo, rho_vovo_1423, X36_ovoo, X37_vooo, X40_voov, X43_vvoo, X46_voov, X47_vovo, X48_vovv, X49_vovo, X52_voov, X55_vvoo, X56_ovoo, X57_oovo, X56_ovoo_3412, X58_vooo, X61_vvoo, X62_oovo, X63_vvoo, X64_vovv, X65_ovoo, X66_ovoo, X67_oovo, X66_ovoo_3412, X68_vooo, X69_ovoo, X70_ooov, X71_voov, X69_ovoo_3412, X72_vooo, X73_vovo, X74_vvoo, X75_ovvo, X76_voov, X77_oovv, X78_oooo, X79_vvoo, X80_oooo, X81_ovoo, X82_oovo, X81_ovoo_3412, X83_vooo, X84_vovo, X85_vovv, X86_vovo, X87_oovo, X88_vvoo, X89_vovv, X90_ovoo, X91_ovoo, X92_oovo, X91_ovoo_3412, X93_vooo, X94_ovoo, X95_ooov, X96_voov, X94_ovoo_3412, X97_vooo, X98_oooo, X99_vvoo, X100_oooo, X101_ovvo, X102_vovv, X103_voov, X104_vovo, X105_oovo, X106_voov, rho_vovo_4123, X107_ovvo, X108_voov, X109_oovv, X110_oovo, X111_ovoo, X112_voov, X113_ovoo, rho_vovo_2143, X114_vovo, X115_vvoo, X116_ovvo, X117_voov, X118_oovv, X119_oooo, X120_vvoo, X121_oooo, X122_ovoo, X123_vovv, X124_vvoo, X125_oovo, X126_ovoo, X127_oovo, X126_ovoo_3412, X128_vooo, X129_vovo, X130_vovv, X131_vovo, X133_ovov, X134_voov, X136_vvoo, X137_ovvo, X138_vvoo, X139_vovo, X140_voov, X139_vovo_3412, X141_voov, X143_ovov, X144_voov, X146_vvoo, X147_ovvo, X148_vvoo, X149_vovo, X150_voov, X149_vovo_3412, X151_voov, X152_vovo, X153_ovvo, X154_voov, X155_ovvo, X156_ovvo, X157_voov, X158_ovov, X159_voov, X161_vovo, X163_oovv, X164_voov, X165_voov, X166_vovo, X167_voov, X168_voov, X169_vovo, X170_ovvo, X171_voov, X172_ovvo, X173_ovvo, X174_voov, X175_ovov, X176_voov, X178_vovo, X180_oovv, X181_voov, X182_voov, X183_vovo, X184_voov, X185_voov, X188_vooo, X189_voov, X190_ooov, X194_vvoo, X197_vooo, X198_voov, X199_ooov, X203_vvoo, X204_oooo, X205_oovv, X206_vvoo, X207_ovoo, X208_oooo, X207_ovoo_2341, X209_ooov, X210_oovo, X211_vooo, X212_ooov, X213_ovoo, X213_ovoo_2341, X214_voov, X215_voov, X216_vooo, X217_oovv, X218_vooo, X219_oovv, X220_ovoo, X221_ooov, X220_ovoo_3412, X222_vooo, X223_ovoo, X224_oovv, X225_vvoo, X226_ooov, X227_voov, X225_vvoo_3412, X228_voov, X229_vooo, X229_vooo_2341, X230_ovoo, X231_oovo, X232_vvoo, X233_vooo, X234_oovv, X233_vooo_2341, X235_vovo, X236_ovoo, X237_ovov, X238_vooo, X239_oovv, X240_vvoo, X241_vooo, X240_vvoo_3412, X242_vovo
!

!
      call mem%alloc(X1_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(cs_vovo, X1_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 two, &
                 F_vv, &
                 wf%n_v, &
                 X1_voov, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X1_voov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%alloc(X2_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(cs_vovo, X2_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -two, &
                 X2_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 F_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X2_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%alloc(X3_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1243(g_vooo, X3_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 -one, &
                 cs_vo, &
                 wf%n_v, &
                 X3_vooo, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X3_vooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%alloc(X4_vvov, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1342(g_vvvo, X4_vvov, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X4_vvov, &
                 wf%n_v**2*wf%n_o, &
                 cs_vo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X4_vvov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%alloc(X5_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(L_voov, X5_vovo, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two, &
                 X5_vovo, &
                 wf%n_v*wf%n_o, &
                 cs_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X5_vovo)
      call mem%alloc(X6_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(cs_vovo, X6_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 -two, &
                 X6_voov, &
                 wf%n_v*wf%n_o, &
                 g_voov, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X6_voov)
      call mem%alloc(X7_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(cs_vovo, X7_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X8_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(g_vvoo, X8_vovo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1432, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 -two, &
                 X7_vovo, &
                 wf%n_v*wf%n_o, &
                 X8_vovo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 rho_vovo_1432, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X7_vovo)
      call mem%dealloc(X8_vovo)
      call add_1432_to_1234(one, rho_vovo_1432, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1432)
      call mem%alloc(X9_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X10_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, X10_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 X10_voov, &
                 wf%n_v*wf%n_o**2, &
                 F_ov, &
                 wf%n_o, &
                 zero, &
                 X9_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X10_voov)
      call mem%alloc(X11_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2341(X9_vooo, X11_ooov, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 one, &
                 cs_vo, &
                 wf%n_v, &
                 X11_ooov, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X11_ooov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X9_vooo)
      call mem%alloc(X12_oo, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 cs_vo, &
                 wf%n_v, &
                 F_ov, &
                 wf%n_o, &
                 zero, &
                 X12_oo, &
                 wf%n_o)
!
      call mem%alloc(X13_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X13_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 X13_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X12_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X13_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X12_oo)
      call mem%alloc(X14_oo, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 F_ov, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X14_oo, &
                 wf%n_o)
!
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -two, &
                 ct_vovo, &
                 wf%n_v**2*wf%n_o, &
                 X14_oo, &
                 wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X14_oo)
      call mem%alloc(X15_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X16_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(ct_vovo, X16_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X16_voov, &
                 wf%n_v*wf%n_o**2, &
                 F_ov, &
                 wf%n_o, &
                 zero, &
                 X15_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X16_voov)
      call mem%alloc(X17_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2341(X15_vooo, X17_ooov, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_o, &
                 -two, &
                 X17_ooov, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1243, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X17_ooov)
      call add_1243_to_1234(one, rho_vovo_1243, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1243)
      call mem%dealloc(X15_vooo)
      call mem%alloc(X18_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X19_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(s_vovo, X19_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 X19_voov, &
                 wf%n_v*wf%n_o**2, &
                 F_ov, &
                 wf%n_o, &
                 zero, &
                 X18_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X19_voov)
      call mem%alloc(X20_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2341(X18_vooo, X20_ooov, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X20_ooov, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X20_ooov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X18_vooo)
      call mem%alloc(X21_oo, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 ct_vo, &
                 wf%n_v, &
                 F_ov, &
                 wf%n_o, &
                 zero, &
                 X21_oo, &
                 wf%n_o)
!
      call mem%alloc(X22_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(s_vovo, X22_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 X22_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X21_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X22_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X21_oo)
      call mem%alloc(X23_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 g_voov, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 X23_ovoo, &
                 wf%n_o)
!
      call mem%alloc(X24_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2341(X23_ovoo, X24_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 -one, &
                 ct_vo, &
                 wf%n_v, &
                 X24_vooo, &
                 wf%n_v, &
                 one, &
                 rho_vovo, &
                 wf%n_v)
!
      call mem%dealloc(X24_vooo)
      call mem%dealloc(X23_ovoo)
      call mem%alloc(X25_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X26_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1342(g_vvoo, X26_voov, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X25_oovo_2341, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X26_voov, &
                 wf%n_v*wf%n_o**2, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X25_oovo_2341, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X26_voov)
      call sort_to_4123(X25_oovo_2341, X25_oovo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(X25_oovo_2341)
      call mem%alloc(X27_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_2314(X25_oovo, X27_ovoo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1432, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 -one, &
                 ct_vo, &
                 wf%n_v, &
                 X27_ovoo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1432, &
                 wf%n_v)
!
      call mem%dealloc(X27_ovoo)
      call add_1432_to_1234(one, rho_vovo_1432, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1432)
      call mem%dealloc(X25_oovo)
      call mem%alloc(X28_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X29_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1243(g_oooo, X29_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X28_ovoo_2341, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_o**3, &
                 wf%n_v, &
                 wf%n_o, &
                 one, &
                 X29_oooo, &
                 wf%n_o**3, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X28_ovoo_2341, &
                 wf%n_o**3)
!
      call mem%dealloc(X29_oooo)
      call sort_to_4123(X28_ovoo_2341, X28_ovoo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X28_ovoo_2341)
      call mem%alloc(X30_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1342(X28_ovoo, X30_ooov, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1324, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X30_ooov, &
                 wf%n_o**3, &
                 zero, &
                 rho_vovo_1324, &
                 wf%n_v)
!
      call mem%dealloc(X30_ooov)
      call add_1324_to_1234(one, rho_vovo_1324, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1324)
      call mem%dealloc(X28_ovoo)
      call mem%alloc(X31_vovv, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call mem%alloc(X31_vovv_2341, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v**3, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 g_vvvv, &
                 wf%n_v**3, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X31_vovv_2341, &
                 wf%n_v**3)
!
      call sort_to_4123(X31_vovv_2341, X31_vovv, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(X31_vovv_2341)
      call mem%alloc(X32_vvvo, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1342(X31_vovv, X32_vvvo, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call mem%alloc(rho_vovo_4132, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X32_vvvo, &
                 wf%n_v, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_4132, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X32_vvvo)
      call add_4132_to_1234(one, rho_vovo_4132, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_4132)
      call mem%dealloc(X31_vovv)
      call mem%alloc(X33_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X34_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1342(g_vvoo, X34_voov, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 X34_voov, &
                 wf%n_v*wf%n_o**2, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X33_ovoo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X34_voov)
      call mem%alloc(X35_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1324(X33_ovoo, X35_oovo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1423, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_o, &
                 one, &
                 X35_oovo, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1423, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X35_oovo)
      call add_1423_to_1234(one, rho_vovo_1423, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1423)
      call mem%dealloc(X33_ovoo)
      call mem%alloc(X36_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 g_voov, &
                 wf%n_v*wf%n_o**2, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X36_ovoo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%alloc(X37_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2314(X36_ovoo, X37_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 X37_vooo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X37_vooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X36_ovoo)
      call mem%alloc(X38_vv, wf%n_v, wf%n_v)
      call mem%alloc(X39_ov, wf%n_o, wf%n_v)
      call sort_to_21(cs_vo, X39_ov, wf%n_v, wf%n_o)
!
      call dgemv('N', &
                 wf%n_v**2, &
                 wf%n_v*wf%n_o, &
                 one, &
                 L_vvov, &
                 wf%n_v**2, &
                 X39_ov, 1, &
                 zero, &
                 X38_vv, 1)
!
      call mem%dealloc(X39_ov)
      call mem%alloc(X40_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, X40_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 one, &
                 X38_vv, &
                 wf%n_v, &
                 X40_voov, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X40_voov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X38_vv)
      call mem%alloc(X41_oo, wf%n_o, wf%n_o)
      call mem%alloc(X42_ov, wf%n_o, wf%n_v)
      call sort_to_21(cs_vo, X42_ov, wf%n_v, wf%n_o)
!
      call dgemv('N', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 -one, &
                 L_ooov, &
                 wf%n_o**2, &
                 X42_ov, 1, &
                 zero, &
                 X41_oo, 1)
!
      call mem%dealloc(X42_ov)
      call mem%alloc(X43_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X43_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 X43_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X41_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X43_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X41_oo)
      call mem%alloc(X44_vv, wf%n_v, wf%n_v)
      call mem%alloc(X45_ov, wf%n_o, wf%n_v)
      call sort_to_21(wf%s1, X45_ov, wf%n_v, wf%n_o)
!
      call dgemv('N', &
                 wf%n_v**2, &
                 wf%n_v*wf%n_o, &
                 one, &
                 L_vvov, &
                 wf%n_v**2, &
                 X45_ov, 1, &
                 zero, &
                 X44_vv, 1)
!
      call mem%dealloc(X45_ov)
      call mem%alloc(X46_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(ct_vovo, X46_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 two, &
                 X44_vv, &
                 wf%n_v, &
                 X46_voov, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X46_voov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X44_vv)
      call mem%alloc(X47_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X48_vovv, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(L_vvov, X48_vovv, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X48_vovv, &
                 wf%n_v**2*wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X47_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X48_vovv)
      call mem%alloc(X49_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X47_vovo, X49_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two, &
                 X49_vovo, &
                 wf%n_v*wf%n_o, &
                 ct_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X49_vovo)
      call mem%dealloc(X47_vovo)
      call mem%alloc(X50_vv, wf%n_v, wf%n_v)
      call mem%alloc(X51_ov, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, X51_ov, wf%n_v, wf%n_o)
!
      call dgemv('N', &
                 wf%n_v**2, &
                 wf%n_v*wf%n_o, &
                 one, &
                 L_vvov, &
                 wf%n_v**2, &
                 X51_ov, 1, &
                 zero, &
                 X50_vv, 1)
!
      call mem%dealloc(X51_ov)
      call mem%alloc(X52_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(s_vovo, X52_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 one, &
                 X50_vv, &
                 wf%n_v, &
                 X52_voov, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X52_voov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X50_vv)
      call mem%alloc(X53_oo, wf%n_o, wf%n_o)
      call mem%alloc(X54_ov, wf%n_o, wf%n_v)
      call sort_to_21(wf%s1, X54_ov, wf%n_v, wf%n_o)
!
      call dgemv('N', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 one, &
                 L_ooov, &
                 wf%n_o**2, &
                 X54_ov, 1, &
                 zero, &
                 X53_oo, 1)
!
      call mem%dealloc(X54_ov)
      call mem%alloc(X55_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(ct_vovo, X55_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -two, &
                 X55_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X53_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X55_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X53_oo)
      call mem%alloc(X56_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X57_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(L_ooov, X57_oovo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(X56_ovoo_3412, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X57_oovo, &
                 wf%n_o**2, &
                 ct_vovo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X56_ovoo_3412, &
                 wf%n_o**2)
!
      call mem%dealloc(X57_oovo)
      call sort_to_3412(X56_ovoo_3412, X56_ovoo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X56_ovoo_3412)
      call mem%alloc(X58_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2314(X56_ovoo, X58_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 -two, &
                 wf%s1, &
                 wf%n_v, &
                 X58_vooo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X58_vooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X56_ovoo)
      call mem%alloc(X59_oo, wf%n_o, wf%n_o)
      call mem%alloc(X60_ov, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, X60_ov, wf%n_v, wf%n_o)
!
      call dgemv('N', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 -one, &
                 L_ooov, &
                 wf%n_o**2, &
                 X60_ov, 1, &
                 zero, &
                 X59_oo, 1)
!
      call mem%dealloc(X60_ov)
      call mem%alloc(X61_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(s_vovo, X61_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 X61_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X59_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X61_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X59_oo)
      call mem%alloc(X62_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X63_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(t_vovo, X63_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X64_vovv, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_vvov, X64_vovv, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v**2, &
                 one, &
                 X63_vvoo, &
                 wf%n_v**2, &
                 X64_vovv, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X62_oovo, &
                 wf%n_o**2)
!
      call mem%dealloc(X63_vvoo)
      call mem%dealloc(X64_vovv)
      call mem%alloc(X65_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_2341(X62_oovo, X65_ovoo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1423, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 -one, &
                 cs_vo, &
                 wf%n_v, &
                 X65_ovoo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1423, &
                 wf%n_v)
!
      call mem%dealloc(X65_ovoo)
      call add_1423_to_1234(one, rho_vovo_1423, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1423)
      call mem%dealloc(X62_oovo)
      call mem%alloc(X66_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X67_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovoo, X67_oovo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X66_ovoo_3412, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X67_oovo, &
                 wf%n_o**2, &
                 t_vovo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X66_ovoo_3412, &
                 wf%n_o**2)
!
      call mem%dealloc(X67_oovo)
      call sort_to_3412(X66_ovoo_3412, X66_ovoo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X66_ovoo_3412)
      call mem%alloc(X68_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2314(X66_ovoo, X68_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 one, &
                 cs_vo, &
                 wf%n_v, &
                 X68_vooo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X68_vooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X66_ovoo)
      call mem%alloc(X69_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X70_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1432(g_ovoo, X70_ooov, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X71_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(t_vovo, X71_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X69_ovoo_3412, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X70_ooov, &
                 wf%n_o**2, &
                 X71_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X69_ovoo_3412, &
                 wf%n_o**2)
!
      call mem%dealloc(X70_ooov)
      call mem%dealloc(X71_voov)
      call sort_to_3412(X69_ovoo_3412, X69_ovoo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X69_ovoo_3412)
      call mem%alloc(X72_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2314(X69_ovoo, X72_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1324, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 one, &
                 cs_vo, &
                 wf%n_v, &
                 X72_vooo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1324, &
                 wf%n_v)
!
      call mem%dealloc(X72_vooo)
      call add_1324_to_1234(one, rho_vovo_1324, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1324)
      call mem%dealloc(X69_ovoo)
      call mem%alloc(X73_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
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
                 X73_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%alloc(X74_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(X73_vovo, X74_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 one, &
                 X74_vvoo, &
                 wf%n_v**2, &
                 t_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v**2)
!
      call mem%dealloc(X74_vvoo)
      call mem%dealloc(X73_vovo)
      call mem%alloc(X75_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
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
                 X75_ovvo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%alloc(X76_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(t_vovo, X76_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X77_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1432(X75_ovvo, X77_oovv, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1432, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X76_voov, &
                 wf%n_v*wf%n_o, &
                 X77_oovv, &
                 wf%n_o**2, &
                 zero, &
                 rho_vovo_1432, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X76_voov)
      call mem%dealloc(X77_oovv)
      call add_1432_to_1234(one, rho_vovo_1432, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1432)
      call mem%dealloc(X75_ovvo)
      call mem%alloc(X78_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o**3, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 g_ooov, &
                 wf%n_o**3, &
                 cs_vo, &
                 wf%n_v, &
                 zero, &
                 X78_oooo, &
                 wf%n_o**3)
!
      call mem%alloc(X79_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(t_vovo, X79_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X80_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1423(X78_oooo, X80_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2, &
                 wf%n_o**2, &
                 wf%n_o**2, &
                 one, &
                 X79_vvoo, &
                 wf%n_v**2, &
                 X80_oooo, &
                 wf%n_o**2, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2)
!
      call mem%dealloc(X79_vvoo)
      call mem%dealloc(X80_oooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X78_oooo)
      call mem%alloc(X81_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X82_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(g_ooov, X82_oovo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(X81_ovoo_3412, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X82_oovo, &
                 wf%n_o**2, &
                 wf%u_aibj, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X81_ovoo_3412, &
                 wf%n_o**2)
!
      call mem%dealloc(X82_oovo)
      call sort_to_3412(X81_ovoo_3412, X81_ovoo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X81_ovoo_3412)
      call mem%alloc(X83_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2314(X81_ovoo, X83_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 -one, &
                 cs_vo, &
                 wf%n_v, &
                 X83_vooo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X83_vooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X81_ovoo)
      call mem%alloc(X84_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X85_vovv, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_vvov, X85_vovv, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X85_vovv, &
                 wf%n_v**2*wf%n_o, &
                 cs_vo, &
                 wf%n_v, &
                 zero, &
                 X84_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X85_vovv)
      call mem%alloc(X86_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X84_vovo, X86_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X86_vovo, &
                 wf%n_v*wf%n_o, &
                 wf%u_aibj, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X86_vovo)
      call mem%dealloc(X84_vovo)
      call mem%alloc(X87_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X88_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(s_vovo, X88_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X89_vovv, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_vvov, X89_vovv, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v**2, &
                 one, &
                 X88_vvoo, &
                 wf%n_v**2, &
                 X89_vovv, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X87_oovo, &
                 wf%n_o**2)
!
      call mem%dealloc(X88_vvoo)
      call mem%dealloc(X89_vovv)
      call mem%alloc(X90_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_2341(X87_oovo, X90_ovoo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1423, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 -one, &
                 ct_vo, &
                 wf%n_v, &
                 X90_ovoo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1423, &
                 wf%n_v)
!
      call mem%dealloc(X90_ovoo)
      call add_1423_to_1234(one, rho_vovo_1423, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1423)
      call mem%dealloc(X87_oovo)
      call mem%alloc(X91_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X92_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovoo, X92_oovo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X91_ovoo_3412, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X92_oovo, &
                 wf%n_o**2, &
                 s_vovo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X91_ovoo_3412, &
                 wf%n_o**2)
!
      call mem%dealloc(X92_oovo)
      call sort_to_3412(X91_ovoo_3412, X91_ovoo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X91_ovoo_3412)
      call mem%alloc(X93_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2314(X91_ovoo, X93_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X93_vooo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X93_vooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X91_ovoo)
      call mem%alloc(X94_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X95_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1432(g_ovoo, X95_ooov, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X96_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(s_vovo, X96_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X94_ovoo_3412, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X95_ooov, &
                 wf%n_o**2, &
                 X96_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X94_ovoo_3412, &
                 wf%n_o**2)
!
      call mem%dealloc(X95_ooov)
      call mem%dealloc(X96_voov)
      call sort_to_3412(X94_ovoo_3412, X94_ovoo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X94_ovoo_3412)
      call mem%alloc(X97_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2314(X94_ovoo, X97_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1324, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X97_vooo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1324, &
                 wf%n_v)
!
      call mem%dealloc(X97_vooo)
      call add_1324_to_1234(one, rho_vovo_1324, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1324)
      call mem%dealloc(X94_ovoo)
      call mem%alloc(X98_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
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
                 X98_oooo, &
                 wf%n_o**3)
!
      call mem%alloc(X99_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(ct_vovo, X99_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X100_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1324(X98_oooo, X100_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1324, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2, &
                 wf%n_o**2, &
                 wf%n_o**2, &
                 two, &
                 X99_vvoo, &
                 wf%n_v**2, &
                 X100_oooo, &
                 wf%n_o**2, &
                 zero, &
                 rho_vovo_1324, &
                 wf%n_v**2)
!
      call mem%dealloc(X99_vvoo)
      call mem%dealloc(X100_oooo)
      call add_1324_to_1234(one, rho_vovo_1324, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1324)
      call mem%dealloc(X98_oooo)
      call mem%alloc(X101_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(X102_vovv, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_vvov, X102_vovv, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X102_vovv, &
                 wf%n_v**2*wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X101_ovvo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X102_vovv)
      call mem%alloc(X103_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X103_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X104_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_2431(X101_ovvo, X104_vovo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 -two, &
                 X103_voov, &
                 wf%n_v*wf%n_o, &
                 X104_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X103_voov)
      call mem%dealloc(X104_vovo)
      call mem%dealloc(X101_ovvo)
      call mem%alloc(X105_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X106_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X106_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 g_ovoo, &
                 wf%n_v*wf%n_o, &
                 X106_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X105_oovo, &
                 wf%n_o**2)
!
      call mem%dealloc(X106_voov)
      call mem%alloc(rho_vovo_4123, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_o, &
                 two, &
                 X105_oovo, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_4123, &
                 wf%n_v*wf%n_o**2)
!
      call add_4123_to_1234(one, rho_vovo_4123, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_4123)
      call mem%dealloc(X105_oovo)
      call mem%alloc(X107_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 g_vvov, &
                 wf%n_v**2*wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X107_ovvo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%alloc(X108_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X108_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X109_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1432(X107_ovvo, X109_oovv, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1432, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 -two, &
                 X108_voov, &
                 wf%n_v*wf%n_o, &
                 X109_oovv, &
                 wf%n_o**2, &
                 zero, &
                 rho_vovo_1432, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X108_voov)
      call mem%dealloc(X109_oovv)
      call add_1432_to_1234(one, rho_vovo_1432, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1432)
      call mem%dealloc(X107_ovvo)
      call mem%alloc(X110_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X111_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1423(g_ooov, X111_ovoo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(X112_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X112_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X111_ovoo, &
                 wf%n_v*wf%n_o, &
                 X112_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X110_oovo, &
                 wf%n_o**2)
!
      call mem%dealloc(X111_ovoo)
      call mem%dealloc(X112_voov)
      call mem%alloc(X113_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(X110_oovo, X113_ovoo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_2143, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_o, &
                 two, &
                 X113_ovoo, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_2143, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X113_ovoo)
      call add_2143_to_1234(one, rho_vovo_2143, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_2143)
      call mem%dealloc(X110_oovo)
      call mem%alloc(X114_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
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
                 X114_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%alloc(X115_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(X114_vovo, X115_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 one, &
                 X115_vvoo, &
                 wf%n_v**2, &
                 s_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v**2)
!
      call mem%dealloc(X115_vvoo)
      call mem%dealloc(X114_vovo)
      call mem%alloc(X116_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
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
                 X116_ovvo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%alloc(X117_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(s_vovo, X117_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X118_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1432(X116_ovvo, X118_oovv, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1432, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X117_voov, &
                 wf%n_v*wf%n_o, &
                 X118_oovv, &
                 wf%n_o**2, &
                 zero, &
                 rho_vovo_1432, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X117_voov)
      call mem%dealloc(X118_oovv)
      call add_1432_to_1234(one, rho_vovo_1432, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1432)
      call mem%dealloc(X116_ovvo)
      call mem%alloc(X119_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o**3, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 g_ooov, &
                 wf%n_o**3, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X119_oooo, &
                 wf%n_o**3)
!
      call mem%alloc(X120_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(s_vovo, X120_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X121_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1423(X119_oooo, X121_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2, &
                 wf%n_o**2, &
                 wf%n_o**2, &
                 one, &
                 X120_vvoo, &
                 wf%n_v**2, &
                 X121_oooo, &
                 wf%n_o**2, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2)
!
      call mem%dealloc(X120_vvoo)
      call mem%dealloc(X121_oooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X119_oooo)
      call mem%alloc(X122_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X123_vovv, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_vvov, X123_vovv, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X124_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(ct_vovo, X124_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 wf%n_v**2, &
                 one, &
                 X123_vovv, &
                 wf%n_v*wf%n_o, &
                 X124_vvoo, &
                 wf%n_v**2, &
                 zero, &
                 X122_ovoo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X123_vovv)
      call mem%dealloc(X124_vvoo)
      call mem%alloc(X125_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1324(X122_ovoo, X125_oovo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_o, &
                 -two, &
                 X125_oovo, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1243, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X125_oovo)
      call add_1243_to_1234(one, rho_vovo_1243, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1243)
      call mem%dealloc(X122_ovoo)
      call mem%alloc(X126_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X127_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1243(g_ooov, X127_oovo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(X126_ovoo_3412, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X127_oovo, &
                 wf%n_o**2, &
                 v_vovo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X126_ovoo_3412, &
                 wf%n_o**2)
!
      call mem%dealloc(X127_oovo)
      call sort_to_3412(X126_ovoo_3412, X126_ovoo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X126_ovoo_3412)
      call mem%alloc(X128_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2314(X126_ovoo, X128_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 -one, &
                 ct_vo, &
                 wf%n_v, &
                 X128_vooo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X128_vooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X126_ovoo)
      call mem%alloc(X129_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X130_vovv, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_vvov, X130_vovv, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X130_vovv, &
                 wf%n_v**2*wf%n_o, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X129_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X130_vovv)
      call mem%alloc(X131_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(X129_vovo, X131_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X131_vovo, &
                 wf%n_v*wf%n_o, &
                 v_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X131_vovo)
      call mem%dealloc(X129_vovo)
      call mem%alloc(X132_vv, wf%n_v, wf%n_v)
      call mem%alloc(X133_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1432(L_ovov, X133_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 one, &
                 cs_vovo, &
                 wf%n_v, &
                 X133_ovov, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 X132_vv, &
                 wf%n_v)
!
      call mem%dealloc(X133_ovov)
      call mem%alloc(X134_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, X134_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 -two, &
                 X132_vv, &
                 wf%n_v, &
                 X134_voov, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X134_voov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X132_vv)
      call mem%alloc(X135_oo, wf%n_o, wf%n_o)
      call mem%alloc(X136_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(cs_vovo, X136_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X137_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1243(L_ovov, X137_ovvo, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v**2*wf%n_o, &
                 one, &
                 X136_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X137_ovvo, &
                 wf%n_o, &
                 zero, &
                 X135_oo, &
                 wf%n_o)
!
      call mem%dealloc(X136_vvoo)
      call mem%dealloc(X137_ovvo)
      call mem%alloc(X138_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X138_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -two, &
                 X138_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X135_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X138_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X135_oo)
      call mem%alloc(X139_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X140_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2134(L_ovov, X140_voov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X139_vovo_3412, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X140_voov, &
                 wf%n_v*wf%n_o, &
                 cs_vovo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X139_vovo_3412, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X140_voov)
      call sort_to_3412(X139_vovo_3412, X139_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X139_vovo_3412)
      call mem%alloc(X141_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(X139_vovo, X141_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two, &
                 X141_voov, &
                 wf%n_v*wf%n_o, &
                 wf%u_aibj, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X141_voov)
      call mem%dealloc(X139_vovo)
      call mem%alloc(X142_vv, wf%n_v, wf%n_v)
      call mem%alloc(X143_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1432(L_ovov, X143_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 one, &
                 ct_vovo, &
                 wf%n_v, &
                 X143_ovov, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 X142_vv, &
                 wf%n_v)
!
      call mem%dealloc(X143_ovov)
      call mem%alloc(X144_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(s_vovo, X144_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 -two, &
                 X142_vv, &
                 wf%n_v, &
                 X144_voov, &
                 wf%n_v*wf%n_o**2, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X144_voov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X142_vv)
      call mem%alloc(X145_oo, wf%n_o, wf%n_o)
      call mem%alloc(X146_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(ct_vovo, X146_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X147_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1243(L_ovov, X147_ovvo, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v**2*wf%n_o, &
                 one, &
                 X146_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X147_ovvo, &
                 wf%n_o, &
                 zero, &
                 X145_oo, &
                 wf%n_o)
!
      call mem%dealloc(X146_vvoo)
      call mem%dealloc(X147_ovvo)
      call mem%alloc(X148_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(s_vovo, X148_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -two, &
                 X148_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X145_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X148_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X145_oo)
      call mem%alloc(X149_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X150_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2134(L_ovov, X150_voov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X149_vovo_3412, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X150_voov, &
                 wf%n_v*wf%n_o, &
                 ct_vovo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X149_vovo_3412, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X150_voov)
      call sort_to_3412(X149_vovo_3412, X149_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X149_vovo_3412)
      call mem%alloc(X151_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(X149_vovo, X151_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two, &
                 X151_voov, &
                 wf%n_v*wf%n_o, &
                 v_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X151_voov)
      call mem%dealloc(X149_vovo)
      call mem%alloc(X152_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X153_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovov, X153_ovvo, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X154_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(cs_vovo, X154_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X153_ovvo, &
                 wf%n_v*wf%n_o, &
                 X154_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X152_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X153_ovvo)
      call mem%dealloc(X154_voov)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two, &
                 X152_vovo, &
                 wf%n_v*wf%n_o, &
                 t_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X152_vovo)
      call mem%alloc(X155_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(X156_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovov, X156_ovvo, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X157_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(cs_vovo, X157_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X156_ovvo, &
                 wf%n_v*wf%n_o, &
                 X157_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X155_ovvo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X156_ovvo)
      call mem%dealloc(X157_voov)
      call mem%alloc(X158_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1243(X155_ovvo, X158_ovov, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(X159_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(t_vovo, X159_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1432, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two, &
                 X158_ovov, &
                 wf%n_v*wf%n_o, &
                 X159_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 rho_vovo_1432, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X158_ovov)
      call mem%dealloc(X159_voov)
      call add_1432_to_1234(one, rho_vovo_1432, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1432)
      call mem%dealloc(X155_ovvo)
      call mem%alloc(X160_oo, wf%n_o, wf%n_o)
      call mem%alloc(X161_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(wf%u_aibj, X161_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v**2*wf%n_o, &
                 one, &
                 g_ovov, &
                 wf%n_o, &
                 X161_vovo, &
                 wf%n_v**2*wf%n_o, &
                 zero, &
                 X160_oo, &
                 wf%n_o)
!
      call mem%dealloc(X161_vovo)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -two, &
                 cs_vovo, &
                 wf%n_v**2*wf%n_o, &
                 X160_oo, &
                 wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X160_oo)
      call mem%alloc(X162_vv, wf%n_v, wf%n_v)
      call mem%alloc(X163_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_ovov, X163_oovv, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X164_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(wf%u_aibj, X164_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v, &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 one, &
                 X163_oovv, &
                 wf%n_v*wf%n_o**2, &
                 X164_voov, &
                 wf%n_v, &
                 zero, &
                 X162_vv, &
                 wf%n_v)
!
      call mem%dealloc(X163_oovv)
      call mem%dealloc(X164_voov)
      call mem%alloc(X165_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(cs_vovo, X165_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_v, &
                 -two, &
                 X165_voov, &
                 wf%n_v*wf%n_o**2, &
                 X162_vv, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1243, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X165_voov)
      call add_1243_to_1234(one, rho_vovo_1243, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1243)
      call mem%dealloc(X162_vv)
      call mem%alloc(X166_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X167_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(cs_vovo, X167_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 g_ovov, &
                 wf%n_v*wf%n_o, &
                 X167_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X166_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X167_voov)
      call mem%alloc(X168_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(X166_vovo, X168_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 -two, &
                 X168_voov, &
                 wf%n_v*wf%n_o, &
                 wf%u_aibj, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X168_voov)
      call mem%dealloc(X166_vovo)
      call mem%alloc(X169_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X170_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovov, X170_ovvo, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X171_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X171_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X170_ovvo, &
                 wf%n_v*wf%n_o, &
                 X171_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X169_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X170_ovvo)
      call mem%dealloc(X171_voov)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two, &
                 X169_vovo, &
                 wf%n_v*wf%n_o, &
                 s_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X169_vovo)
      call mem%alloc(X172_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(X173_ovvo, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_to_1423(g_ovov, X173_ovvo, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X174_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X174_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X173_ovvo, &
                 wf%n_v*wf%n_o, &
                 X174_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X172_ovvo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X173_ovvo)
      call mem%dealloc(X174_voov)
      call mem%alloc(X175_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_1243(X172_ovvo, X175_ovov, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(X176_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(s_vovo, X176_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1432, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 two, &
                 X175_ovov, &
                 wf%n_v*wf%n_o, &
                 X176_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 rho_vovo_1432, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X175_ovov)
      call mem%dealloc(X176_voov)
      call add_1432_to_1234(one, rho_vovo_1432, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1432)
      call mem%dealloc(X172_ovvo)
      call mem%alloc(X177_oo, wf%n_o, wf%n_o)
      call mem%alloc(X178_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_1432(v_vovo, X178_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v**2*wf%n_o, &
                 one, &
                 g_ovov, &
                 wf%n_o, &
                 X178_vovo, &
                 wf%n_v**2*wf%n_o, &
                 zero, &
                 X177_oo, &
                 wf%n_o)
!
      call mem%dealloc(X178_vovo)
!
      call dgemm('N', 'N', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 -two, &
                 ct_vovo, &
                 wf%n_v**2*wf%n_o, &
                 X177_oo, &
                 wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X177_oo)
      call mem%alloc(X179_vv, wf%n_v, wf%n_v)
      call mem%alloc(X180_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_ovov, X180_oovv, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X181_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(v_vovo, X181_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v, &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 one, &
                 X180_oovv, &
                 wf%n_v*wf%n_o**2, &
                 X181_voov, &
                 wf%n_v, &
                 zero, &
                 X179_vv, &
                 wf%n_v)
!
      call mem%dealloc(X180_oovv)
      call mem%dealloc(X181_voov)
      call mem%alloc(X182_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(ct_vovo, X182_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1243, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_v, &
                 -two, &
                 X182_voov, &
                 wf%n_v*wf%n_o**2, &
                 X179_vv, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1243, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X182_voov)
      call add_1243_to_1234(one, rho_vovo_1243, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1243)
      call mem%dealloc(X179_vv)
      call mem%alloc(X183_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X184_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(ct_vovo, X184_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 g_ovov, &
                 wf%n_v*wf%n_o, &
                 X184_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X183_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X184_voov)
      call mem%alloc(X185_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(X183_vovo, X185_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 -two, &
                 X185_voov, &
                 wf%n_v*wf%n_o, &
                 v_vovo, &
                 wf%n_v*wf%n_o, &
                 one, &
                 rho_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X185_voov)
      call mem%dealloc(X183_vovo)
      call mem%alloc(X186_ov, wf%n_o, wf%n_v)
      call mem%alloc(X187_ov, wf%n_o, wf%n_v)
      call sort_to_21(wf%s1, X187_ov, wf%n_v, wf%n_o)
!
      call dgemv('N', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 -one, &
                 L_ovov, &
                 wf%n_v*wf%n_o, &
                 X187_ov, 1, &
                 zero, &
                 X186_ov, 1)
!
      call mem%dealloc(X187_ov)
      call mem%alloc(X188_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X189_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, X189_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X189_voov, &
                 wf%n_v*wf%n_o**2, &
                 X186_ov, &
                 wf%n_o, &
                 zero, &
                 X188_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X189_voov)
      call mem%dealloc(X186_ov)
      call mem%alloc(X190_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2341(X188_vooo, X190_ooov, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X190_ooov, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X190_ooov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X188_vooo)
      call mem%alloc(X191_ov, wf%n_o, wf%n_v)
      call mem%alloc(X192_ov, wf%n_o, wf%n_v)
      call sort_to_21(wf%s1, X192_ov, wf%n_v, wf%n_o)
!
      call dgemv('N', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 -one, &
                 L_ovov, &
                 wf%n_v*wf%n_o, &
                 X192_ov, 1, &
                 zero, &
                 X191_ov, 1)
!
      call mem%dealloc(X192_ov)
      call mem%alloc(X193_oo, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X191_ov, &
                 wf%n_o, &
                 zero, &
                 X193_oo, &
                 wf%n_o)
!
      call mem%dealloc(X191_ov)
      call mem%alloc(X194_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X194_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 X194_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X193_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X194_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X193_oo)
      call mem%alloc(X195_ov, wf%n_o, wf%n_v)
      call mem%alloc(X196_ov, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, X196_ov, wf%n_v, wf%n_o)
!
      call dgemv('T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 -one, &
                 L_ovov, &
                 wf%n_v*wf%n_o, &
                 X196_ov, 1, &
                 zero, &
                 X195_ov, 1)
!
      call mem%dealloc(X196_ov)
      call mem%alloc(X197_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X198_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1243(t_vovo, X198_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X198_voov, &
                 wf%n_v*wf%n_o**2, &
                 X195_ov, &
                 wf%n_o, &
                 zero, &
                 X197_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X198_voov)
      call mem%dealloc(X195_ov)
      call mem%alloc(X199_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2341(X197_vooo, X199_ooov, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 X199_ooov, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X199_ooov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X197_vooo)
      call mem%alloc(X200_ov, wf%n_o, wf%n_v)
      call mem%alloc(X201_ov, wf%n_o, wf%n_v)
      call sort_to_21(ct_vo, X201_ov, wf%n_v, wf%n_o)
!
      call dgemv('T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 -one, &
                 L_ovov, &
                 wf%n_v*wf%n_o, &
                 X201_ov, 1, &
                 zero, &
                 X200_ov, 1)
!
      call mem%dealloc(X201_ov)
      call mem%alloc(X202_oo, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 X200_ov, &
                 wf%n_o, &
                 zero, &
                 X202_oo, &
                 wf%n_o)
!
      call mem%dealloc(X200_ov)
      call mem%alloc(X203_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1342(t_vovo, X203_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 wf%n_o, &
                 one, &
                 X203_vvoo, &
                 wf%n_v**2*wf%n_o, &
                 X202_oo, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v**2*wf%n_o)
!
      call mem%dealloc(X203_vvoo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X202_oo)
      call mem%alloc(X204_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X205_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1324(g_ovov, X205_oovv, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X206_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(t_vovo, X206_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_o**2, &
                 wf%n_o**2, &
                 wf%n_v**2, &
                 one, &
                 X205_oovv, &
                 wf%n_o**2, &
                 X206_vvoo, &
                 wf%n_v**2, &
                 zero, &
                 X204_oooo, &
                 wf%n_o**2)
!
      call mem%dealloc(X205_oovv)
      call mem%dealloc(X206_vvoo)
      call mem%alloc(X207_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X208_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1324(X204_oooo, X208_oooo, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X207_ovoo_2341, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'T', &
                 wf%n_o**3, &
                 wf%n_v, &
                 wf%n_o, &
                 one, &
                 X208_oooo, &
                 wf%n_o, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X207_ovoo_2341, &
                 wf%n_o**3)
!
      call mem%dealloc(X208_oooo)
      call sort_to_4123(X207_ovoo_2341, X207_ovoo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X207_ovoo_2341)
      call mem%dealloc(X204_oooo)
      call mem%alloc(X209_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1342(X207_ovoo, X209_ooov, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1324, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X209_ooov, &
                 wf%n_o**3, &
                 zero, &
                 rho_vovo_1324, &
                 wf%n_v)
!
      call mem%dealloc(X209_ooov)
      call add_1324_to_1234(one, rho_vovo_1324, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1324)
      call mem%dealloc(X207_ovoo)
      call mem%alloc(X210_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
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
                 X210_oovo, &
                 wf%n_o)
!
      call mem%alloc(X211_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 one, &
                 t_vovo, &
                 wf%n_v*wf%n_o, &
                 X210_oovo, &
                 wf%n_o**2, &
                 zero, &
                 X211_vooo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X210_oovo)
      call mem%alloc(X212_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_2341(X211_vooo, X212_ooov, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X212_ooov, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X212_ooov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X211_vooo)
      call mem%alloc(X213_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X213_ovoo_2341, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
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
                 X213_ovoo_2341, &
                 wf%n_v*wf%n_o**2)
!
      call sort_to_4123(X213_ovoo_2341, X213_ovoo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X213_ovoo_2341)
      call mem%alloc(X214_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(X215_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(t_vovo, X215_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X216_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2134(X213_ovoo, X216_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X215_voov, &
                 wf%n_v*wf%n_o, &
                 X216_vooo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X214_voov, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X215_voov)
      call mem%dealloc(X216_vooo)
      call mem%dealloc(X213_ovoo)
      call mem%alloc(X217_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_2341(X214_voov, X217_oovv, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%alloc(rho_vovo_1324, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X217_oovv, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1324, &
                 wf%n_v)
!
      call mem%dealloc(X217_oovv)
      call add_1324_to_1234(one, rho_vovo_1324, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1324)
      call mem%dealloc(X214_voov)
      call mem%alloc(X218_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X219_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_ovov, X219_oovv, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X219_oovv, &
                 wf%n_v*wf%n_o**2, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X218_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X219_oovv)
      call mem%alloc(X220_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X221_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_3241(X218_vooo, X221_ooov, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X220_ovoo_3412, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 one, &
                 X221_ooov, &
                 wf%n_o**2, &
                 t_vovo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X220_ovoo_3412, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X221_ooov)
      call sort_to_3412(X220_ovoo_3412, X220_ovoo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X220_ovoo_3412)
      call mem%dealloc(X218_vooo)
      call mem%alloc(X222_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_2314(X220_ovoo, X222_vooo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_o**3, &
                 wf%n_o, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 X222_vooo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X222_vooo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X220_ovoo)
      call mem%alloc(X223_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X224_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_ovov, X224_oovv, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X224_oovv, &
                 wf%n_v*wf%n_o**2, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X223_ovoo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X224_oovv)
      call mem%alloc(X225_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X226_ooov, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_3142(X223_ovoo, X226_ooov, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X227_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1423(t_vovo, X227_voov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X225_vvoo_3412, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_o**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 X226_ooov, &
                 wf%n_o**2, &
                 X227_voov, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X225_vvoo_3412, &
                 wf%n_o**2)
!
      call mem%dealloc(X226_ooov)
      call mem%dealloc(X227_voov)
      call sort_to_3412(X225_vvoo_3412, X225_vvoo, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(X225_vvoo_3412)
      call mem%dealloc(X223_ovoo)
      call mem%alloc(X228_voov, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_to_1342(X225_vvoo, X228_voov, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1423, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('T', 'T', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 wf%n_v, &
                 one, &
                 X228_voov, &
                 wf%n_v, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1423, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X228_voov)
      call add_1423_to_1234(one, rho_vovo_1423, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1423)
      call mem%dealloc(X225_vvoo)
      call mem%alloc(X229_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X229_vooo_2341, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
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
                 X229_vooo_2341, &
                 wf%n_v*wf%n_o**2)
!
      call sort_to_4123(X229_vooo_2341, X229_vooo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X229_vooo_2341)
      call mem%alloc(X230_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X231_oovo, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_2314(X229_vooo, X231_oovo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N', &
                 wf%n_o, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v, &
                 one, &
                 ct_vo, &
                 wf%n_v, &
                 X231_oovo, &
                 wf%n_o, &
                 zero, &
                 X230_ovoo, &
                 wf%n_o)
!
      call mem%dealloc(X231_oovo)
      call mem%dealloc(X229_vooo)
      call mem%alloc(X232_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_1324(t_vovo, X232_vvoo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1324, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v**2, &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 one, &
                 X232_vvoo, &
                 wf%n_v**2, &
                 X230_ovoo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 rho_vovo_1324, &
                 wf%n_v**2)
!
      call mem%dealloc(X232_vvoo)
      call add_1324_to_1234(one, rho_vovo_1324, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1324)
      call mem%dealloc(X230_ovoo)
      call mem%alloc(X233_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X234_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_ovov, X234_oovv, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(X233_vooo_2341, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 one, &
                 X234_oovv, &
                 wf%n_v*wf%n_o**2, &
                 wf%s1, &
                 wf%n_v, &
                 zero, &
                 X233_vooo_2341, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X234_oovv)
      call sort_to_4123(X233_vooo_2341, X233_vooo, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X233_vooo_2341)
      call mem%alloc(X235_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X236_ovoo, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_to_2134(X233_vooo, X236_ovoo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one, &
                 wf%u_aibj, &
                 wf%n_v*wf%n_o, &
                 X236_ovoo, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X235_vovo, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X236_ovoo)
      call mem%dealloc(X233_vooo)
      call mem%alloc(X237_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_to_2341(X235_vovo, X237_ovov, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v**2*wf%n_o, &
                 wf%n_o, &
                 -one, &
                 ct_vo, &
                 wf%n_v, &
                 X237_ovov, &
                 wf%n_o, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X237_ovov)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X235_vovo)
      call mem%alloc(X238_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X239_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_to_1342(g_ovov, X239_oovv, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 X239_oovv, &
                 wf%n_v*wf%n_o**2, &
                 ct_vo, &
                 wf%n_v, &
                 zero, &
                 X238_vooo, &
                 wf%n_v*wf%n_o**2)
!
      call mem%dealloc(X239_oovv)
      call mem%alloc(X240_vvoo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(X241_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_to_1432(X238_vooo, X241_vooo, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(X240_vvoo_3412, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T', &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_o**2, &
                 one, &
                 X241_vooo, &
                 wf%n_v*wf%n_o, &
                 wf%u_aibj, &
                 wf%n_v*wf%n_o, &
                 zero, &
                 X240_vvoo_3412, &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X241_vooo)
      call sort_to_3412(X240_vvoo_3412, X240_vvoo, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(X240_vvoo_3412)
      call mem%dealloc(X238_vooo)
      call mem%alloc(X242_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_to_2314(X240_vvoo, X242_vovo, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_vovo_1342, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v, &
                 wf%n_v*wf%n_o**2, &
                 wf%n_o, &
                 one, &
                 wf%s1, &
                 wf%n_v, &
                 X242_vovo, &
                 wf%n_v, &
                 zero, &
                 rho_vovo_1342, &
                 wf%n_v)
!
      call mem%dealloc(X242_vovo)
      call add_1342_to_1234(one, rho_vovo_1342, rho_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_vovo_1342)
      call mem%dealloc(X240_vvoo)
!
   end subroutine jacobian_qed_ccsd_electronic_s2_qed_ccsd
