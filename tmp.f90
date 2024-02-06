   subroutine jacobian_qed_ccsd_photon_s2_qed_ccsd(wf, rho_vovo, cs_vovo)
!!
!! Generated function
!!
      implicit none
!
      class(qed_ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout) :: rho_vovo
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: cs_vovo
!
      call daxpy(wf%n_v**2*wf%n_o**2, two*wf%qed%frequencies(wf%mode), cs_vovo, 1, rho_vovo, 1)
!
   end subroutine jacobian_qed_ccsd_photon_s2_qed_ccsd
