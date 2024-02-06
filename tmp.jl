let
    func = FortranFunction(("rho_vovo", ["v", "o", "v", "o"]))
    cs_vovo = ("cs_vovo", true)
    ω = ("wf%qed%frequencies(wf%mode)", false)
    update_code!(func, ein"aibj,->aibj", 2//1, [cs_vovo, ω])
    finalize_eT_function(func, "jacobian_qed_ccsd_photon_s2", "qed_ccsd")
end

