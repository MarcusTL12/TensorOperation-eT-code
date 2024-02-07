let
    func = FortranFunction(("rho_vovo", ["v", "o", "v", "o"]))
    outperms = [[1, 2, 3, 4], [3, 4, 1, 2]]
    cs_vovo = ("cs_vovo", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    g_oooo = ("g_oooo", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    g_ovov = ("g_ovov", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    t_vovo = ("t_vovo", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    ct_vovo = ("ct_vovo", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    s_vovo = ("s_vovo", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    update_code!(func, ein"akbl,kilj->aibj", 1//1, [cs_vovo, g_oooo], outperms)
    update_code!(func, ein"akbl,kcld,cidj->aibj", 1//1, [cs_vovo, g_ovov, t_vovo], outperms)
    update_code!(func, ein"cidj,kcld,akbl->aibj", 1//1, [cs_vovo, g_ovov, t_vovo], outperms)
    update_code!(func, ein"akbl,kcld,cidj->aibj", 1//1, [ct_vovo, g_ovov, s_vovo], outperms)
    update_code!(func, ein"cidj,kcld,akbl->aibj", 1//1, [ct_vovo, g_ovov, s_vovo], outperms)
    finalize_eT_function(func, "jacobian_qed_ccsd_electronic_s2_sym", "qed_ccsd")
end

