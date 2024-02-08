let
    func = FortranFunction(("E", String[]))
    d_oo = ("d_oo", true)
    γ₁ = ("γ₁", true)
    d_ov = ("d_ov", true)
    s_vo = ("s_vo", true)
    L_J_ov = ("L_J_ov", true)
    u_vovo = ("u_vovo", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    update_code!(func, ein"ii,->", 2//1, [d_oo, γ₁])
    update_code!(func, ein"ia,ai->", 2//1, [d_ov, s_vo])
    update_code!(func, ein"zia,zjb,aibj->", 1//1, [L_J_ov, L_J_ov, u_vovo])
    finalize_eT_function(func, "correlation", "qed_ccsd")
end
