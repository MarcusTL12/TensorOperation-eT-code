include("src/omeinsum_impl.jl")

s = let
    func = FortranFunction(("omega_vv", ["v", "v"]))
    outperms = [[1, 2], [2, 1]]
    L_J_vv = ("L_J_vv", true)
    γ_vv = ("gamma_vv", true)
    update_code!(func, ein"zac,zbd,cd->ab", 1//1, [L_J_vv, L_J_vv, γ_vv])
    finalize_eT_function(func, "omega_1_ai", "qed_ccsd_2")
end

println(s)
