let
    func = FortranFunction(("omega_cor_ai", ["v", "o"]))
    L_ovov = ("L_ovov", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    Rᴬ_vo = ("Ra_vo", true)
    Rᴮ_vo = ("Rb_vo", true)
    ζ = ("zeta", true)
    update_code!(func, ein"jbkc,ai,bj,ck,->ai", 4//1, [L_ovov, Rᴬ_vo, Rᴬ_vo, Rᴮ_vo, ζ])    
    update_code!(func, ein"jbkc,aj,bi,ck,->ai", -2//1, [L_ovov, Rᴬ_vo, Rᴬ_vo, Rᴮ_vo, ζ])   
    update_code!(func, ein"jbkc,aj,ck,bi,->ai", -2//1, [L_ovov, Rᴬ_vo, Rᴬ_vo, Rᴮ_vo, ζ])   
    update_code!(func, ein"jbkc,bi,ck,aj,->ai", -2//1, [L_ovov, Rᴬ_vo, Rᴬ_vo, Rᴮ_vo, ζ])   
    update_code!(func, ein"jbkc,bj,ck,ai,->ai", 2//1, [L_ovov, Rᴬ_vo, Rᴬ_vo, Rᴮ_vo, ζ])    
    finalize_eT_function(func, "omega_singles_X3N5", "scc2")
end