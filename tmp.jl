let
    reset_state()
    rho_vo = "rho_vo" => ("v", "o") => (1:2...,)
    cγ = Sym("cγ")
    d_vo = "d_vo" => ("v", "o") => (1:2...,)
    cs_vo = "cs_vo" => ("v", "o") => (1:2...,)
    d_oo = "d_oo" => ("o", "o") => (1:2...,)
    d_vv = "d_vv" => ("v", "v") => (1:2...,)
    ct_vo = "ct_vo" => ("v", "o") => (1:2...,)
    γ₁ = Sym("γ₁")
    cv_vovo = "cv_vovo" => ("v", "o", "v", "o") => (1:4...,)
    d_ov = "d_ov" => ("o", "v") => (1:2...,)
    s_vo = "s_vo" => ("v", "o") => (1:2...,)
    cu_vovo = "cu_vovo" => ("v", "o", "v", "o") => (1:4...,)
    u_vovo = "u_vovo" => ("v", "o", "v", "o") => (1:4...,)
    @tensor opt=(a=>10χ,b=>10χ,i=>χ,j=>χ) backend=eTbackend begin
        rho_vo[a,i] += cγ*d_vo[a,i]
        rho_vo[a,i] += 2*cs_vo[a,i]*d_oo[j,j]
        rho_vo[a,i] += -cs_vo[a,j]*d_oo[j,i]
        rho_vo[a,i] += cs_vo[b,i]*d_vv[a,b]
        rho_vo[a,i] += -ct_vo[a,j]*d_oo[j,i]*γ₁
        rho_vo[a,i] += ct_vo[b,i]*d_vv[a,b]*γ₁
        rho_vo[a,i] += 2*cv_vovo[a,i,b,j]*d_ov[j,b]
        rho_vo[a,i] += -ct_vo[a,j]*d_ov[j,b]*s_vo[b,i]
        rho_vo[a,i] += -ct_vo[b,i]*d_ov[j,b]*s_vo[a,j]
        rho_vo[a,i] += 2*ct_vo[b,j]*d_ov[j,b]*s_vo[a,i]
        rho_vo[a,i] += 2*cu_vovo[a,i,b,j]*d_ov[j,b]*γ₁
        rho_vo[a,i] += cγ*d_ov[j,b]*u_vovo[a,i,b,j]
    end
    finalize_eT_function("jacobian_qed_ccsd_bilinear_t1", "qed_ccsd")
end

