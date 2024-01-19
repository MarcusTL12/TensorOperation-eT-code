let
    reset_state()
    rho_vovo = output_tensor("rho_vovo", ("v", "o", "v", "o"))
    ct_vovo = input_tensor("ct_vovo", ("v", "o", "v", "o"))
    d_oo = input_tensor("d_oo", ("o", "o"))
    d_vv = input_tensor("d_vv", ("v", "v"))
    cs_vo = input_tensor("cs_vo", ("v", "o"))
    s_vo = noio_tensor("wf%s1", ("v", "o"))
    cs_vovo = input_tensor("cs_vovo", ("v", "o", "v", "o"))
    γ₁ = Sym("wf%s0")
    cγ = input_scalar("cs")
    s_vovo = input_tensor("s_vovo", ("v", "o", "v", "o"))
    d_ov = input_tensor("d_ov", ("o", "v"))
    v_vovo = input_tensor("v_vovo", ("v", "o", "v", "o"))
    ct_vo = input_tensor("ct_vo", ("v", "o"))
    t_vovo = input_tensor("t_vovo", ("v", "o", "v", "o"))
    cv_vovo = input_tensor("cv_vovo", ("v", "o", "v", "o"))
    @tensor opt=(a=>10χ,b=>10χ,c=>10χ,i=>χ,j=>χ,k=>χ) backend=eTbackend begin
        rho_vovo[a,i,b,j] += -2*ct_vovo[a,i,b,k]*d_oo[k,j]
        rho_vovo[a,i,b,j] += 2*ct_vovo[a,i,c,j]*d_vv[b,c]
        rho_vovo[a,i,b,j] += cs_vo[a,i]*d_vv[b,c]*s_vo[c,j]
        rho_vovo[a,i,b,j] += -cs_vo[a,i]*d_oo[k,j]*s_vo[b,k]
        rho_vovo[a,i,b,j] += -cs_vo[a,k]*d_oo[k,i]*s_vo[b,j]
        rho_vovo[a,i,b,j] += cs_vo[c,i]*d_vv[a,c]*s_vo[b,j]
        rho_vovo[a,i,b,j] += -2*cs_vovo[a,i,b,k]*d_oo[k,j]*γ₁
        rho_vovo[a,i,b,j] += 2*cs_vovo[a,i,c,j]*d_vv[b,c]*γ₁
        rho_vovo[a,i,b,j] += cγ*d_vv[a,c]*s_vovo[b,j,c,i]
        rho_vovo[a,i,b,j] += -cγ*d_oo[k,i]*s_vovo[a,k,b,j]
        rho_vovo[a,i,b,j] += -4*cs_vovo[a,i,b,k]*d_ov[k,c]*s_vo[c,j]
        rho_vovo[a,i,b,j] += -4*cs_vovo[a,i,c,j]*d_ov[k,c]*s_vo[b,k]
        rho_vovo[a,i,b,j] += -2*cs_vo[a,k]*d_ov[k,c]*s_vovo[b,j,c,i]
        rho_vovo[a,i,b,j] += -2*cs_vo[c,i]*d_ov[k,c]*s_vovo[a,k,b,j]
        rho_vovo[a,i,b,j] += cs_vo[a,i]*d_ov[k,c]*v_vovo[b,j,c,k]
        rho_vovo[a,i,b,j] += -ct_vo[a,k]*d_ov[k,c]*t_vovo[b,j,c,i]
        rho_vovo[a,i,b,j] += -ct_vo[c,i]*d_ov[k,c]*t_vovo[a,k,b,j]
        rho_vovo[a,i,b,j] += 2*cv_vovo[a,i,c,k]*d_ov[k,c]*s_vo[b,j]
        rho_vovo[a,i,b,j] += -cs_vo[a,k]*d_ov[k,c]*t_vovo[b,j,c,i]*γ₁
        rho_vovo[a,i,b,j] += -cs_vo[c,i]*d_ov[k,c]*t_vovo[a,k,b,j]*γ₁
        rho_vovo[a,i,b,j] += -ct_vo[a,k]*d_ov[k,c]*s_vo[b,j]*s_vo[c,i]
        rho_vovo[a,i,b,j] += -ct_vo[c,i]*d_ov[k,c]*s_vo[a,k]*s_vo[b,j]
        rho_vovo[a,i,b,j] += -2*ct_vovo[a,i,b,k]*d_ov[k,c]*s_vo[c,j]*γ₁
        rho_vovo[a,i,b,j] += -2*ct_vovo[a,i,c,j]*d_ov[k,c]*s_vo[b,k]*γ₁
        rho_vovo[a,i,b,j] += -ct_vo[a,k]*d_ov[k,c]*s_vovo[b,j,c,i]*γ₁
        rho_vovo[a,i,b,j] += -ct_vo[c,i]*d_ov[k,c]*s_vovo[a,k,b,j]*γ₁
        rho_vovo[a,i,b,j] += -cγ*d_ov[k,c]*s_vo[a,k]*t_vovo[b,j,c,i]
        rho_vovo[a,i,b,j] += -cγ*d_ov[k,c]*s_vo[c,i]*t_vovo[a,k,b,j]
    end
    finalize_eT_function("jacobian_qed_ccsd_bilinear_s2", "qed_ccsd")
end

