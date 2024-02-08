using OMEinsum

include("omeinsum_impl.jl")

function test_code(code, outname, names=make_tmp_inputnames(code),
    inputperms=make_trivial_inputperms(code),
    outperms=make_trivial_outperm(code))
    cost = make_ov_cost(code)

    # optcode = optimize_code(code, cost,
    #     TreeSA(sc_target=28, βs=0.1:0.1:10, ntrials=2, niters=100, sc_weight=3.0))

    optcode = optimize_code(code, cost, GreedyMethod())

    @show optcode

    steps = walk_einsums(optcode)

    display(steps)

    choices, perms, outperm, sortcost = optimize_choices(cost, steps, inputperms, outperms)

    @show sortcost perms outperm

    println()

    func = FortranFunction(outname)

    dimdict = make_ov_dimdict(code)

    prefactor = 1

    make_code!(func, choices, names, perms, outname[1], outperm, dimdict, prefactor, steps)
    println()

    println(finalize_eT_function(func, "test", "ccs"))
end

function test1()
    test_code(ein"ik,jk->ij", ("X", ["o", "o"]))
end

function test2()
    code = ein"ak,kcld,dj,blci->aibj"

    inputperms = make_trivial_inputperms(code)

    push!(inputperms[2], [3, 4, 1, 2])
    push!(inputperms[4], [3, 4, 1, 2])

    outperms = make_trivial_outperm(code)
    push!(outperms, [3, 4, 1, 2])

    test_code(code, ("X", ["v", "o", "v", "o"]),
        [("A", true), ("B", true), ("C", true), ("D", true)],
        inputperms, outperms)
end

function test3()
    test_code(ein"ai,aijj->", ("X", String[]))
end

function test4()
    test_code(EinCode(((1, 3), (2, 3)), (1, 2)), ("X", ["g", "g"]))
end

function test5()
    test_code(ein"ai,bj,jb->ai", ("X", ["v", "o"]))
end

function test6()
    test_code(ein"ci,ackd,bjdk->aibj", ("X", ["v", "o", "v", "o"]))
end

function test7()
    test_code(ein",ak,,kibj->aibj", ("X", ["v", "o", "v", "o"]))
end

function test7_5()
    test_code(ein"cl,,ak,,kibj,cl->aibj", ("X", ["v", "o", "v", "o"]))
end

function test8()
    test_code(ein"ak,bj,ik,jb->ai", ("X", ["v", "o"]))
end

function test9()
    test_code(ein"ai,jb->aibj", ("X", ["v", "o", "v", "o"]))
end

function test10()
    test_code(ein"aibj,bi->ja", ("X", ["o", "v"]))
end

function test11()
    test_code(ein"(((abcd,ia),jb),kc),ld->ijkl", ("X", ["o", "o", "o", "o"]))
end

function test12()
    code = ein"aicj,bc->aibj"

    A = ("A", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    B = ("B", true)

    outperms = [[1, 2, 3, 4], [3, 4, 1, 2]]

    func = FortranFunction(("X", ["v", "o", "v", "o"]))

    update_code!(func, code, -3, [A, B], outperms)

    open("tmp.f90", "w") do io
        println(io, finalize_eT_function(func, "rho_test", "qed_ccsd"))
    end
end

function test13()
    test_code(ein"ak,ik,bj,jb,cl,cl->ai", ("X", ["v", "o"]))
end

function test14()
    test_code(ein"ak,ik,bj,jb,->ai", ("X", ["v", "o"]))
end

function test15()
    test_code(ein"ai,->ai", ("X", ["v", "o"]))
end

function test16()
    test_code(ein"ia,bj,bj->ai", ("X", ["v", "o"]))
end

function test17()
    test_code(ein",bj,bj->", ("X", String[]))
end

function test18()
    test_code(ein",iijj->", ("X", String[]))
end

function test18()
    test_code(ein",ijji->", ("X", String[]))
end

function test19()
    h_oo = ("h_oo", true)
    g_oooo = ("g_oooo", true)

    func = FortranFunction("E")

    update_code!(func, ein"ii->", 2, [h_oo])
    update_code!(func, ein"iijj->", 2, [g_oooo])
    update_code!(func, ein"ijji->", -1, [g_oooo])

    println(finalize_eT_function(func, "hf_energy", "ccs"))
end

function test20()
    #  cs_vo[a,i]*d_ov[k,c]*u_vovo[b,j,c,k]
    # -ct_vo[a,k]*d_ov[k,c]*t_vovo[b,j,c,i]*γ₁
    # -ct_vo[c,i]*d_ov[k,c]*t_vovo[a,k,b,j]*γ₁

    ct = ("ct_vo", true)
    cs = ("cs_vo", true)
    d_ov = ("d_ov", true)
    t_vovo = ("t_vovo", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    u_vovo = ("wf%u_aibj", false, [[1, 2, 3, 4], [3, 4, 1, 2]])
    γ = ("wf%s0", false)

    outperms = [[1, 2, 3, 4], [3, 4, 1, 2]]

    func = FortranFunction(("rho_vovo", ["v", "o", "v", "o"]))

    update_code!(func, ein"ai,kc,bjck->aibj", 1, [cs, d_ov, u_vovo], outperms)
    update_code!(func, ein"ak,kc,bjci,->aibj", -1, [ct, d_ov, t_vovo, γ], outperms)
    update_code!(func, ein"ci,kc,akbj,->aibj", -1, [ct, d_ov, t_vovo, γ], outperms)

    open("tmp.f90", "w") do io
        println(io, finalize_eT_function(func, "rho_test", "qed_ccsd"))
    end
end

function test21()
    # rho[] += 2*F_ov[i,a]*cs_vo[a,i]
    # rho[] += 2*L_ovov[i,a,j,b]*cs_vovo[a,i,b,j]
    # rho[] += 2*L_ovov[i,a,j,b]*ct_vo[a,i]*s_vo[b,j]

    ct = ("ct_vo", true)
    cs = ("cs_vo", true)
    F_ov = ("wf%fock_ia", false)
    s_vo = ("wf%s1", false)
    L_ovov = ("L_ovov", true, [[1, 2, 3, 4], [3, 4, 1, 2]])
    cs_vovo = ("cs_vovo", true, [[1, 2, 3, 4], [3, 4, 1, 2]])

    func = FortranFunction("rho")

    update_code!(func, ein"ia,ai->", 2, [F_ov, cs])
    update_code!(func, ein"iajb,aibj->", 2, [L_ovov, cs_vovo])
    update_code!(func, ein"iajb,ai,bj->", 2, [L_ovov, ct, s_vo])

    open("tmp.f90", "w") do io
        println(io, finalize_eT_function(func, "rho_test", "qed_ccsd"))
    end
end

function test22()
    test_code(ein"(aibj,jb,ac),ci->", ("X", ["v", "o"]))
end

function test23()
    A = ("A", true)
    B = ("B", true)

    func = FortranFunction("X")

    update_code!(func, ein"jb,aibj,ia->", 2, [B, A, B])

    open("tmp.f90", "w") do io
        println(io, finalize_eT_function(func, "rho_test", "qed_ccsd"))
    end
end

function test24()
    test_code(ein",iajb->aibj", ("X", ["v", "o"]))
end

function test25()
    test_code(ein"kaib,k->aib", ("X", ["v", "o", "v"]))
end

function test26()
    test_code(ein"jbkc,ai,bj,ck,->ai", ("X", ["v", "o"]))
end
