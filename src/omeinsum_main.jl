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
    steps = make_scaled_steps(steps)

    display(steps)

    choices, perms, outperm, sortcost = optimize_choices(cost, steps, inputperms, outperms)

    @show sortcost perms outperm

    println()

    func = FortranFunction(IOBuffer(), outname,
        Tuple{String,Vector{String}}[], Tens[], Ref(0), Ref(false))

    dimdict = make_ov_dimdict(code)

    prefactor = 1

    make_code!(func, choices, names, perms, outname[1], outperm, dimdict, prefactor, steps)
    println()

    # println(String(take!(func.code_body)))
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

    inputperms = make_trivial_inputperms(code)

    push!(inputperms[1], [3, 4, 1, 2])

    outperms = make_trivial_outperm(code)
    push!(outperms, [3, 4, 1, 2])

    test_code(code, [("A", true), ("B", true)], "X", inputperms, outperms)
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
