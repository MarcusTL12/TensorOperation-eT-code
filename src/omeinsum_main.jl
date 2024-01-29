using OMEinsum

include("omeinsum_impl.jl")

function test_code(code, names=make_tmp_inputnames(code),
    outname="X",
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

    func = FortranFunction(IOBuffer(), ("", String[]),
        Tuple{String,Vector{String}}[], Tens[], Ref(0), Ref(false))

    dimdict = make_ov_dimdict(code)

    prefactor = 1

    make_code!(func, choices, names, perms, outname, outperm, dimdict, prefactor, steps)
    println()

    println(String(take!(func.code_body)))
end

function test1()
    test_code(ein"ik,jk->ij")
end

function test2()
    code = ein"ak,kcld,dj,blci->aibj"

    inputperms = make_trivial_inputperms(code)

    push!(inputperms[2], [3, 4, 1, 2])
    push!(inputperms[4], [3, 4, 1, 2])

    outperms = make_trivial_outperm(code)
    push!(outperms, [3, 4, 1, 2])

    test_code(code, [("A", true), ("B", true), ("C", true), ("D", true)], "X",
        inputperms, outperms)
end

function test3()
    test_code(ein"ai,aijj->")
end

function test4()
    test_code(EinCode(((1, 3), (2, 3)), (1, 2)))
end

function test5()
    test_code(ein"ai,bj,jb->ai")
end

function test6()
    test_code(ein"ci,ackd,bjdk->aibj")
end

function test7()
    test_code(ein",ak,,kibj->aibj")
end

function test7_5()
    test_code(ein"cl,,ak,,kibj,cl->aibj")
end

function test8()
    test_code(ein"ak,bj,ik,jb->ai")
end

function test9()
    test_code(ein"ai,jb->aibj")
end

function test10()
    test_code(ein"aibj,bi->ja")
end

function test11()
    test_code(ein"(((abcd,ia),jb),kc),ld->ijkl")
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
    test_code(ein"ak,ik,bj,jb,cl,cl->ai")
end

function test14()
    test_code(ein"ak,ik,bj,jb,->ai")
end

function test15()
    test_code(ein"ai,->ai")
end

function test16()
    test_code(ein"ia,bj,bj->ai")
end

function test17()
    test_code(ein",bj,bj->")
end

function test18()
    test_code(ein",iijj->")
end

function test18()
    test_code(ein",ijji->")
end
