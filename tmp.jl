include("src/omeinsum_impl.jl")

let
    func = FortranFunction(("X", ["v", "v"]))
    A_vovo = ("A_vovo", true)
    B_oo = ("B_oo", true)
    update_code!(func, ein"aibj,ij->ab", 1, [A_vovo, B_oo])
    finalize_eT_function(func, "test_routine", "ccsd")
end
