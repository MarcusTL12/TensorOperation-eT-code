using TupleTools
using SymPy

const eTbackend = TensorOperations.Backend{:eTbackend}

include("eT_utils.jl")

n_intermediates::Int = 0

code_body::IOBuffer = IOBuffer()
output_parameters::Vector{Pair} = []
input_parameters::Vector{Pair} = []
local_variables::Vector{Tuple{String,Int}} = Tuple{String,Int}[]
n_integers::Int = 0
use_ddot::Bool = false

function reset_state()
    global n_intermediates
    n_intermediates = 0
    global use_ddot
    use_ddot = false
    global n_integers
    n_integers = 0
    take!(code_body)
    empty!(output_parameters)
    empty!(input_parameters)
    empty!(local_variables)
end

function input_scalar(s)
    push!(input_parameters, s => () => ())
    Sym(s)
end

function input_tensor(name, structure)
    tens = name => structure => (1:length(structure)...,)
    push!(input_parameters, tens)
    tens
end

function output_scalar(s)
    push!(output_parameters, s => () => ())
    Sym(s)
end

function output_tensor(name, structure)
    tens = name => structure => (1:length(structure)...,)
    push!(output_parameters, tens)
    tens
end

function noio_tensor(name, structure)
    name => structure => (1:length(structure)...,)
end

TensorOperations.scalartype(::Pair) = Float64
function TensorOperations.tensorscalar(p::Pair)
    n, (s, _) = p
    Sym(n)
end

function TensorOperations.tensorscalar(p::Sym)
    p
end

function Base.:*(p::Sym, ::TensorOperations.VectorInterface.One)
    p
end

function Base.:*(::TensorOperations.VectorInterface.One, p::Sym)
    p
end

function Base.:*(::Sym, ::TensorOperations.VectorInterface.Zero)
    Sym(0)
end

function Base.:*(::TensorOperations.VectorInterface.Zero, ::Sym)
    Sym(0)
end

function get_intermediate_name(structure)
    global n_intermediates
    n_intermediates += 1
    if isempty(structure)
        "X$(n_intermediates)"
    else
        "X$(n_intermediates)_" * prod(String, structure)
    end
end

function finalize_eT_function(routine_name, wf_type="ccs")
    io = IOBuffer()

    routine_name = "$(routine_name)_$(wf_type)"

    sort!(output_parameters)
    sort!(input_parameters)

    print(io, "   subroutine $routine_name(wf")

    for param in output_parameters
        param_str = if param isa Sym
            string(param)
        else
            param[1]
        end
        print(io, ", $param_str")
    end

    for param in input_parameters
        param_str = if param isa Sym
            string(param)
        else
            param[1]
        end
        print(io, ", $param_str")
    end

    println(io, ")")

    println(
        io,
        """
!!
!! Generated function
!!
      implicit none
!"""
    )

    println(io, "      class($wf_type), intent(in) :: wf\n!")

    for (params, intent) in ((output_parameters, "inout"), (input_parameters, "in"))
        structure_dict = Dict()

        for (name, (structure, _)) in params
            if !haskey(structure_dict, structure)
                structure_dict[structure] = String[]
            end
            push!(structure_dict[structure], name)
        end

        structure_sort = [(length(s), s, names) for (s, names) in structure_dict]
        sort!(structure_sort)

        for (_, structure, names) in structure_sort
            print(io, "      real(dp)")
            if length(structure) > 0
                print(io, ", dimension(")
                isfirst = true
                for d in structure
                    if isfirst
                        isfirst = false
                    else
                        print(io, ",")
                    end

                    print(io, eT_dim_dict[d])
                end
                print(io, ")")
            end

            print(io, ", intent($intent) :: ")

            isfirst = true
            for name in names
                if isfirst
                    isfirst = false
                else
                    print(io, ", ")
                end
                print(io, name)
            end
            println(io, "")
        end

        println(io, "!")
    end

    dim_dict = Dict()

    for (name, n) in local_variables
        if !haskey(dim_dict, n)
            dim_dict[n] = String[]
        end
        push!(dim_dict[n], name)
    end

    dim_sort = sort!(collect(dim_dict))

    for (n, names) in dim_sort
        print(io, "      real(dp)")
        if n > 0
            print(io, ", dimension(")
            isfirst = true
            for _ in 1:n
                if isfirst
                    isfirst = false
                else
                    print(io, ",")
                end

                print(io, ":")
            end
            print(io, "), allocatable")
        end

        print(io, " :: ")

        isfirst = true
        for name in names
            if isfirst
                isfirst = false
            else
                print(io, ", ")
            end
            print(io, name)
        end
        println(io, "\n!")
    end

    if n_integers > 0
        print(io, "      integer :: ")
        isfirst = true
        for i in 1:n_integers
            if isfirst
                isfirst = false
            else
                print(io, ", ")
            end
            print(io, "i$i")
        end
    end

    println(io, "")

    if use_ddot
        println(io, "      real(dp), external :: ddot")
    end

    println(io, "!")

    print(io, String(take!(code_body)))

    println(io, "!\n   end subroutine $routine_name")

    String(take!(io))
end

function TensorOperations.tensoradd_structure(pC, A::Pair, conjA)
    nA, (sA, _) = A
    return sA[collect(linearize(pC))]
end

function TensorOperations.tensorfree!(C::Pair)
    nC, (sC, _) = C

    println("Deallocating tensor $nC")

    if !isempty(sC)
        println(code_body, "      call mem%dealloc($nC)")
    end
end

function eT_alloc(A::Pair)
    nA, sA = A

    println("Allocating intermediate $nA with size $sA")

    if !isempty(sA)
        print(code_body, "      call mem%alloc($nA")
        for d in sA
            print(code_body, ", ", eT_dim_dict[d])
        end
        println(code_body, ")")
    end

    if (nA, length(sA)) ∉ local_variables
        push!(local_variables, (nA, length(sA)))
    end
end

function TensorOperations.tensoralloc_add(TC, pC, A::Pair, conjA, istemp=false,
    backend::TensorOperations.Backend...)
    if isempty(linearize(pC))
        name = get_intermediate_name(())

        eT_alloc(name => ())

        Sym(name)
    else
        structure = TensorOperations.tensoradd_structure(pC, A, conjA)

        name = get_intermediate_name(structure)

        eT_alloc(name => structure)

        name => structure => (1:length(structure)...,)
    end
end

function get_dimstr(structure)
    counts = Dict{String,Int}()
    for x in structure
        counts[x] = get(counts, x, 0) + 1
    end

    dimbuf = IOBuffer()

    isfirst = true
    for (d, c) in counts
        if isfirst
            isfirst = false
        else
            print(dimbuf, "*")
        end
        print(dimbuf, eT_dim_dict[d])
        if c > 1
            print(dimbuf, "**", c)
        end
    end

    dimstr = String(take!(dimbuf))
end

function TensorOperations.tensoradd!(C, pC,
    A, conjA,
    α, β,
    backend::eTbackend)

    println("\nAdding: C = α * A + β * C")
    @show A C α β pC
    println()

    p = linearize(pC)

    nA, (sA, lpA) = A
    nC, (sC, lpC) = C

    # @assert issorted(lpA)
    # @assert issorted(lpC)

    A_loc = (nA, length(sA))
    C_loc = (nC, length(sC))

    dimstr = get_dimstr(sA)

    if iszero(β) && !isone(α)
        println(code_body, "      call zero_array($nC, $dimstr)")
    elseif !isone(β) && !isone(α)
        println(code_body, "      call dscal($dimstr, $(make_eT_num(β)), $nC, 1)")
    end

    if issorted(p)
        if isone(α) && iszero(β)
            println(code_body, "      call dcopy($dimstr, $nA, 1, $nC, 1)")
        else
            println(code_body, "      call daxpy($dimstr, $(make_eT_num(α)), $nA, 1, $nC, 1)")
        end
    else
        if isone(α) && iszero(β)
            pstr = prod(string, p)
            print(code_body, "      call sort_to_$pstr($nA, $nC")
            for d in sA
                print(code_body, ", ", eT_dim_dict[d])
            end
            println(code_body, ")")
        else
            pstr = prod(string, invperm(p))
            pstr_sort = prod(string, 1:length(p))
            print(code_body, "      call add_$(pstr)_to_$pstr_sort($(make_eT_num(α)), $nA, $nC")
            for d in sC
                print(code_body, ", ", eT_dim_dict[d])
            end
            println(code_body, ")")
        end
    end

    C
end

function TensorOperations.tensortrace!(C, pC,
    A, pA, conjA,
    α, β,
    backend::eTbackend)

    println("\nTracing C = α * tr(A) + β * C")
    @show A C α β pA pC
    println()

    n_loop_indices = length(linearize(pC)) + length(pA[1])

    global n_integers
    n_integers = max(n_integers, n_loop_indices)

    nA, (sA, lpA) = A

    nC, (sC, lpC) = if C isa Sym
        string(C), ((), ())
    else
        C
    end

    @assert issorted(lpA)
    @assert issorted(lpC)

    A_loc = (nA, length(sA))
    C_loc = (nC, length(sC))

    used_inds = 0
    input_inds = zeros(Int, length(sA))

    for i in linearize(pC)
        used_inds += 1
        input_inds[i] = used_inds
    end

    for (i, j) in zip(pA...)
        used_inds += 1
        input_inds[i] = input_inds[j] = used_inds
    end

    input_ind_dims = ["" for _ in 1:n_loop_indices]
    for (i, j) in enumerate(input_inds)
        input_ind_dims[j] = eT_dim_dict[sA[i]]
    end

    tab_level = 2

    dimstr = get_dimstr(TupleTools.getindices(sA, linearize(pC)))

    if !iszero(β) && !isone(β)
        if isempty(linearize(pC))
            print(code_body, "      $nC = $nC * ", make_eT_num(β))
        else
            print(code_body, "      call dscal($dimstr, $(make_eT_num(β)), $nC, 1)")
        end
    elseif iszero(β)
        if isempty(linearize(pC))
            println(code_body, "      $nC = zero")
        else
            println(code_body, "      call zero_array($nC, $dimstr)")
        end
    end

    println(code_body, "!")

    print(code_body, "!\$omp parallel do schedule(static)")
    if n_loop_indices > 1
        print(code_body, " collapse($n_loop_indices)")
    end
    print(code_body, " private(")

    isfirst = true
    for i in 1:n_loop_indices
        if isfirst
            isfirst = false
        else
            print(code_body, ",")
        end
        print(code_body, "i$i")
    end
    println(code_body, ")")

    for i in reverse(1:n_loop_indices)
        println(code_body, "   "^tab_level, "do i$i = 1, ", input_ind_dims[i])
        tab_level += 1
    end

    out_string = IOBuffer()

    print(out_string, "$nC")
    if !isempty(linearize(pC))
        print(out_string, "(")
        isfirst = true
        for i in linearize(pC)
            if isfirst
                isfirst = false
            else
                print(out_string, ",")
            end
            print(out_string, "i$(input_inds[i])")
        end
        print(out_string, ")")
    end

    out_string = String(take!(out_string))

    print(code_body, "   "^tab_level, out_string, " = ")

    print(code_body, out_string, " + ")

    if !isone(α)
        if isone(-α)
            print(code_body, "-")
        else
            print(code_body, make_eT_num(α), "*")
        end
    end

    print(code_body, "$nA(")

    isfirst = true
    for i in input_inds
        if isfirst
            isfirst = false
        else
            print(code_body, ",")
        end
        print(code_body, "i$i")
    end
    println(code_body, ")")

    for i in 1:n_loop_indices
        tab_level -= 1
        println(code_body, "   "^tab_level, "end do")
    end

    println(code_body, "!\$omp end parallel do")

    println(code_body, "!")

    C
end

function compose_perms(pA, lpA)
    lin = TupleTools.getindices(linearize(pA), lpA)

    r1 = 1:length(pA[1])
    r2 = length(pA[1])+1:length(lin)

    lin[r1], lin[r2]
end

function compose_iperms(pA, lpA)
    lin = TupleTools.getindices(linearize(pA), invperm(lpA))

    r1 = 1:length(pA[1])
    r2 = length(pA[1])+1:length(lin)

    lin[r1], lin[r2]
end

function TensorOperations.tensorcontract_structure(pC,
    A::Union{Pair,Sym}, pA, conjA,
    B::Union{Pair,Sym}, pB, conjB)
    nA, (sA, lpA) = if A isa Sym
        string(A), ((), ())
    else
        A
    end

    nB, (sB, lpB) = if B isa Sym
        string(B), ((), ())
    else
        B
    end

    pA = compose_perms(pA, lpA)
    pB = compose_perms(pB, lpB)

    lA = length(pA[1])
    sC = map(n -> n <= lA ? sA[pA[1][n]] : sB[pB[2][n-lA]], linearize(pC))

    if issorted(linearize(pC))
        sC => linearize(pC), true
    elseif issorted(linearize(reverse(pC)))
        map(n -> n <= lA ? sA[pA[1][n]] : sB[pB[2][n-lA]], linearize(reverse(pC))) => linearize(pC), true
    else
        sC => linearize(pC), false
    end
end

function TensorOperations.tensoralloc_contract(TC, pC, A::Pair, pA, conjA, B::Pair, pB, conjB,
    istemp=false, backend::TensorOperations.Backend...)
    if isempty(linearize(pC))
        name = get_intermediate_name(())

        eT_alloc(name => ())

        Sym(name)
    else
        (structure, lpC), _ = TensorOperations.tensorcontract_structure(pC, A, pA, conjA, B, pB, conjB)
        name = get_intermediate_name(structure)

        eT_alloc(name => structure)

        name => structure => lpC
    end
end

function sort_tensor_input(A, pA, backend)
    nA, (sA, _) = if A isa Sym
        string(A), ((), ())
    else
        A
    end

    if issorted(linearize(pA))
        println("$A is sorted")
        A, false, false
    elseif issorted(linearize(reverse(pA)))
        println("$A is transposed")
        A, true, false
    else
        do_trans = linearize(get_pT(pA)) < linearize(pA)

        if do_trans
            pA = get_pT(pA)
        end

        println("$A not sorted, need to allocate tmp storage")

        sA_sort = TupleTools.getindices(sA, linearize(pA))

        nA_sort = get_intermediate_name(sA_sort)

        eT_alloc(nA_sort => sA_sort)

        A_sort = nA_sort => sA_sort => (1:length(linearize(pA))...,)

        TensorOperations.tensoradd!(A_sort, pA, nA => (sA, (1:length(linearize(pA))...,)), :N, 1, 0, backend)

        A_sort, do_trans, true
    end
end

function get_pT(pA)
    lin = linearize(reverse(pA))

    r1 = 1:length(pA[1])
    r2 = length(pA[1])+1:length(lin)

    lin[r1], lin[r2]
end

function invT(A, pA)
    nA, (sA, lpA) = A

    r1 = 1:length(pA[1])
    r2 = length(pA[1])+1:length(linearize(pA))

    lpA1 = lpA[r1]
    lpA2 = lpA[r2]

    nA => sA => linearize((lpA2, lpA1))
end

function TensorOperations.tensorcontract!(C, pC,
    A, pA, conjA,
    B, pB, conjB,
    α, β,
    backend::eTbackend)

    println("\nContracting C = α * contract(A, B) + β * C")
    @show A B C α β pA pB pC
    println()

    rpA = reverse(pA)
    rpB = reverse(pB)
    indCinoBA = let N₁ = TensorOperations.numout(pA), N₂ = TensorOperations.numin(pB)
        map(n -> ifelse(n > N₁, n - N₁, n + N₂), linearize(pC))
    end
    tpC = TensorOperations.trivialpermutation(pC)
    rpC = (TupleTools.getindices(indCinoBA, tpC[1]),
        TupleTools.getindices(indCinoBA, tpC[2]))

    if TensorOperations.tensorcontract_structure(pC, A, pA, conjA, B, pB, conjB)[2]
        eT_contract(C, pC, A, pA, B, pB, α, β, backend)
    elseif TensorOperations.tensorcontract_structure(rpC, B, rpB, conjB, A, rpA, conjA)[2]
        println("Swapping multiplication order")
        eT_contract(C, rpC, B, rpB, A, rpA, α, β, backend)
    elseif linearize(rpC) < linearize(pC)
        println("Swapping multiplication order")
        eT_contract(C, rpC, B, rpB, A, rpA, α, β, backend)
    else
        eT_contract(C, pC, A, pA, B, pB, α, β, backend)
    end
end

function eT_contract(C, pC,
    A, pA,
    B, pB,
    α, β, backend::eTbackend)

    nA, (sA, lpA) = if A isa Sym
        string(A), ((), ())
    else
        A
    end

    nB, (sB, lpB) = if B isa Sym
        string(B), ((), ())
    else
        B
    end

    nC, (sC, lpC) = if C isa Sym
        string(C), ((), ())
    else
        C
    end

    A_loc = (nA, length(sA))
    B_loc = (nB, length(sB))
    C_loc = (nC, length(sC))

    pA = compose_iperms(pA, lpA)
    pB = compose_iperms(pB, lpB)
    pC = compose_iperms(pC, lpC)

    @show pA pB pC

    if length(linearize(pA)) > length(linearize(pB)) && !issorted(pA[2])
        println("Swapping permutation order of contraction indices")
        P = TupleTools.sortperm(pA[2])
        pA = (pA[1], TupleTools.getindices(pA[2], P))
        pB = (TupleTools.getindices(pB[1], P), pB[2])
    elseif length(linearize(pA)) < length(linearize(pB)) && !issorted(pB[1])
        println("Swapping permutation order of contraction indices")
        P = TupleTools.sortperm(pB[1])
        pA = (pA[1], TupleTools.getindices(pA[2], P))
        pB = (TupleTools.getindices(pB[1], P), pB[2])
    end

    @show pA pB pC

    A_sort, A_T, A_alloc = sort_tensor_input(A, pA, backend)

    nA_sort, (sA_sort, _) = if A_sort isa Sym
        string(A_sort), ((), ())
    else
        A_sort
    end

    B_sort, B_T, B_alloc = sort_tensor_input(B, pB, backend)

    nB_sort, (sB_sort, _) = if B_sort isa Sym
        string(B_sort), ((), ())
    else
        B_sort
    end

    outdim1 = get_dimstr(TupleTools.getindices(sA, pA[1]))
    outdim2 = get_dimstr(TupleTools.getindices(sB, pB[2]))
    contdim = get_dimstr(TupleTools.getindices(sA, pA[2]))

    lda = get_dimstr(TupleTools.getindices(sA, pA[A_T ? 2 : 1]))
    ldb = get_dimstr(TupleTools.getindices(sB, pB[B_T ? 2 : 1]))
    ldc = outdim1

    A_T_str = A_T ? "'T'" : "'N'"
    B_T_str = B_T ? "'T'" : "'N'"

    # Check if we get away with ddot or add

    return_early = false

    if outdim1 == "" && outdim2 == "" && contdim == ""
        println("Scalar scalar multiplication")

        print(code_body, "      $nC = ")

        if !iszero(β)
            if !isone(β)
                print(code_body, make_eT_num(β), "*")
            end
            print(code_body, "$nC + ")
        end

        if !isone(α)
            print(code_body, make_eT_num(α), "*")
        end

        println(code_body, "$nA_sort*$nB_sort")
        return_early = true
    elseif outdim1 == "" && outdim2 == ""
        global use_ddot
        use_ddot = true
        println("Using ddot")
        print(code_body, "      $nC = ")
        if !iszero(β)
            if !isone(β)
                print(code_body, "$(make_eT_num(β)) * ")
            end
            print(code_body, "$nC + ")
        end

        if !isone(α)
            print(code_body, "$(make_eT_num(α)) * ")
        end
        println(code_body, "ddot($contdim, $nA_sort, 1, $nB_sort, 1)")
        return_early = true
    elseif contdim == "" && outdim2 == ""
        println("Using tensoradd")

        TensorOperations.tensoradd!(C, pC, A_sort, :N, Sym(nB) * α, β, backend)
        return_early = true
    elseif contdim == "" && outdim1 == ""
        println("Using tensoradd")

        TensorOperations.tensoradd!(C, pC, B_sort, :N, Sym(nA) * α, β, backend)
        return_early = true
    end

    if return_early
        if A_alloc
            TensorOperations.tensorfree!(A_sort)
        end

        if B_alloc
            TensorOperations.tensorfree!(B_sort)
        end

        return C
    end

    # Checking if output must be sorted

    ipC = invperm(linearize(pC))

    explicit_sort_output = false

    β_out = β

    C_pre = if issorted(linearize(pC))
        println("Output is sorted")
        C
    else
        println("Output not sorted, need to allocate tmp output!")

        explicit_sort_output = true

        β = 0

        nC_pre = nC * "_" * prod(string, ipC)

        sC_pre = TupleTools.getindices(sC, ipC)

        eT_alloc(nC_pre => sC_pre)

        nC_pre => sC_pre => (1:length(sC_pre)...,)
    end

    nC_out, (sC_out, _) = C_pre

    # Doing the contraction

    if contdim == ""
        println("Using dger")

        if !isone(β)
            dimstr = get_dimstr(sC_out)
            if iszero(β)
                println(code_body, "      call zero_array($nC_out, $dimstr)")
            else
                println(code_body, "      call dscal($dimstr, $(make_eT_num(β)), $nC_out, 1)")
            end
        end

        println(
            code_body,
            """
    !
          call dger($outdim1, &
                    $outdim2, &
                    $(make_eT_num(α)), &
                    $nA_sort, 1, &
                    $nB_sort, 1, &
                    $nC_out, $ldc)
    !""")
    elseif outdim1 == ""
        println("Using left dgemv")

        n = B_T ? contdim : outdim2

        T_str = B_T ? "'N'" : "'T'"

        println(
            code_body,
            """
    !
          call dgemv($T_str, &
                     $ldb, &
                     $n, &
                     $(make_eT_num(α)), &
                     $nB_sort, &
                     $ldb, &
                     $nA_sort, 1, &
                     $(make_eT_num(β)), &
                     $nC_out, 1)
    !""")
    elseif outdim2 == ""
        println("Using right dgemv")

        n = A_T ? outdim1 : contdim

        println(
            code_body,
            """
    !
          call dgemv($A_T_str, &
                     $lda, &
                     $n, &
                     $(make_eT_num(α)), &
                     $nA_sort, &
                     $lda, &
                     $nB_sort, 1, &
                     $(make_eT_num(β)), &
                     $nC_out, 1)
    !""")
    else
        println("Using dgemm")
        println(
            code_body,
            """
    !
          call dgemm($A_T_str, $B_T_str, &
                     $outdim1, &
                     $outdim2, &
                     $contdim, &
                     $(make_eT_num(α)), &
                     $nA_sort, &
                     $lda, &
                     $nB_sort, &
                     $ldb, &
                     $(make_eT_num(β)), &
                     $nC_out, &
                     $ldc)
    !""")
    end

    if A_alloc
        TensorOperations.tensorfree!(A_sort)
    end

    if B_alloc
        TensorOperations.tensorfree!(B_sort)
    end

    if explicit_sort_output
        println("Sorting output")

        TensorOperations.tensoradd!(C, pC, C_pre, :N, 1, β_out, backend)

        TensorOperations.tensorfree!(C_pre)
    end

    C
end
