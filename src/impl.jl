using TupleTools

const eTbackend = TensorOperations.Backend{:eTbackend}

include("eT_utils.jl")

n_intermediates::Int = 0

code_body::IOBuffer = IOBuffer()
output_parameters::Vector{Pair} = Pair[]
input_parameters::Vector{Pair} = Pair[]
local_variables::Vector{Tuple{String,Int}} = Tuple{String,Int}[]
use_ddot::Bool = false

function reset_state()
    global n_intermediates
    n_intermediates = 0
    global use_ddot
    use_ddot = false
    take!(code_body)
    empty!(output_parameters)
    empty!(input_parameters)
    empty!(local_variables)
end

TensorOperations.scalartype(::Pair) = Float64
function TensorOperations.tensorscalar(p::Pair)
    n, (s, _) = p
    n
end

function Base.:*(p::Pair{String,Pair{Tuple{},Tuple{}}}, s::String)
    n, _ = p
    "$n*$s"
end

function Base.:*(a::Int, b::String)
    "$(make_eT_num(a))*$b"
end

function Base.:+(a::String, b::String)
    "$a+$b"
end

function Base.:-(a::String, b::String)
    "$a-$b"
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

    print(io, "   subroutine $routine_name(wf")

    for (param, _) in output_parameters
        print(io, ", $param")
    end

    for (param, _) in input_parameters
        print(io, ", $param")
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

    for (params, intent) in ((output_parameters, "out"), (input_parameters, "in"))
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
        println(io, "")
    end

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
    structure = TensorOperations.tensoradd_structure(pC, A, conjA)
    name = get_intermediate_name(structure)

    eT_alloc(name => structure)

    name => structure => (1:length(structure)...,)
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

    nA, (sA, _) = A
    nC, (sC, _) = C

    A_loc = (nA, length(sA))
    C_loc = (nC, length(sC))

    if A_loc ∉ local_variables && A ∉ input_parameters
        push!(input_parameters, A)
    end

    if C_loc ∉ local_variables && C ∉ output_parameters
        push!(output_parameters, C)
    end

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
    A::Pair, pA, conjA,
    B::Pair, pB, conjB)
    nA, (sA, lpA) = A
    nB, (sB, lpB) = B

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
    (structure, lpC), _ = TensorOperations.tensorcontract_structure(pC, A, pA, conjA, B, pB, conjB)
    name = get_intermediate_name(structure)

    eT_alloc(name => structure)

    name => structure => lpC
end

function sort_tensor_input(A, pA, backend)
    nA, (sA, _) = A

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

        nA_sort = nA * "_" * prod(string, linearize(pA))

        sA_sort = TupleTools.getindices(sA, linearize(pA))

        eT_alloc(nA_sort => sA_sort)

        A_sort = nA_sort => sA_sort => (1:length(linearize(pA))...,)

        TensorOperations.tensoradd!(A_sort, pA, A, :N, 1, 0, backend)

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
        eT_contract(C, rpC, invT(B, pB), rpB, invT(A, pA), rpA, α, β, backend)
    elseif linearize(rpC) < linearize(pC)
        println("Swapping multiplication order")
        eT_contract(C, rpC, invT(B, pB), rpB, invT(A, pA), rpA, α, β, backend)
    else
        eT_contract(C, pC, A, pA, B, pB, α, β, backend)
    end
end

function eT_contract(C, pC,
    A, pA,
    B, pB,
    α, β, backend::eTbackend)

    nA, (sA, lpA) = A
    nB, (sB, lpB) = B
    nC, (sC, lpC) = C

    A_loc = (nA, length(sA))
    B_loc = (nB, length(sB))
    C_loc = (nC, length(sC))

    if A_loc ∉ local_variables && A ∉ input_parameters
        push!(input_parameters, A)
    end

    if B_loc ∉ local_variables && B ∉ input_parameters
        push!(input_parameters, B)
    end

    if C_loc ∉ local_variables && C ∉ output_parameters
        push!(output_parameters, C)
    end

    pA = compose_iperms(pA, lpA)
    pB = compose_iperms(pB, lpB)
    pC = compose_iperms(pC, lpC)

    @show pA pB pC

    A_sort, A_T, A_alloc = sort_tensor_input(A, pA, backend)

    nA_sort, (sA_sort, _) = A_sort

    B_sort, B_T, B_alloc = sort_tensor_input(B, pB, backend)

    nB_sort, (sB_sort, _) = B_sort

    # Checking if output must be sorted

    ipC = invperm(linearize(pC))

    explicit_sort_output = false

    C_pre = if issorted(linearize(pC))
        println("Output is sorted")
        C
    else
        println("Output not sorted, need to allocate tmp output!")

        explicit_sort_output = true

        nC_pre = nC * "_" * prod(string, ipC)

        sC_pre = TupleTools.getindices(sC, ipC)

        eT_alloc(nC_pre => sC_pre)

        nC_pre => sC_pre => :N
    end

    nC_out, (sC_out, _) = C_pre

    # Doing the contraction

    outdim1 = get_dimstr(TupleTools.getindices(sA, pA[1]))
    outdim2 = get_dimstr(TupleTools.getindices(sB, pB[2]))
    contdim = get_dimstr(TupleTools.getindices(sA, pA[2]))

    lda = get_dimstr(TupleTools.getindices(sA_sort, pA[A_T ? 2 : 1]))
    ldb = get_dimstr(TupleTools.getindices(sB_sort, pB[B_T ? 2 : 1]))
    ldc = outdim1

    A_T_str = A_T ? "'T'" : "'N'"
    B_T_str = B_T ? "'T'" : "'N'"

    # Bodge to use dgemm for everything
    # TODO: recognize when to use dgemv, ddot and dger

    if outdim1 == "" && outdim2 == ""
        global use_ddot
        use_ddot = true
        println("Using ddot")
        print(code_body, "      $nC_out = ")
        if !iszero(β)
            if !isone(β)
                print(code_body, "$(make_eT_num(β)) * ")
            end
            print(code_body, "$nC_out + ")
        end

        if !isone(α)
            print(code_body, "$(make_eT_num(α)) * ")
        end
        println(code_body, "ddot($contdim, $nA_sort, 1, $nB_sort, 1)")
    elseif contdim == ""
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
                    $nC_out, 1)
    !""")
    elseif outdim1 == ""
        println("Using left dgemv")

        m = B_T ? outdim2 : contdim
        n = B_T ? contdim : outdim2

        T_str = B_T ? "'N'" : "'T'"

        println(
            code_body,
            """
    !
          call dgemv($T_str, &
                     $m, &
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

        m = A_T ? outdim1 : contdim
        n = A_T ? contdim : outdim1

        println(
            code_body,
            """
    !
          call dgemv($A_T_str, &
                     $m, &
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

        TensorOperations.tensoradd!(C, pC, C_pre, :N, 1, 0, backend)

        TensorOperations.tensorfree!(C_pre)
    end

    C
end
