const eTbackend = TensorOperations.Backend{:eTbackend}

include("eT_utils.jl")

TensorOperations.scalartype(::Pair) = Float64
function TensorOperations.tensorscalar(p::Pair)
    n, s = p
    n
end

n_intermediates::Int = 0

code_body::IOBuffer = IOBuffer()
output_parameters::Vector{Tuple{String,Int,String}} = Tuple{String,Int,String}[]
input_parameters::Vector{Tuple{String,Int,String}} = Tuple{String,Int,String}[]
local_variables::Vector{Tuple{String,Int,String}} = Tuple{String,Int,String}[]

function reset_state()
    global n_intermediates
    n_intermediates = 0
    take!(code_body)
    empty!(output_parameters)
    empty!(input_parameters)
    empty!(local_variables)
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

function finalize_eT_function()
    println(String(take!(code_body)))
end

function TensorOperations.tensoradd_structure(pC, A::Pair, conjA)
    nA, sA = A
    return sA[collect(linearize(pC))]
end

function TensorOperations.tensorfree!(C::Pair)
    nC, sC = C

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
end

function TensorOperations.tensoralloc_add(TC, pC, A::Pair, conjA, istemp=false,
    backend::TensorOperations.Backend...)
    structure = TensorOperations.tensoradd_structure(pC, A, conjA)
    name = get_intermediate_name(structure)

    eT_alloc(name => structure)

    name => structure
end

function TensorOperations.tensoradd!(C, pC,
    A, conjA,
    α, β,
    backend::eTbackend)

    println("\nAdding: C = α * A + β * C")
    @show A C α β pC
    println()

    p = linearize(pC)

    nA, sA = A
    nC, sC = C

    counts = Dict{String,Int}()
    for x in sA
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

    if iszero(β) && !isone(α)
        println(code_body, "      call zero_array($nC, $dimstr)")
    elseif !isone(β) && !isone(α)
        println(code_body, "      call dscal($dimstr, $(make_eT_num(β)), $nA, 1)")
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
            println(code_body, "      call sort_to_$pstr($nA, $nC")
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

function TensorOperations.tensorcontract_type(TC, pC, A::Pair, pA, conjA,
    B::Pair, pB, conjB, backend...)
    Pair{String,NTuple{sum(length.(pC)),Symbol}}
end

function TensorOperations.tensorcontract_structure(pC,
    A::Pair, pA, conjA,
    B::Pair, pB, conjB)
    nA, sA = A
    nB, sB = B

    return let lA = length(pA[1])
        map(n -> n <= lA ? sA[pA[1][n]] : sB[pB[2][n-lA]], linearize(pC))
    end
end

function TensorOperations.tensoralloc_contract(TC, pC, A::Pair, pA, conjA, B::Pair, pB, conjB,
    istemp=false, backend::TensorOperations.Backend...)
    structure = TensorOperations.tensorcontract_structure(pC, A, pA, conjA, B, pB, conjB)
    name = get_intermediate_name(structure)

    eT_alloc(name => structure)

    name => structure
end

function TensorOperations.tensorcontract!(C, pC,
    A, pA, conjA,
    B, pB, conjB,
    α, β,
    backend::eTbackend)

    println("\nContracting C = α * contract(A, B) + β * C")
    @show A B C α β pA pB pC
    println()

    C
end
