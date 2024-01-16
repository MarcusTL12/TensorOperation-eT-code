using Strided

const eTbackend = TensorOperations.Backend{:eTbackend}

TensorOperations.scalartype(::Pair) = Float64

n_intermediates::Int = 0

function get_intermediate_name(structure)
    global n_intermediates
    n_intermediates += 1
    "X$(n_intermediates)_" * prod(String, structure)
end

function reset_state()
    global n_intermediates
    n_intermediates = 0
end

function TensorOperations.tensoradd_structure(pC, A::Pair, conjA)
    nA, sA = A
    return sA[collect(linearize(pC))]
end

function TensorOperations.tensoralloc_add(TC, pC, A::Pair, conjA, istemp=false,
    backend::TensorOperations.Backend...)
    structure = TensorOperations.tensoradd_structure(pC, A, conjA)
    name = get_intermediate_name(structure)

    println("Allocating intermediate $name with size $structure")

    name => structure
end

function TensorOperations.tensoradd!(C, pC,
    A, conjA,
    α, β,
    backend::eTbackend)

    println("\nAdding: C = α * A + β * C")
    @show A C α β pC
    println()

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

    println("Allocating intermediate $name with size $structure")

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
