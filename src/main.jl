using TensorOperations

include("impl.jl")

function test()
    A = randn(10, 10)
    B = randn(10, 10)

    @tensor backend=eTbackend C[i, j] := A[i, k] * B[k, j]

    nothing
end

function test_sym()
    A = "A" => (:o, :v)
    B = "B" => (:v, :o)
    C = "C" => (:o, :o)

    reset_state()

    @tensor backend=eTbackend C[j, i] += A[i, k] * B[k, j]
end

function test2()
    g = "g" => (:v, :v, :o, :o)
    t1 = "t1" => (:v, :o)
    h_ov = "h" => (:o, :v)

    reset_state()

    @tensor backend=eTbackend C[a, i] := g[a, b, j, j] * t1[b, i]
end

function test3()
    A = "A" => (:v, :o, :v, :o)
    B = "B" => (:o, :v, :o, :v)
    C = "C" => (:o, :o, :o, :o)

    reset_state()

    @tensor backend=eTbackend C[i, j, k, l] = A[a, j, b, i] * B[k, b, l, a]
end

function test4()
    A = "A" => ("v", "o", "v", "o")
    B = "B" => ("v", "v", "o", "o")
    C = "C" => ("v", "o", "v", "o")

    reset_state()

    @tensor backend=eTbackend C[a, i, b, j] = A[a, i, b, j] - 2 * B[a, b, j, i]
end

function test5()
    A = "A" => (:o, :v)
    B = "B" => (:v, :o)
    C = "C" => (:o, :v)
    D = "D" => (:o, :v)

    reset_state()

    @tensor backend=eTbackend D[i, a] = A[i, b] * B[b, j] * C[j, a]
end
