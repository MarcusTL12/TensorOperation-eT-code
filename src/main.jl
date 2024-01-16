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
