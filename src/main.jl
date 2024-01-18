using TensorOperations

include("impl.jl")

function test()
    A = randn(10, 10)
    B = randn(10, 10)

    @tensor backend=eTbackend C[i, j] := A[i, k] * B[k, j]

    nothing
end

function test_sym()
    A = "A" => ("o", "v") => :N
    B = "B" => ("v", "o") => :N
    C = "C" => ("o", "o") => :N

    reset_state()

    @tensor backend=eTbackend C[j, i] += A[i, k] * B[k, j]
end

function test2()
    g = "g" => ("v", "v", "o", "o") => :N
    t1 = "t1" => ("v", "o") => :N
    h_ov = "h" => ("o", "v") => :N

    C = "C" => ("v", "o")

    reset_state()

    @tensor backend=eTbackend C[a, i] = g[a, b, j, j] * t1[b, i]
end

function test3()
    A = "A" => ("v", "o", "v", "o") => :N
    B = "B" => ("o", "v", "o", "v") => :N
    C = "C" => ("o", "o", "o", "o") => :N

    reset_state()

    @tensor backend=eTbackend C[i, j, k, l] = A[a, j, b, i] * B[k, b, l, a]
end

function test4()
    A = "A" => ("v", "o", "v", "o") => :N
    B = "B" => ("v", "v", "o", "o") => :N
    C = "C" => ("v", "o", "v", "o") => :N

    reset_state()

    @tensor backend=eTbackend C[a, i, b, j] = A[a, i, b, j] - 2 * B[a, b, j, i]
end

function test5()
    A = "A" => ("o", "v") => :N
    B = "B" => ("v", "o") => :N
    C = "C" => ("o", "v") => :N
    D = "D" => ("o", "v") => :N

    reset_state()

    @tensor backend=eTbackend D[i, a] = A[i, b] * B[b, j] * C[j, a]
end

function test6()
    A = "A" => ("o", "v") => :N
    B = "B" => ("v", "o") => :N
    C = "C" => ("o", "v") => :N
    D = "D" => ("o", "v") => :N

    reset_state()

    @tensor opt=(a, b) backend=eTbackend D[i, a] = A[i, a] * B[b, j] * C[j, b]
end

function test7()
    A = "A" => ("v", "o") => :N
    B = "B" => ("v", "o") => :N
    C = "C" => ("o", "v") => :N
    D = "D" => ("o", "v") => :N

    γ = "wf%s0" => () => :N

    reset_state()

    @tensor opt=(a, b) backend=eTbackend D[i, a] = γ * A[a, i] * B[b, j] * C[j, b]
end

function test8()
    A = "A" => ("v", "o", "v", "o") => :N
    B = "B" => ("v", "o", "v", "v") => :N
    C = "C" => ("o", "o", "v", "v") => :N
    D = "D" => ("o", "v") => :N

    reset_state()

    @tensor opt=(a, b, c, d, i, j, k) backend=eTbackend D[i, a] = A[a, j, b, k] * B[b, k, c, d] * C[i, j, c, d]
end

function test9()
    A = "A" => ("v", "o", "v", "o") => :N
    B = "B" => ("o", "v", "v", "v") => :N
    C = "C" => ("o", "o", "v", "v") => :N
    D = "D" => ("o", "v") => :N

    reset_state()

    @tensor opt=(a, b, c, d, i, j, k) backend=eTbackend D[i, a] = A[a, j, b, k] * B[k, b, c, d] * C[i, j, c, d]
end
