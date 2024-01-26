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
    A = "A" => ("o", "v") => (1, 2)
    B = "B" => ("v", "o") => (1, 2)
    C = "C" => ("o", "v") => (1, 2)
    D = "D" => ("o", "v") => (1, 2)

    reset_state()

    @tensor opt=(a, b) backend=eTbackend D[i, a] = A[i, a] * B[b, j] * C[j, b]
end

function test6_5()
    A = "A" => ("o", "v") => (1, 2)
    B = "B" => ("v", "o") => (1, 2)
    C = "C" => ("v", "o") => (1, 2)
    D = "D" => ("o", "v") => (1, 2)

    reset_state()

    @tensor opt=(a, b) backend=eTbackend D[i, a] = A[i, a] * B[b, j] * C[b, j]
end

function test7()
    A = "A" => ("v", "o") => (1, 2)
    B = "B" => ("v", "o") => (1, 2)
    C = "C" => ("o", "v") => (1, 2)
    D = "D" => ("o", "v") => (1, 2)

    γ = "wf%s0" => () => ()

    reset_state()

    @tensor opt=(a, b) backend=eTbackend D[i, a] = γ * A[a, i] * B[b, j] * C[j, b]
end

function test7_5()
    A = "A" => ("v", "o") => (1, 2)
    B = "B" => ("v", "o") => (1, 2)
    C = "C" => ("o", "v") => (1, 2)
    B2 = "B2" => ("v", "o") => (1, 2)
    C2 = "C2" => ("v", "o") => (1, 2)
    D = "D" => ("o", "v") => (1, 2)

    γ = "wf%s0" => () => ()

    reset_state()

    @tensor opt=(a, b, i, j) backend=eTbackend D[i, a] = γ * A[a, i] * (2 * B[b, j] * C[j, b] + B2[b, j] * C2[b, j])
end

function test8()
    A = "A" => ("v", "o", "v", "o") => (1:4...,)
    B = "B" => ("v", "o", "v", "v") => (1:4...,)
    C = "C" => ("o", "o", "v", "v") => (1:4...,)
    D = "D" => ("o", "v") => (1, 2)

    reset_state()

    @tensor opt=(a, b, c, d, i, j, k) backend=eTbackend D[i, a] = A[a, j, b, k] * B[b, k, c, d] * C[i, j, c, d]
end

function test9()
    A = "A" => ("v", "o", "v", "o") => (1:4...,)
    B = "B" => ("o", "v", "v", "v") => (1:4...,)
    C = "C" => ("o", "o", "v", "v") => (1:4...,)
    D = "D" => ("o", "v") => (1, 2)

    reset_state()

    @tensor opt=(a, b, c, d, i, j, k) backend=eTbackend D[i, a] = A[a, j, b, k] * B[k, b, c, d] * C[i, j, c, d]
end

function test10()
    A = "A" => ("v", "v", "o") => (1:3...,)
    B = "B" => ("v", "v", "o", "o") => (1:4...,)
    C = "C" => ("o", "o", "o") => (1, 2, 3)

    reset_state()

    @tensor opt=(a, b, i, j, k) backend=eTbackend C[i, j, k] = B[a, b, k, j] * A[a, b, i]
end

function test11()
    A = "A" => ("v", "o") => (1, 2)
    B = "B" => ("v", "o") => (1, 2)
    C = "C" => ("v", "o", "o", "v") => (1, 2, 3, 4)

    reset_state()

    @tensor opt=(a, b, i, j) backend=eTbackend C[a, i, j, b] = A[a, i] * B[b, j]
end

function test12()
    A = "A" => ("v", "o", "v", "o") => (1, 2, 3, 4)
    B = "B" => ("v", "o") => (1, 2)
    C = "C" => ("v", "o") => (1, 2)

    reset_state()

    @tensor opt=(a, b, i, j) backend=eTbackend C[a, i] = A[a, i, b, j] * B[b, j]
end

function test13()
    A = "A" => ("v", "o", "v", "o") => (1, 2, 3, 4)
    B = "B" => ("v", "o", "v") => (1, 2, 3)
    C = "C" => ("o",) => (1,)

    reset_state()

    @tensor opt=(a, b, i, j) backend=eTbackend C[j] = B[a, i, b] * A[a, i, b, j]
end

function test14()
    A = "A" => ("v", "v", "v", "v") => (1, 2, 3, 4)
    B = "B" => ("v", "o") => (1, 2)
    C = "C" => ("o", "o", "o", "o") => (1, 2, 3, 4)

    reset_state()

    @tensor opt=(a, b, c, d, i, j, k, l) backend=eTbackend begin
        C[i, j, k, l] = A[a, b, c, d] * B[a, i] * B[b, j] * B[c, k] * B[d, l]
    end
end

function test14()
    A1 = "A1" => ("v", "o", "v", "o") => (1, 2, 3, 4)
    A2 = "A2" => ("v", "o", "v", "o") => (1, 2, 3, 4)
    B = "B" => ("v", "o") => (1, 2)
    C = "C" => ("v", "o") => (1, 2)

    reset_state()

    @tensor opt=(a, b, c, i, j, k) backend=eTbackend begin
        C[a, i] = A1[a, i, b, j] * A2[b, j, c, k] * B[c, k]
    end
end

function test15()
    A = "A" => ("v", "o", "v", "o") => (1, 2, 3, 4)
    B = "B" => ("v", "o", "v") => (1, 2, 3)
    C = "C" => ("o",) => (1,)
    D = "D" => () => ()

    reset_state()

    @tensor opt=(a, b, c, i, j, k) backend=eTbackend begin
        D[] = A[a, i, b, j] * B[a, i, b] * C[j]
    end
end

# Here it is good to sort after the contraction
function test16()
    A = "A" => ("o", "o", "v", "v", "v", "v") => (1:6...,)
    B = "B" => ("o", "o", "v", "v", "v", "v") => (1:6...,)
    C = "C" => ("o", "o", "o", "o") => (1, 2, 3, 4)

    reset_state()

    @tensor opt=(i, j, k, l, a, b, c, d) backend=eTbackend begin
        C[i, j, k, l] = A[i, j, a, b, c, d] * B[l, k, a, b, c, d]
    end
end

# Here it is good to sort before the contraction
function test17()
    A = "A" => ("o", "o", "v") => (1, 2, 3)
    B = "B" => ("o", "o", "v") => (1, 2, 3)
    C = "C" => ("o", "o", "o", "o") => (1, 2, 3, 4)

    reset_state()

    @tensor opt=(i, j, k, l, a) backend=eTbackend begin
        C[i, j, k, l] = A[i, j, a] * B[l, k, a]
    end
end

function test18()
    F = "F" => ("g", "g") => (1, 2)
    h = "h" => ("g", "g") => (1, 2)
    g_ggoo = "g_ggoo" => ("g", "g", "o", "o") => (1, 2, 3, 4)
    g_goog = "g_goog" => ("g", "o", "o", "g") => (1, 2, 3, 4)

    reset_state()

    @tensor backend=eTbackend begin
        F[p, q] = h[p, q] + 2 * g_ggoo[q, p, i, i] - g_goog[q, i, i, p]
    end
end

function test19()
    F = "F" => ("v", "o") => (1, 2)
    h = "h" => ("v", "o") => (1, 2)
    d = "d" => ("o", "o") => (1, 2)

    reset_state()

    @tensor backend=eTbackend begin
        F[a, i] = h[a, i] * d[j, j]
    end
end

function test20()
    F = "F" => ("v", "o") => (1, 2)
    d = "d" => ("v", "o") => (1, 2)

    γ = "wf%s0" => () => ()

    reset_state()

    @tensor backend=eTbackend begin
        F[a, i] = γ[] * d[a, i]
    end
end

function test21()
    E = "E" => () => ()
    d = "d" => ("o", "o") => (1, 2)

    reset_state()
    
    γ = input_scalar("gamma")

    @tensor backend=eTbackend begin
        E[] = 2 * d[i, j] * 3 * d[j, i]
    end
end

function test22()
    A = "A" => ("v", "o", "v", "o") => (1, 2, 3, 4)
    B = "B" => ("o", "v") => (1, 2)
    C = "C" => ("v", "o") => (1, 2)

    reset_state()

    @tensor opt=(a=>10χ,b=>10χ,i=>χ,j=>χ) backend=eTbackend begin
        C[a, i] = 2 * A[a, i, b, j] * B[j, b]
    end
end

function test23()
    A = "A" => ("v", "o") => (1, 2)
    B = "B" => ("v", "v") => (1, 2)
    C = "C" => ("v", "o") => (1, 2)

    reset_state()

    @tensor opt=(a=>10χ,b=>10χ,i=>χ,j=>χ) backend=eTbackend begin
        C[a, i] += A[b, i] * B[a, b]
    end
end

function test24()
    reset_state()

    rho_vovo = output_tensor("rho_vovo", ("v", "o", "v", "o"))
    cs_vovo = input_tensor("cs_vovo", ("v", "o", "v", "o"))
    d_vv = input_tensor("d_vv", ("v", "v"))

    @tensor backend = eTbackend begin
        rho_vovo[a,i,b,j] += 2*cs_vovo[a,i,c,j]*d_vv[b,c]
    end
end

function test25()
    reset_state()

    rho_vovo = output_tensor("rho_vovo", ("v", "o", "v", "o"))
    ct_vo = input_tensor("ct_vo", ("v", "o"))
    s_vo = input_tensor("s_vo", ("v", "o"))
    g_ovov = input_tensor("g_ovov", ("o", "v", "o", "v"))
    t_vovo = input_tensor("t_vovo", ("v", "o", "v", "o"))

    @tensor opt=(a=>10χ,b=>10χ,c=>10χ,d=>10χ,i=>χ,j=>χ,k=>χ,l=>χ) backend=eTbackend begin
        rho_vovo[a,i,b,j] += ct_vo[a,k]*g_ovov[k,c,l,d]*s_vo[d,j]*t_vovo[b,l,c,i]
    end
end
