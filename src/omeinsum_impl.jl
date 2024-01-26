using OMEinsum
using StaticArrays
using Permutations

function get_partition_perm(inds, continds)
    function isleftind(i)
        !((inds[1] ∈ continds) ⊻ (i ∈ continds))
    end

    order = zeros(Int, length(inds))

    j = 1

    for (i, ind) in enumerate(inds)
        if isleftind(ind)
            order[i] = j
            j += 1
        end
    end

    for i in eachindex(order)
        if iszero(order[i])
            order[i] = j
            j += 1
        end
    end

    invperm(order)
end

function make_ov_cost(code::Union{StaticEinCode{Char},DynamicEinCode{Char}})
    cost = uniformsize(code, 10)
    for i in keys(cost)
        if i in "abcdefg"
            cost[i] *= 10
        end
    end
    cost
end

make_ov_cost(code::AbstractEinsum) = uniformsize(code, 10)

const Tens{T} = Tuple{Tuple{Bool,Int64},Vector{T}} where {T}
const Step{T} = Pair{Tens{T},NTuple{2,Tens{T}}}

function walk_einsums!(steps, code::DynamicNestedEinsum, n_intermediates=Ref(0))
    if code.tensorindex == -1
        name = (false, n_intermediates[])
        n_intermediates[] += 1

        inds = getiyv(code.eins)

        arg1, arg2 = code.args

        name1 = walk_einsums!(steps, arg1, n_intermediates)
        name2 = walk_einsums!(steps, arg2, n_intermediates)

        inds1, inds2 = getixsv(code.eins)

        push!(steps, (name, inds) => ((name1, inds1), (name2, inds2)))

        name
    else
        (true, code.tensorindex)
    end
end

function walk_einsums(code::DynamicNestedEinsum{T}) where {T}
    steps = Step{T}[]
    walk_einsums!(steps, code)
    steps
end

function walk_einsums(code)
    walk_einsums(code.eins)
end

function Base.show(io::IO, tens::Tens)
    (is_inp, ind), inds = tens

    if is_inp
        print(io, 'I')
    else
        print(io, 'X')
    end

    print(io, ind)

    if !isempty(inds)
        print(io, '_', OMEinsum._join(inds))
    end
end

function get_trace_inds(inds)
    T = eltype(inds)
    seen = T[]
    trace_inds = T[]

    for ind in inds
        if ind ∈ seen && ind ∉ trace_inds
            push!(trace_inds, ind)
        else
            push!(seen, ind)
        end
    end

    trace_inds
end

function get_permutation(from, to)
    [findfirst(==(x), from) for x in to]
end

function splitat(vec, i)
    (@view vec[1:i]), (@view vec[i+1:end])
end

function partition_inds(inds, continds)
    n_exinds = length(inds) - length(continds)

    if !isempty(continds) && inds[1] ∈ continds
        reverse(splitat(inds, length(continds))), false
    else
        splitat(inds, n_exinds), true
    end
end

function parse_contraction((out, (in1, in2))::Step)
    @assert isempty(get_trace_inds(in1[2])) "Trace not implemented yet"
    @assert isempty(get_trace_inds(in2[2])) "Trace not implemented yet"

    continds = intersect(in1[2], in2[2])

    perm1 = get_partition_perm(in1[2], continds)
    perm2 = get_partition_perm(in2[2], continds)

    @show perm1 perm2

    sorted_inds1 = view(copy(in1[2]), perm1)
    sorted_inds2 = view(copy(in2[2]), perm2)

    @show String(sorted_inds1) String(sorted_inds2)

    (exinds1, continds1), trans1 = partition_inds(sorted_inds1, continds)
    (exinds2, continds2), trans2 = partition_inds(sorted_inds2, continds)

    @show String(exinds1) String(exinds2) String(continds1) String(continds2) trans1 trans2

    # Just permute left continds
    if continds1 != continds2
        contperm = get_permutation(continds1, continds2)
        largeperm = collect(1:length(perm1))
        if trans1
            contperm .+= length(exinds1)
            largeperm[length(exinds1)+1:end] .= contperm
        else
            largeperm[1:length(continds)] .= contperm
        end

        permute!(perm1, largeperm)
    end

    println()

    sorted_inds1 = view(copy(in1[2]), perm1)
    sorted_inds2 = view(copy(in2[2]), perm2)

    @show String(sorted_inds1) String(sorted_inds2)

    (exinds1, continds1), trans1 = partition_inds(sorted_inds1, continds)
    (exinds2, continds2), trans2 = partition_inds(sorted_inds2, continds)
    @show String(exinds1) String(exinds2) String(continds1) String(continds2) trans1 trans2

    println()

    @show perm1 perm2

    nothing
end

function get_num_choices((out, (in1, in2))::Step)
    @assert isempty(get_trace_inds(in1[2])) "Trace not implemented yet"
    @assert isempty(get_trace_inds(in2[2])) "Trace not implemented yet"

    continds = intersect(in1[2], in2[2])

    perm1 = get_partition_perm(in1[2], continds)
    perm2 = get_partition_perm(in2[2], continds)

    sorted_inds1 = view(copy(in1[2]), perm1)
    sorted_inds2 = view(copy(in2[2]), perm2)

    (exinds1, continds1), _ = partition_inds(sorted_inds1, continds)
    (exinds2, continds2), _ = partition_inds(sorted_inds2, continds)

    mulorder = if isempty(exinds1) || isempty(exinds2)
        1
    else
        2
    end

    contperm = if continds1 != continds2
        2
    else
        1
    end

    (mulorder, contperm, length(exinds1), length(exinds2))
end

function get_num_choices(steps::Vector)
    get_num_choices.(steps)
end

function compute_sorting_complexity(costdict, inds, perm)
    if isempty(inds)
        return 0
    end

    cost = prod(costdict[i] for i in inds)

    contiguous_penalty = 1

    if issorted(perm)
        0
    elseif perm[1] == 1
        cost
    else
        contiguous_penalty * cost
    end
end

function compute_sorting_complexity(costdict::Dict{T,C},
    choices::Vector{NTuple{4,Vector{Int}}},
    input_perms::Vector{Vector{Int}},
    steps::Vector{Step{T}}) where {T,C}
    intermediates_order = Dict{Int,Vector{T}}()

    total_cost = 0

    for ((mulorder, whichcontperm, ex1perm, ex2perm),
        ((nameout, indsout), ((name1, _inds1), (name2, _inds2)))) in zip(choices, steps)
        inds1 = if name1[1]
            @view _inds1[input_perms[name1[2]]]
        else
            intermediates_order[name1[2]]
        end

        inds2 = if name2[1]
            @view _inds2[input_perms[name2[2]]]
        else
            intermediates_order[name2[2]]
        end

        @assert isempty(get_trace_inds(inds1)) "Trace not implemented yet"
        @assert isempty(get_trace_inds(inds2)) "Trace not implemented yet"

        continds = intersect(inds1, inds2)

        perm1 = get_partition_perm(inds1, continds)
        perm2 = get_partition_perm(inds2, continds)

        partitioned_inds1 = @view inds1[perm1]
        partitioned_inds2 = @view inds2[perm2]

        (exinds1, continds1), trans1 = partition_inds(partitioned_inds1, continds)
        (exinds2, continds2), trans2 = partition_inds(partitioned_inds2, continds)

        if continds1 != continds2
            if whichcontperm[1] == 1
                contperm = get_permutation(continds1, continds2)
                largeperm = collect(1:length(perm1))
                if trans1
                    contperm .+= length(exinds1)
                    largeperm[length(exinds1)+1:end] .= contperm
                else
                    largeperm[1:length(continds)] .= contperm
                end
                permute!(perm1, largeperm)
            else
                contperm = get_permutation(continds2, continds1)
                largeperm = collect(1:length(perm2))
                if trans2
                    contperm .+= length(exinds2)
                    largeperm[length(exinds2)+1:end] .= contperm
                else
                    largeperm[1:length(continds)] .= contperm
                end
                permute!(perm2, largeperm)
            end
        end

        largeperm = collect(1:length(perm1))
        if trans1
            largeperm[1:length(exinds1)] .= ex1perm
        else
            ex1perm_off = ex1perm .+ length(continds)
            largeperm[length(continds)+1:end] .= ex1perm_off
        end
        permute!(perm1, largeperm)

        largeperm = collect(1:length(perm2))
        if trans2
            largeperm[1:length(exinds2)] .= ex2perm
        else
            ex2perm_off = ex2perm .+ length(continds)
            largeperm[length(continds)+1:end] .= ex2perm_off
        end
        permute!(perm2, largeperm)

        sorted_inds1 = @view inds1[perm1]
        sorted_inds2 = @view inds2[perm2]

        (exinds1, continds1), _ = partition_inds(sorted_inds1, continds)
        (exinds2, continds2), _ = partition_inds(sorted_inds2, continds)

        total_cost += compute_sorting_complexity(costdict, inds1, perm1)
        total_cost += compute_sorting_complexity(costdict, inds2, perm2)

        outorder = if mulorder[1] == 1
            [exinds1; exinds2]
        else
            [exinds2; exinds1]
        end

        intermediates_order[nameout[2]] = outorder
    end

    outinds = last(steps)[1][2]
    actual_outinds = intermediates_order[0]

    outperm = get_permutation(actual_outinds, outinds)

    total_cost += compute_sorting_complexity(costdict, outinds, outperm)

    total_cost
end

function choices_iter(num_choices::NTuple{4,Int})
    Iterators.product(((p.data for p in PermGen(n)) for n in num_choices)...)
end

function choices_iter(num_choices::Vector{NTuple{4,Int}})
    Iterators.product((choices_iter(nc) for nc in num_choices)...)
end

function optimize_choices(costdict::Dict{T,C}, steps::Vector{Step{T}},
    inputperms::Vector{Vector{Vector{Int}}}) where {T,C}
    best_choice = nothing
    best_perms = nothing
    best_score = nothing

    for permchoices in Iterators.product(inputperms...)
        permchoices_vec = collect(permchoices)
        for choices in choices_iter(get_num_choices(steps))
            choices_vec = collect(choices)

            score = compute_sorting_complexity(costdict, choices_vec, permchoices_vec, steps)

            if isnothing(best_score) || score < best_score
                best_score = score
                best_choice = choices_vec
                best_perms = permchoices_vec
            end
        end
    end

    best_choice, best_perms, best_score
end

function make_trivial_inputperms(code)
    [[collect(eachindex(ix))] for ix in getixsv(code)]
end

function make_code(choices::Vector{NTuple{4,Vector{Int}}},
    input_perms::Vector{Vector{Int}},
    steps::Vector{Step{T}}) where {T}
    intermediates_order = Dict{Int,Vector{T}}()

    for ((mulorder, whichcontperm, ex1perm, ex2perm),
        ((nameout, indsout), ((name1, _inds1), (name2, _inds2)))) in zip(choices, steps)
        inds1 = if name1[1]
            @view _inds1[input_perms[name1[2]]]
        else
            intermediates_order[name1[2]]
        end

        inds2 = if name2[1]
            @view _inds2[input_perms[name2[2]]]
        else
            intermediates_order[name2[2]]
        end

        @assert isempty(get_trace_inds(inds1)) "Trace not implemented yet"
        @assert isempty(get_trace_inds(inds2)) "Trace not implemented yet"

        continds = intersect(inds1, inds2)

        perm1 = get_partition_perm(inds1, continds)
        perm2 = get_partition_perm(inds2, continds)

        partitioned_inds1 = @view inds1[perm1]
        partitioned_inds2 = @view inds2[perm2]

        (exinds1, continds1), trans1 = partition_inds(partitioned_inds1, continds)
        (exinds2, continds2), trans2 = partition_inds(partitioned_inds2, continds)

        if continds1 != continds2
            if whichcontperm[1] == 1
                contperm = get_permutation(continds1, continds2)
                largeperm = collect(1:length(perm1))
                if trans1
                    contperm .+= length(exinds1)
                    largeperm[length(exinds1)+1:end] .= contperm
                else
                    largeperm[1:length(continds)] .= contperm
                end
                permute!(perm1, largeperm)
            else
                contperm = get_permutation(continds2, continds1)
                largeperm = collect(1:length(perm2))
                if trans2
                    contperm .+= length(exinds2)
                    largeperm[length(exinds2)+1:end] .= contperm
                else
                    largeperm[1:length(continds)] .= contperm
                end
                permute!(perm2, largeperm)
            end
        end

        largeperm = collect(1:length(perm1))
        if trans1
            largeperm[1:length(exinds1)] .= ex1perm
        else
            ex1perm_off = ex1perm .+ length(continds)
            largeperm[length(continds)+1:end] .= ex1perm_off
        end
        permute!(perm1, largeperm)

        largeperm = collect(1:length(perm2))
        if trans2
            largeperm[1:length(exinds2)] .= ex2perm
        else
            ex2perm_off = ex2perm .+ length(continds)
            largeperm[length(continds)+1:end] .= ex2perm_off
        end
        permute!(perm2, largeperm)

        sorted_inds1 = @view inds1[perm1]
        sorted_inds2 = @view inds2[perm2]

        (exinds1, continds1), _ = partition_inds(sorted_inds1, continds)
        (exinds2, continds2), _ = partition_inds(sorted_inds2, continds)

        outorder = if mulorder[1] == 1
            [exinds1; exinds2]
        else
            [exinds2; exinds1]
        end

        sorting1 = false
        sorting2 = false

        if !issorted(perm1)
            sorting1 = true
            println("Allocating   $((name1, collect(sorted_inds1)))")
            println("Sorting      $((name1, inds1)) -> $((name1, collect(sorted_inds1)))")
            if !name1[1]
                println("Deallocating $((name1, inds1))")
            end
        end

        if !issorted(perm2)
            sorting2 = true
            println("Allocating   $((name2, collect(sorted_inds2)))")
            println("Sorting      $((name2, inds2)) -> $((name2, collect(sorted_inds2)))")
            if !name2[1]
                println("Deallocating $((name2, inds2))")
            end
        end

        if nameout[2] != 0 || outorder != indsout
            println("Allocating   $((nameout, outorder))")
        end

        if mulorder[1] == 1
            print("Contracting  $((name1, collect(sorted_inds1))) * $((name2, collect(sorted_inds2)))")
        else
            print("Contracting  $((name2, collect(sorted_inds2))) * $((name1, collect(sorted_inds1)))")
        end

        println(" -> $((nameout, outorder))")

        if !name1[1] || sorting1
            println("Deallocating $((name1, collect(sorted_inds1)))")
        end

        if !name2[1] || sorting2
            println("Deallocating $((name2, collect(sorted_inds2)))")
        end

        println()

        intermediates_order[nameout[2]] = outorder
    end

    outname = last(steps)[1][1]
    outinds = last(steps)[1][2]
    actual_outinds = intermediates_order[0]

    outperm = get_permutation(actual_outinds, outinds)

    if !issorted(outperm)
        println("Sorting     $((outname, actual_outinds)) -> $((outname, outinds))")
        println("Deallocating $((outname, actual_outinds))")
    end
end
