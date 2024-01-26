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

function make_ov_dimdict(code::Union{StaticEinCode{Char},DynamicEinCode{Char}})
    dimdict = Dict{Char,String}()
    cost = uniformsize(code, 10)
    for i in keys(cost)
        dimdict[i] =
            if i in "pqrstuvw"
                "g"
            elseif i in "ijklmno"
                "o"
            elseif i in "abcdefgh"
                "v"
            end
    end
    dimdict
end

function make_ov_dimdict(code::Union{StaticEinCode{T},DynamicEinCode{T}}) where {T}
    dimdict = Dict{T,String}()
    cost = uniformsize(code, 10)
    for i in keys(cost)
        dimdict[i] = "g"
    end
    dimdict
end

function make_tmp_inputnames(code)
    curname = 'A'
    [
        begin
            name = string(curname)
            curname += 1
            (name, true)
        end for _ in getixsv(code)
    ]
end

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

    contiguous_penalty = 2

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
    output_perm::Vector{Int},
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
    outinds = @view outinds[output_perm]
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
    inputperms::Vector{Vector{Vector{Int}}},
    outputperms::Vector{Vector{Int}}) where {T,C}
    best_choice = nothing
    best_perms = nothing
    best_outperm = nothing
    best_score = nothing

    i = 0

    for outpermchoice in outputperms
        for permchoices in Iterators.product(inputperms...)
            permchoices_vec = collect(permchoices)
            for choices in choices_iter(get_num_choices(steps))
                choices_vec = collect(choices)

                score = compute_sorting_complexity(costdict, choices_vec,
                    permchoices_vec, outpermchoice, steps)

                if isnothing(best_score) || score < best_score
                    best_score = score
                    best_choice = choices_vec
                    best_perms = permchoices_vec
                    best_outperm = outpermchoice
                    @show best_score
                end

                i += 1

                if score == 0
                    println("Total number of checked permutations: $i")
                    return best_choice, best_perms, best_outperm, best_score
                end
            end
        end
    end

    println("Total number of checked permutations: $i")

    best_choice, best_perms, best_outperm, best_score
end

function make_trivial_inputperms(code)
    [[collect(eachindex(ix))] for ix in getixsv(code)]
end

function make_trivial_outperm(code)
    [collect(eachindex(getiyv(code)))]
end

function get_dimstr(dims)
    dimbuf = IOBuffer()

    isfirst = true
    for d in dims
        if isfirst
            isfirst = false
        else
            print(dimbuf, ", ")
        end
        print(dimbuf, eT_dim_dict[d])
    end

    String(take!(dimbuf))
end

function get_compact_dimstr(dims)
    if isempty(dims)
        return "1"
    end

    counts = Dict{String,Int}()
    for x in dims
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

    String(take!(dimbuf))
end

function get_dims(dimdict, inds)
    [dimdict[i] for i in inds]
end

eT_dim_dict::Dict{String,String} = Dict{String,String}([
    "v" => "wf%n_v",
    "o" => "wf%n_o",
    "g" => "wf%n_mo",
])

struct FortranFunction
    code_body::IOBuffer
    output_param::Tuple{String,Vector{String}}
    input_parameters::Vector{Tuple{String,Vector{String}}}
    local_variables::Vector{Tuple{String,Vector{String}}}
    n_integers::Ref{Int}
    use_ddot::Ref{Bool}
end

function get_intermediate_name!(func::FortranFunction, dims)
    i = length(func.local_variables) + 1
    name = "X$i"
    push!(func.local_variables, (name, dims))
    name
end

function make_code!(func::FortranFunction,
    choices::Vector{NTuple{4,Vector{Int}}},
    input_names::Vector{Tuple{String,Bool}},
    input_perms::Vector{Vector{Int}},
    output_name::String,
    output_perm::Vector{Int},
    dimdict::Dict{T,String},
    steps::Vector{Step{T}}) where {T}

    intermediates_order = Dict{Int,Vector{T}}()

    name_translation = Dict{Tens{T},String}()

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

        if name1[1]
            name_translation[(name1, collect(inds1))] = input_names[name1[2]][1]
        end

        if name2[1]
            name_translation[(name2, collect(inds2))] = input_names[name2[2]][1]
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

        (exinds1, continds1), trans1 = partition_inds(sorted_inds1, continds)
        (exinds2, continds2), trans2 = partition_inds(sorted_inds2, continds)

        outorder = if mulorder[1] == 1
            [exinds1; exinds2]
        else
            [exinds2; exinds1]
        end

        inds1 = collect(inds1)
        inds2 = collect(inds2)

        sorted_inds1 = collect(sorted_inds1)
        sorted_inds2 = collect(sorted_inds2)

        sorting1 = false
        dims1 = get_dims(dimdict, inds1)
        sorted_dims1 = get_dims(dimdict, sorted_inds1)

        if !issorted(perm1)
            sorting1 = true
            println("Allocating   $((name1, sorted_inds1))")
            intermediate_name = get_intermediate_name!(func, sorted_dims1)
            name_translation[(name1, sorted_inds1)] = intermediate_name
            println(func.code_body,
                "      call mem%alloc($intermediate_name, $(get_dimstr(sorted_dims1)))")

            println("Sorting      $((name1, inds1)) -> $((name1, sorted_inds1))")
            old_name = name_translation[(name1, inds1)]
            println(func.code_body,
                "      call sort_to_$(prod(string, perm1))($old_name, \
                $intermediate_name, $(get_dimstr(dims1)))")

            if !name1[1]
                println("Deallocating $((name1, inds1))")
                println(func.code_body,
                    "      call mem%dealloc($old_name)")
            end
        end

        sorting2 = false
        dims2 = get_dims(dimdict, inds2)
        sorted_dims2 = get_dims(dimdict, sorted_inds2)

        if !issorted(perm2)
            sorting2 = true
            println("Allocating   $((name2, sorted_inds2))")
            intermediate_name = get_intermediate_name!(func, sorted_dims2)
            name_translation[(name2, sorted_inds2)] = intermediate_name
            println(func.code_body,
                "      call mem%alloc($intermediate_name, $(get_dimstr(sorted_dims2)))")

            println("Sorting      $((name2, inds2)) -> $((name2, sorted_inds2))")
            old_name = name_translation[(name2, inds2)]
            println(func.code_body,
                "      call sort_to_$(prod(string, perm2))($old_name, \
                $intermediate_name, $(get_dimstr(dims2)))")

            if !name2[1]
                println("Deallocating $((name2, inds2))")
                println(func.code_body,
                    "      call mem%dealloc($old_name)")
            end
        end

        outdims = get_dims(dimdict, outorder)

        allocating_output = false

        if nameout[2] != 0 || outorder != (@view indsout[output_perm])
            allocating_output = true
            println("Allocating   $((nameout, outorder))")
            intermediate_name = get_intermediate_name!(func, outdims)
            name_translation[(nameout, outorder)] = intermediate_name
            println(func.code_body,
                "      call mem%alloc($intermediate_name, $(get_dimstr(outdims)))")
        end

        left_tens, left_T, left_exinds, right_tens, right_T, right_exinds =
            if mulorder[1] == 1
                (name1, sorted_inds1), !trans1, exinds1,
                (name2, sorted_inds2), trans2, exinds2
            else
                (name2, sorted_inds2), !trans2, exinds2,
                (name1, sorted_inds1), trans1, exinds1
            end

        left_exdims = get_dims(dimdict, left_exinds)
        right_exdims = get_dims(dimdict, right_exinds)
        contdims = get_dims(dimdict, continds1)

        m = get_compact_dimstr(left_exdims)
        n = get_compact_dimstr(right_exdims)
        k = get_compact_dimstr(contdims)

        lda = left_T ? k : m
        ldb = right_T ? n : k
        ldc = m

        left_name = name_translation[left_tens]
        right_name = name_translation[right_tens]
        out_name = get(name_translation, (nameout, outorder), output_name)

        α = "one"
        β = allocating_output ? "zero" : "one"

        println("Contracting  $left_tens * $right_tens \
        -> $((nameout, outorder))")

        println(
            func.code_body,
            """
!
      call dgemm('$(left_T ? 'T' : 'N')', '$(right_T ? 'T' : 'N')', &
         $m, &
         $n, &
         $k, &
         $α, &
         $left_name, &
         $lda, &
         $right_name, &
         $ldb, &
         $β, &
         $out_name, &
         $ldc)
!"""
        )

        if !name1[1] || sorting1
            println("Deallocating $((name1, sorted_inds1))")
            old_name = name_translation[(name1, sorted_inds1)]
            println(func.code_body,
                "      call mem%dealloc($old_name)")
        end

        if !name2[1] || sorting2
            println("Deallocating $((name2, sorted_inds2))")
            old_name = name_translation[(name2, sorted_inds2)]
            println(func.code_body,
                "      call mem%dealloc($old_name)")
        end

        println()

        intermediates_order[nameout[2]] = outorder
    end

    outname = last(steps)[1][1]
    outinds = last(steps)[1][2]
    outinds = outinds[output_perm]
    actual_outinds = intermediates_order[0]
    outdims = get_dims(dimdict, actual_outinds)

    outperm = get_permutation(actual_outinds, outinds)

    if !issorted(outperm)
        println("Sorting      $((outname, actual_outinds)) -> $((outname, outinds))")
        old_name = name_translation[(outname, actual_outinds)]
        println(func.code_body,
            "      call add_$(prod(string, invperm(outperm)))_to_\
            $(prod(string, 1:length(outperm)))(one, $old_name, \
            $output_name, $(get_dimstr(outdims)))")

        println("Deallocating $((outname, actual_outinds))")
        println(func.code_body,
            "      call mem%dealloc($old_name)")
    end
end
