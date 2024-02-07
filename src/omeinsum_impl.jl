using OMEinsum
using StaticArrays
using Permutations
using SymPy

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

function make_ov_cost(code::Union{StaticEinCode{Char},DynamicEinCode{Char},StaticNestedEinsum{Char}})
    cost = uniformsize(code, 10)
    for i in keys(cost)
        if i in "abcdefg"
            cost[i] *= 10
        end
    end
    cost
end

make_ov_cost(code::AbstractEinsum) = uniformsize(code, 10)

function make_ov_dimdict(code::Union{StaticEinCode{Char},DynamicEinCode{Char},StaticNestedEinsum{Char}})
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

function make_ov_dimdict(code::Union{StaticEinCode{T},DynamicEinCode{T},StaticNestedEinsum{T}}) where {T}
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
    sort_scalars_first!(steps)
    make_scaled_steps(steps)
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

function sort_scalars_first!(steps::Vector{Step{T}}) where {T}
    for i in eachindex(steps)
        out, (in1, in2) = steps[i]
        if isempty(out[2])
            last_j = i
            for j in i-1:-1:1
                (out_j, in) = steps[j]
                if out_j != in1 && out_j != in2
                    steps[j+1] = out_j => in
                    last_j = j
                else
                    break
                end
            end
            steps[last_j] = out => (in1, in2)
        end
    end
end

function bypass_scalars!(steps::Vector{Step{T}}) where {T}
    done = false
    while !done
        done = true
        for i in eachindex(steps)
            out, (in1, in2) = steps[i]
            otherin = if isempty(in1[2])
                in2
            elseif isempty(in2[2])
                in1
            end

            if !isnothing(otherin)
                j = findfirst(((x, _),) -> x == otherin, steps)

                if !isnothing(j)
                    _, in = steps[j]
                    steps[j] = out => in
                    deleteat!(steps, i)
                    done = false
                    break
                end
            end
        end
    end
end

function get_scalars(code::AbstractEinsum, inputnames)
    scalar = 1

    for ((name, _), inds) in zip(inputnames, getixsv(code))
        if isempty(inds)
            scalar *= Sym(name)
        end
    end

    scalar
end

const ScaledStep{T} = Pair{Tens{T},Tuple{Vector{Tens{T}},Tens{T},Tens{T}}}

function make_scaled_steps(steps::Vector{Step{T}}) where {T}
    newstepsdict = Dict([out => (Tens{T}[], in1, in2) for (out, (in1, in2)) in steps])

    for (out, (in1, in2)) in steps
        scalar, otherin = if isempty(in1[2]) && haskey(newstepsdict, in2)
            in1, in2
        elseif isempty(in2[2]) && haskey(newstepsdict, in1)
            in2, in1
        else
            nothing, nothing
        end

        if !isnothing(otherin)
            scalars, in... = newstepsdict[otherin]
            push!(scalars, scalar)
            delete!(newstepsdict, otherin)
            newstepsdict[out] = (scalars, in...)
        end
    end

    ScaledStep{T}[out => newstepsdict[out] for (out, _) in steps if haskey(newstepsdict, out)]
end

function partition_trace_inds(inds)
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

    nontrace_inds = [ind for ind in inds if ind ∉ trace_inds]

    nontrace_inds, trace_inds
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

function get_num_choices((out, (in1, in2))::Step)
    inds1, _ = partition_trace_inds(in1[2])
    inds2, _ = partition_trace_inds(in2[2])

    continds = intersect(inds1, inds2)

    perm1 = get_partition_perm(inds1, continds)
    perm2 = get_partition_perm(inds2, continds)

    sorted_inds1 = view(copy(inds1), perm1)
    sorted_inds2 = view(copy(inds2), perm2)

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

function get_num_choices((out, (_, in1, in2))::ScaledStep)
    get_num_choices(out => (in1, in2))
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
    steps::Vector{ScaledStep{T}}) where {T,C}
    intermediates_order = Dict{Int,Vector{T}}()

    total_cost = 0

    for ((mulorder, whichcontperm, ex1perm, ex2perm),
        ((nameout, indsout), (_, (name1, _inds1), (name2, _inds2)))) in zip(choices, steps)
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

        inds1, _ = partition_trace_inds(inds1)
        inds2, _ = partition_trace_inds(inds2)

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

function optimize_choices(costdict::Dict{T,C}, steps::Vector{ScaledStep{T}},
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
    local_variables::Vector{Tuple{String,Int}}
    n_integers::Ref{Int}
    use_ddot::Ref{Bool}

    function FortranFunction(outname)
        new(IOBuffer(), outname,
            Tuple{String,Vector{String}}[], Tens[], Ref(0), Ref(false))
    end

    function FortranFunction(outname::String)
        new(IOBuffer(), (outname, String[]),
            Tuple{String,Vector{String}}[], Tens[], Ref(0), Ref(false))
    end
end

function get_intermediate_name!(func::FortranFunction, dims)
    i = length(func.local_variables) + 1
    name = "X$i"
    push!(func.local_variables, (name, length(dims)))
    name
end

function make_code!(func::FortranFunction,
    choices::Vector{NTuple{4,Vector{Int}}},
    input_names::Vector{Tuple{String,Bool}},
    input_perms::Vector{Vector{Int}},
    output_name::String,
    output_perm::Vector{Int},
    dimdict::Dict{T,String},
    scalar,
    steps::Vector{ScaledStep{T}}) where {T}

    intermediates_order = Dict{Int,Vector{T}}()

    name_translation = Dict{Tens{T},String}()

    has_output = false

    for ((mulorder, whichcontperm, ex1perm, ex2perm),
        ((nameout, indsout), (scalars, (name1, _inds1), (name2, _inds2)))) in
        zip(choices, steps)
        full_inds1 = if name1[1]
            name, is_inp = input_names[name1[2]]
            if is_inp
                dims = get_dims(dimdict, _inds1)
                if (name, dims) ∉ func.input_parameters
                    push!(func.input_parameters, (name, dims))
                end
            end
            @view _inds1[input_perms[name1[2]]]
        else
            intermediates_order[name1[2]]
        end

        full_inds2 = if name2[1]
            name, is_inp = input_names[name2[2]]
            if is_inp
                dims = get_dims(dimdict, _inds2)
                if (name, dims) ∉ func.input_parameters
                    push!(func.input_parameters, (name, dims))
                end
            end
            @view _inds2[input_perms[name2[2]]]
        else
            intermediates_order[name2[2]]
        end

        if name1[1]
            name_translation[(name1, collect(full_inds1))] = input_names[name1[2]][1]
        end

        if name2[1]
            name_translation[(name2, collect(full_inds2))] = input_names[name2[2]][1]
        end

        for (s_name, s_inds) in scalars
            if s_name[1]
                i_name, is_inp = input_names[s_name[2]]
                name_translation[(s_name, s_inds)] = input_names[s_name[2]][1]
                name_entry = (i_name, String[])
                if is_inp && name_entry ∉ func.input_parameters
                    push!(func.input_parameters, name_entry)
                end
            end
        end

        inds1, trace_inds1 = partition_trace_inds(full_inds1)
        inds2, trace_inds2 = partition_trace_inds(full_inds2)

        tracing1 = false
        if !isempty(trace_inds1)
            tracing1 = true
            dims1 = get_dims(dimdict, inds1)

            intermediate_name = get_intermediate_name!(func, dims1)
            name_translation[(name1, inds1)] = intermediate_name
            if !isempty(inds1)
                println("Allocating   $((name1, inds1))")
                println(func.code_body,
                    "      call mem%alloc($intermediate_name, $(get_dimstr(dims1)))")
            end

            println("Tracing      $((name1, collect(full_inds1))) -> $((name1, inds1))")
            old_name = name_translation[(name1, full_inds1)]
            make_trace_code!(func, dimdict, old_name, full_inds1, intermediate_name, inds1, trace_inds1, true)
        end

        tracing2 = false
        if !isempty(trace_inds2)
            tracing2 = true
            dims2 = get_dims(dimdict, inds2)

            intermediate_name = get_intermediate_name!(func, dims2)
            name_translation[(name2, inds2)] = intermediate_name
            if !isempty(inds2)
                println("Allocating   $((name2, inds2))")
                println(func.code_body,
                    "      call mem%alloc($intermediate_name, $(get_dimstr(dims2)))")
            end

            println("Tracing      $((name1, collect(full_inds2))) -> $((name1, inds2))")
            old_name = name_translation[(name2, full_inds2)]
            make_trace_code!(func, dimdict, old_name, full_inds2, intermediate_name, inds2, trace_inds2, true)
        end

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

        @assert !(sorting1 && tracing1) "Simultaneous sorting and tracing"
        @assert !(sorting2 && tracing2) "Simultaneous sorting and tracing"

        # Doing contraction

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

        left_name = name_translation[left_tens]
        right_name = name_translation[right_tens]

        if isempty(left_exinds) && isempty(continds)
            name_translation[(nameout, outorder)] = right_name
            scalar *= Sym(left_name)
        elseif isempty(right_exinds) && isempty(continds)
            name_translation[(nameout, outorder)] = left_name
            scalar *= Sym(right_name)
        else
            outdims = get_dims(dimdict, outorder)

            allocating_output = false

            if (nameout[2] != 0 || outorder != (@view indsout[output_perm])) &&
               (!isempty(inds1) || !isempty(inds2))
                allocating_output = true
                intermediate_name = get_intermediate_name!(func, outdims)
                name_translation[(nameout, outorder)] = intermediate_name

                if !isempty(outorder)
                    println("Allocating   $((nameout, outorder))")
                    println(func.code_body,
                        "      call mem%alloc($intermediate_name, $(get_dimstr(outdims)))")
                end
            end

            out_name = get(name_translation, (nameout, outorder), output_name)

            α = scalar * if !isempty(scalars)
                prod(Sym(name_translation[s]) for s in scalars)
            else
                1
            end
            scalar = 1

            β = allocating_output ? "zero" : "one"

            if !allocating_output
                has_output = true
            end

            println("Contracting  $left_tens * $right_tens \
            -> $((nameout, outorder))")

            func.use_ddot[] |= simplify_dgemm(func.code_body, left_T, right_T, left_exdims,
                right_exdims, contdims, α, left_name, right_name, β, out_name)

            # Contraction done

            if !isempty(sorted_inds1) && (!name1[1] || sorting1 || tracing1)
                println("Deallocating $((name1, sorted_inds1))")
                old_name = name_translation[(name1, sorted_inds1)]
                println(func.code_body,
                    "      call mem%dealloc($old_name)")
            end

            if !isempty(sorted_inds2) && (!name2[1] || sorting2 || tracing2)
                println("Deallocating $((name2, sorted_inds2))")
                old_name = name_translation[(name2, sorted_inds2)]
                println(func.code_body,
                    "      call mem%dealloc($old_name)")
            end
        end

        println()

        intermediates_order[nameout[2]] = outorder
    end

    outname, outinds = last(steps)[1]
    outinds = outinds[output_perm]
    actual_outinds = intermediates_order[0]
    outdims = get_dims(dimdict, outinds)

    outperm = get_permutation(actual_outinds, outinds)

    if !issorted(outperm)
        println("Sorting      $((outname, actual_outinds)) -> $((outname, outinds))")
        old_name = get(name_translation, (outname, actual_outinds), output_name)
        println(func.code_body,
            "      call add_$(prod(string, invperm(outperm)))_to_\
            $(prod(string, 1:length(outperm)))($(make_eT_num(scalar)), $old_name, \
            $output_name, $(get_dimstr(outdims)))")

        if !any(name == old_name for (name, _) in input_names)
            println("Deallocating $((outname, actual_outinds))")
            println(func.code_body,
                "      call mem%dealloc($old_name)")
        end
    elseif !has_output
        println("Output is scaled input")
        old_name = name_translation[(outname, actual_outinds)]
        if !isempty(outinds)
            dimstr = get_compact_dimstr(outdims)
            println(func.code_body,
                "      call daxpy($dimstr, $(make_eT_num(scalar)), $old_name, 1, $output_name, 1)")
        else
            if isone(scalar)
                println(func.code_body,
                    "      $output_name = $output_name + $old_name")
            elseif isone(-scalar)
                println(func.code_body,
                    "      $output_name = $output_name - $old_name")
            else
                println(func.code_body,
                    "      $output_name = $output_name + $(make_eT_num(scalar)) * $old_name")
            end
        end
    end
end

function simplify_dgemm(io::IO, leftT, rightT, leftexdims, rightexdims,
    contdims, α, leftname, rightname, β, outname)
    m = get_compact_dimstr(leftexdims)
    n = get_compact_dimstr(rightexdims)
    k = get_compact_dimstr(contdims)

    lda = leftT ? k : m
    ldb = rightT ? n : k
    ldc = m

    beta = false

    if β == "one"
        beta = true
    elseif β != "zero"
        throw("β is $β")
    end

    use_ddot = false

    if k == "1"
        if m == "1" && n == "1"
            println("Scalar scalar multiplication")
            print(io, "      $outname = ")
            if beta
                print(io, "$outname")
            end
            if isone(α)
                print(io, " + ")
            elseif isone(-α)
                print(io, " - ")
            else
                print(io, " + $(make_eT_num(α)) * ")
            end

            println(io, "$leftname * $rightname")
        elseif n == "1"
            println("daxpy")
            if !beta
                println(io, "      call zero_array($outname, $m)")
            end
            println(io,
                "      call daxpy($m, $(make_eT_num(α * Sym(rightname))), \
                $leftname, 1, $outname, 1)")
        elseif m == "1"
            println("daxpy")
            if !beta
                println(io, "      call zero_array($outname, $n)")
            end
            println(io,
                "      call daxpy($n, $(make_eT_num(α * Sym(leftname))), \
                $rightname, 1, $outname, 1)")
        else
            println("dger")
            println(
                io,
                """
!
      call dger($m, &
         $n, &
         $(make_eT_num(α)), &
         $leftname, 1, &
         $rightname, 1, &
         $outname, &
         $ldc)
!"""
            )
        end
    elseif m == "1" && n == "1"
        println("ddot")
        use_ddot = true
        print(io, "      $outname = ")
        if beta
            print(io, "$outname")

            if isone(α)
                print(io, " + ")
            elseif isone(-α)
                print(io, " - ")
            else
                print(io, " + $(make_eT_num(α)) * ")
            end
        else
            if isone(-α)
                print(io, " -")
            elseif !isone(α)
                print(io, "$(make_eT_num(α)) * ")
            end
        end

        println(io, "ddot($k, $leftname, 1, $rightname, 1)")
    elseif n == "1"
        println("dgemv")
        println(
            io,
            """
!
      call dgemv('$(leftT ? 'T' : 'N')', &
         $(leftT ? k : m), &
         $(leftT ? m : k), &
         $(make_eT_num(α)), &
         $leftname, &
         $lda, &
         $rightname, 1, &
         $β, &
         $outname, 1)
!"""
        )
    elseif m == "1"
        println("dgemv")
        println(
            io,
            """
!
      call dgemv('$(rightT ? 'N' : 'T')', &
         $(rightT ? n : k), &
         $(rightT ? k : n), &
         $(make_eT_num(α)), &
         $rightname, &
         $ldb, &
         $leftname, 1, &
         $β, &
         $outname, 1)
!"""
        )
    else
        println("dgemm")

        println(
            io,
            """
!
      call dgemm('$(leftT ? 'T' : 'N')', '$(rightT ? 'T' : 'N')', &
         $m, &
         $n, &
         $k, &
         $(make_eT_num(α)), &
         $leftname, &
         $lda, &
         $rightname, &
         $ldb, &
         $β, &
         $outname, &
         $ldc)
!"""
        )
    end

    use_ddot
end

function make_eT_num(x)
    make_eT_num(Sym(x))
end

function make_eT_num(x::Sym)
    if iszero(x)
        return "zero"
    end

    pref, rest = Sym.(x.o.primitive())

    if isone(rest)
        make_eT_num_words(pref)
    elseif isone(-rest)
        "-" * make_eT_num_words(pref)
    else
        rest_str = replace(string(rest), "^" => "**")

        if startswith(rest_str, "-")
            rest_str = rest_str[2:end]
            "-"
        else
            ""
        end *
        if !isone(pref)
            "$(make_eT_num_words(pref))*$rest_str"
        else
            rest_str
        end
    end
end

function make_eT_num_words(x)
    isneg = false
    if x isa Float64 || x isa Integer
        isneg = x < 0
        if isneg
            x = -x
        end
        x = round(x, digits=14)
    end
    if isneg
        "-"
    else
        ""
    end *
    if iszero(x)
        "zero"
    elseif isone(x)
        "one"
    elseif x == 2
        "two"
    elseif x == 3
        "three"
    elseif x == 4
        "four"
    elseif x == 1 // 2
        "half"
    else
        "$(x)d0"
    end
end

function make_trace_code!(func::FortranFunction, dimdict, from_name, from_inds,
    to_name, to_inds, trace_inds, set_zero)
    n_loop_indices = length(to_inds) + length(trace_inds)

    func.n_integers[] = max(func.n_integers[], n_loop_indices)

    used_inds = 0
    input_inds = zeros(Int, length(from_inds))
    output_inds = zeros(Int, length(to_inds))

    index_dims = String[]

    seen_inds = eltype(from_inds)[]

    for i in from_inds
        if i ∉ seen_inds
            push!(seen_inds, i)
            used_inds += 1
            push!(index_dims, eT_dim_dict[dimdict[i]])
            for j in eachindex(input_inds)
                if from_inds[j] == i
                    input_inds[j] = used_inds
                end
            end

            for j in eachindex(output_inds)
                if to_inds[j] == i
                    output_inds[j] = used_inds
                end
            end
        end
    end

    outdims = get_dims(dimdict, to_inds)

    if set_zero
        if isempty(to_inds)
            println(func.code_body, "      $to_name = zero")
        else
            println(func.code_body,
                "      call zero_array($to_name, $(get_compact_dimstr(outdims)))")
        end
    end

    tab_level = 2

    println(func.code_body, "!")
    # print(func.code_body, "!\$omp parallel do schedule(static)")
    # if n_loop_indices > 1
    #     print(func.code_body, " collapse($n_loop_indices)")
    # end
    # print(func.code_body, " private(")

    # isfirst = true
    # for i in 1:n_loop_indices
    #     if isfirst
    #         isfirst = false
    #     else
    #         print(func.code_body, ",")
    #     end
    #     print(func.code_body, "i$i")
    # end
    # println(func.code_body, ")")

    for i in reverse(1:n_loop_indices)
        println(func.code_body, "   "^tab_level, "do i$i = 1, ", index_dims[i])
        tab_level += 1
    end

    out_string = IOBuffer()

    print(out_string, "$to_name")
    if !isempty(to_inds)
        print(out_string, "(")
        isfirst = true
        for i in output_inds
            if isfirst
                isfirst = false
            else
                print(out_string, ",")
            end
            print(out_string, "i$i")
        end
        print(out_string, ")")
    end

    out_string = String(take!(out_string))

    print(func.code_body, "   "^tab_level, "$out_string = $out_string + $from_name(")

    isfirst = true
    for i in input_inds
        if isfirst
            isfirst = false
        else
            print(func.code_body, ",")
        end
        print(func.code_body, "i$i")
    end
    println(func.code_body, ")")

    for _ in 1:n_loop_indices
        tab_level -= 1
        println(func.code_body, "   "^tab_level, "end do")
    end

    # println(func.code_body, "!\$omp end parallel do")
    println(func.code_body, "!")
end

function finalize_eT_function(func::FortranFunction, routine_name, wf_type)
    io = IOBuffer()

    routine_name = "$(routine_name)_$(wf_type)"

    sort!(func.input_parameters)

    print(io, "   subroutine $routine_name(wf")

    print(io, ", ", func.output_param[1])

    for param in func.input_parameters
        print(io, ", ", param[1])
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

    print(io, "      real(dp), ")

    if !isempty(func.output_param[2])
        print(io, "dimension(")
        isfirst = true
        for d in func.output_param[2]
            if isfirst
                isfirst = false
            else
                print(io, ",")
            end

            print(io, eT_dim_dict[d])
        end
        print(io, "), ")
    end

    println(io, "intent(inout) :: ", func.output_param[1], "\n!")

    structure_dict = Dict()

    for (name, structure) in func.input_parameters
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

        print(io, ", intent(in) :: ")

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

    dim_dict = Dict()

    for (name, n) in func.local_variables
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

    if !isempty(func.local_variables)
        println(io, "!")
    end

    if func.n_integers[] > 0
        print(io, "      integer :: ")
        isfirst = true
        for i in 1:func.n_integers[]
            if isfirst
                isfirst = false
            else
                print(io, ", ")
            end
            print(io, "i$i")
        end
        println(io, "\n!")
    end

    if func.use_ddot[]
        println(io, "      real(dp), external :: ddot\n!")
    end

    print(io, String(take!(func.code_body)))

    println(io, "!\n   end subroutine $routine_name")

    String(take!(io))
end

function make_single_tensor_code!(func::FortranFunction, code, name,
    dimdict, prefactor)
    from_inds = only(getixsv(code))
    to_inds = getiyv(code)

    inds, trace_inds = partition_trace_inds(from_inds)

    @assert to_inds == inds "Needs to sort trace"

    if name[2]
        dims = get_dims(dimdict, from_inds)
        if (name[1], dims) ∉ func.input_parameters
            push!(func.input_parameters, (name[1], dims))
        end
    end

    make_trace_code!(func, dimdict, name[1], from_inds, func.output_param[1],
        to_inds, trace_inds, false)
end

function update_code!(func::FortranFunction, code, prefactor, names_perms,
    outperms=make_trivial_outperm(code))

    names = [(name, isinp) for (name, isinp, _...) in names_perms]
    inputperms = make_trivial_inputperms(code)
    for (i, inp) in enumerate(names_perms)
        if length(inp) >= 3
            inputperms[i] = inp[3]
        end
    end

    dimdict = make_ov_dimdict(code)

    if length(getixsv(code)) == 1
        return make_single_tensor_code!(func, code, names[1],
            dimdict, prefactor)
    end

    cost = make_ov_cost(code)

    optcode = optimize_code(code, cost, GreedyMethod())

    steps = walk_einsums(optcode)

    choices, perms, outperm, _ =
        optimize_choices(cost, steps, inputperms, outperms)

    make_code!(func, choices, names, perms, func.output_param[1], outperm,
        dimdict, prefactor, steps)
end
