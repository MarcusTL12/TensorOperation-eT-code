using Printf

eT_dim_dict::Dict{String,String} = Dict{String,String}([
    "v" => "wf%n_v",
    "o" => "wf%n_o",
    "g" => "wf%n_mo",
])

function make_eT_num(x)
    if x isa Sym
        return make_eT_num_sym(x)
    end

    if x isa String
        return x
    end

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

function make_eT_num_sym_single_term(x::Sym)
    if iszero(x)
        return "zero"
    end

    pref, rest = Sym.(x.o.primitive())

    if isone(rest)
        make_eT_num_words(pref)
    else
        rest_str = replace(string(rest), "^" => "**")

        if !isone(pref)
            "$(make_eT_num_words(pref))*$rest_str"
        else
            rest_str
        end
    end
end

function make_eT_num_sym(x::Sym)
    terms = x.o.as_terms()[1]

    io = IOBuffer()

    isfirst = true
    for (term, _) in terms
        term_sym = Sym(term)
        if isfirst
            isfirst = false
        else
            print(io, "+")
        end
        print(io, make_eT_num_sym_single_term(term_sym))
    end
    String(take!(io))
end
