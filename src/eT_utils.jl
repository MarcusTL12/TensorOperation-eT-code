using Printf

eT_dim_dict::Dict{String,String} = Dict{String,String}([
    "v" => "wf%n_v",
    "o" => "wf%n_o",
    "g" => "wf%n_mo",
])

function make_eT_num(x)
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
