# Severo: a software package for analysis and exploration of single-cell RNA-seq datasets.
# Copyright (c) 2021 imec vzw.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version, and Additional Terms
# (see below).

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Affero General Public License for more details.
function is_short_function_def(ex)
    ex.head === :(=) || return false
    while length(ex.args) >= 1 && isa(ex.args[1], Expr)
        (ex.args[1].head === :call) && return true
        (ex.args[1].head === :where || ex.args[1].head === :(::)) || return false
        ex = ex.args[1]
    end
    return false
end

argument_names(v::Vector) = map(argument_names, v)
argument_names(x::Symbol) = x
function argument_names(ex::Expr)
    if ex.head === :(::)
        ex.args[1]
    elseif ex.head === :kw
        argument_names(ex.args[1])
    elseif ex.head === :(...)
        Expr(:(...), argument_names(ex.args[1]))
    else
        error("unknown argument syntax: $ex")
    end
end

"""
    @partial fdef

Automatically creates a partial version of the function leaving out the first argument.
"""
macro partial(fdef)
    if isa(fdef, Expr) && (fdef.head === :function || is_short_function_def(fdef))
        header = fdef.args[1]

        headn = if header.head === :(::)
            header.args[1]
        else
            header
        end

        ex = headn = copy(headn)
        while ex.head !== :call
            ex.head === :where || error("unknown function header: $ex")
            ex.args[1] = copy(ex.args[1])
            ex = ex.args[1]
        end

        fargs = ex.args
        argnames, kwnames = if isa(fargs[2], Expr) && fargs[2].head == :parameters
            deleteat!(fargs, 3)
            argument_names(fargs[3:end]), argument_names(fargs[2].args)
        else
            deleteat!(fargs, 2)
            argument_names(fargs[2:end]), Symbol[]
        end

        f = esc(ex.args[1])
        argnames = esc.(argnames)
        kwnames = esc.(kwnames)

        docstr =
        """
            $headn

        Partial version of $header
        """

        quote
            @doc $docstr
            $(esc(headn)) = x -> $f(x, $(argnames...); $(kwnames...))
            Base.@__doc__ $(esc(fdef))
        end
    else
        error("invalid syntax; @partial must be used with a function definition")
    end
end
