macro time_calls(ex)
    @assert ex.head == :function
    body = ex.args[2].args

    d = gensym()
    new_body = map(body) do ex
        ex isa LineNumberNode && return ex
        is_call = ex isa Expr && ((ex.head == :(=) && ex.args[2].head == :call) || ex.head == :call)
        is_call || return esc(ex)

        if ex.head != :call
            f = String(ex.args[2].args[1])
            quote
                local elapsedtime = time_ns()
                $(esc(ex))
                elapsedtime = time_ns() - elapsedtime
                push!($d, $f => elapsedtime/1e9)
            end
        else
            f = String(ex.args[1])
            quote
                local elapsedtime = time_ns()
                local val = $(esc(ex))
                elapsedtime = time_ns() - elapsedtime
                push!($d, $f => elapsedtime/1e9)
                val
            end
        end
    end
    new_body = Expr(:block, new_body...)

    body = quote
        $d = Dict{String, Float64}()
        local val = quote end
        $d, val
    end
    body.args[4].args[1].args[2] = new_body

    Expr(:function, esc(ex.args[1]), body)
end
macro mem_calls(ex)
    @assert ex.head == :function
    body = ex.args[2].args

    d = gensym()
    new_body = map(body) do ex
        ex isa LineNumberNode && return ex
        is_call = ex isa Expr && ((ex.head == :(=) && ex.args[2].head == :call) || ex.head == :call)
        is_call || return esc(ex)

        if ex.head != :call
            f = String(ex.args[2].args[1])
            quote
                mempeak()
                $(esc(ex))
                peak = mempeak()
                push!($d, $f => peak)
            end
        else
            f = String(ex.args[1])
            quote
                mempeak()
                local val = $(esc(ex))
                peak = mempeak()
                push!($d, $f => peak)
                val
            end
        end
    end
    new_body = Expr(:block, new_body...)

    body = quote
        $d = Dict{String, Float64}()
        local val = quote end
        $d, val
    end
    body.args[4].args[1].args[2] = new_body

    Expr(:function, esc(ex.args[1]), body)
end

