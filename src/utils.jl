function make_unique!(out::AbstractVector{T}, names::AbstractVector{T}, sep::AbstractString=".") where {T <: AbstractString}
    seen = Dict{T, Int64}()
    n = length(names)
    dup = falses(n)

    @inbounds for i in 1:n
        x = names[i]
        if !(x in keys(seen))
            push!(seen, x=>1)
        else
            dup[i] = true
        end
    end

    @inbounds for i in 1:n
        x = names[i]
        if dup[i]
            cnt = seen[x]

            y = string(x, sep, cnt)
            while (y in keys(seen)) && cnt <= n
                cnt += 1
                y = string(x, sep, cnt)
            end

            out[i] = y
            seen[x] = cnt + 1
            push!(seen, y=>1)
        else
            out[i] = x
        end
    end

    out
end

make_unique(names::AbstractVector{T}, sep::AbstractString=".") where {T <: AbstractString} = make_unique!(similar(names), names, sep)

function counting_sort(v::AbstractVector, M::Integer)
    counts = zeros(Int64, M)
    @inbounds for x in v
        counts[x] += 1
    end

    ax = axes(v, 1)
    ix = similar(Vector{eltype(ax)}, ax)

    tot = 1
    @inbounds for (i,c) in enumerate(counts)
        counts[i] = tot
        tot += c
    end

    @inbounds for (i,x) in enumerate(v)
        j = counts[x]
        ix[j] = i
        counts[x] += 1
    end

    ix,counts
end

function rep_each(x::AbstractVector{Tv}, each::AbstractVector{Ti}) where {Tv, Ti <: Integer}
    @assert length(x) == length(each)
    r = similar(x, sum(each))

    idx = 1
    for j in eachindex(x)
        @inbounds v = x[j]
        @inbounds for i in 1:each[j]
            r[idx] = v
            idx += 1
        end
    end

    r
end
