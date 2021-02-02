# copyright imec - evaluation license - not for distribution

import Statistics: quantile

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

function count_labels(v::AbstractVector{<:Integer})
    min, max = extrema(v)
    @assert min >= 1
    max
end

function count_map(v::AbstractVector{<:Integer}, M::Integer)
    counts = zeros(Int64, M)
    @inbounds for x in v
        counts[x] += 1
    end

    counts
end

function counting_sort(v::AbstractVector{<:Integer}, M::Integer)
    counts = count_map(v, M)

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

function _cut!(labels::AbstractVector, v::AbstractVector, breaks::AbstractVector; right=true)
    lower = first(breaks)
    upper = last(breaks)

    @inbounds for i in eachindex(v)
        x = v[i]
        labels[i] = if lower <= x <= upper
            idx = searchsortedlast(breaks, x)
            if right && x == breaks[idx]
                idx -= 1
            end
            idx
        else
            0
        end
    end

    labels
end

function cut!(labels::AbstractVector, v::AbstractVector, breaks::AbstractVector; right=true)
    if issorted(breaks)
        breaks = sort(breaks)
    end

    _cut!(labels, v, breaks)
    labels
end

cut(v::AbstractVector, breaks::AbstractVector; right=true) = cut!(similar(v, Int64), v, breaks; right=right)
function cut(v::AbstractVector, nbreaks::Int64; method=:width, right=true)
    breaks = if method == :width
        min_v, max_v = extrema(v)
        dx = max_v - min_v

        breaks = collect(range(min_v, max_v, length=nbreaks+1))
        breaks[1] -= dx/1000
        breaks[end] += dx/1000
        breaks
    elseif method == :frequency
        quantile(v, range(0, 1.0, length=nbreaks))
    else
        error("unknown binning method: $method")
    end

    breaks, _cut!(similar(v, Int64), v, breaks; right=right)
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

function tiedrank(x::AbstractVector)
	n = length(x)
	J = sortperm(x)

	rk = zeros(size(J))
	i = 1
	while i <= n
		j = i
		@inbounds while (j < n) && x[J[j]] == x[J[j+1]]
			j += 1
		end

		@inbounds for k in i:j
			rk[J[k]] = (i + j) / 2.
		end

		i = j + 1
	end

	rk
end
