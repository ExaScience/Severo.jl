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
