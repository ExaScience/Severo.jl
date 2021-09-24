# copyright imec - evaluation license - not for distribution

using Libdl

function determine_blas_interface(libblas::AbstractString)
    Libdl.dlopen(libblas) do lib
        suffixes = [
            # LP64-mangling suffixes: No underscore, single underscore, double underscore
            "", "_", "__",
            # ILP64-mangling suffixes: the same, but with the suffix prepended as well:
            "64_", "_64_", "_64__", "__64___"
        ]
        s = findfirst(suffixes) do suff
            dlsym(lib, string("isamax", suff); throw_error=false) !== nothing
        end
        s !== nothing || error("unable to find isamax symbol in $libopenblas")

        symname = string("isamax", suffixes[s])
        symaddr = dlsym(lib, symname)

        # Purposefully incorrect length that is `3` on 32-bit, but massively negative on 64-bit
        n = 0xffffffff00000003 % Int64
        X = [1.f0, 2.f0, 1.f0]
        incx = 1

        max_idx = ccall(symaddr, Int64, (Ref{Int64}, Ptr{Float32}, Ref{Int64}), n, X, incx)
        max_idx = max_idx & 0xffffffff # return value might be Int32 if LP64
        if max_idx == 0
            # the `isamax()` implementation saw `N < 0`, ergo it's a 64-bit library
            Int64
        elseif max_idx == 2
            # the `isamax()` implementation saw `N == 3`, ergo it's a 32-bit library
            Int32
        else
            error("failed to detect LP64 or ILP64 interface  in $libopenblas")
        end
    end
end

using OpenBLAS_jll
BlasInt = determine_blas_interface(libopenblas)

deps = quote
    const BlasInt = $(BlasInt)
end

remove_line_numbers(x) = x
function remove_line_numbers(ex::Expr)
    if ex.head == :macrocall
        ex.args[2] = nothing
    else
        ex.args = [remove_line_numbers(arg) for arg in ex.args if !(arg isa LineNumberNode)]
    end
    return ex
end

# only update deps.jl if it has changed.
deps_str = string(remove_line_numbers(deps))

if !isfile("deps.jl") || deps_str != read("deps.jl", String)
    write("deps.jl", deps_str)
end
