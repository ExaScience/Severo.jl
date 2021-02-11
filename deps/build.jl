# copyright imec - evaluation license - not for distribution

using Libdl
using LinearAlgebra

function determine_flags()
    # XXX assume BLAS and LAPACK have the same vendor
    LDFLAGS = String[]
    CFLAGS = String[]

    if BLAS.BlasInt === Int64
        push!(CFLAGS, "-DBLAS64")
    end

    if BLAS.determine_vendor() === :openblas64
        push!(CFLAGS, "-DHAVE_OPENBLAS64")
    else
        push!(CFLAGS, "-DHAVE_F77_UNDERSCORE")
    end

    for libname in [BLAS.libblas, BLAS.liblapack]
        lib = nothing
        try
            lib = Libdl.dlopen(libname)
            libpath = dlpath(lib)
            push!(LDFLAGS, libpath)
        finally
            Libdl.dlclose(lib)
        end
    end

    join(unique(CFLAGS), " "), join(unique(LDFLAGS), " ")
end

function build()
    cflags, ldflags = determine_flags()

    p = pwd()
    cd(dirname(@__FILE__))
    if Sys.iswindows()
        try
            run(`mingw32-make LAPACK_CFLAGS="$cflags" LAPACK_LDFLAGS="$ldflags"`)
        catch
            error("failed to build using mingw32-make")
        end
    elseif Sys.isbsd() && !Sys.isapple()  # e.g. FreeBSD
        run(`gmake LAPACK_CFLAGS="$cflags" LAPACK_LDFLAGS="$ldflags"`)
    else
        run(`make LAPACK_CFLAGS="$cflags" LAPACK_LDFLAGS="$ldflags"`)
    end
    cd(p)
end

build()

const libcell = find_library(["libcell"], [dirname(@__FILE__)])
if libcell == ""
		error("libcell could not be found, compilation must have failed")
end

deps = quote
    const libcell = $libcell
    const BlasInt = $(BLAS.BlasInt)
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
# allows users to call Pkg.build("MPI") without triggering another round of precompilation
deps_str = string(remove_line_numbers(deps))

if !isfile("deps.jl") || deps_str != read("deps.jl", String)
    write("deps.jl", deps_str)
end
