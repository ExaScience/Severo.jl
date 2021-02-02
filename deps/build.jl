# copyright imec - evaluation license - not for distribution

using Libdl

function build()
    p = pwd()
    cd(dirname(@__FILE__))
    if Sys.iswindows()
        try
            run(`mingw32-make`)
        catch
            error("failed to build using mingw32-make")
        end
    elseif Sys.isbsd() && !Sys.isapple()  # e.g. FreeBSD
        run(`gmake`)
    else
        run(`make`)
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
