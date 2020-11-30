function show(io::IO, em::LinearEmbedding; dims::AbstractVector{<:Integer}=1:5, nfeatures::Integer=20)
    loadings = em.basis
    for dim in dims
        x = loadings.array[:,dim]
        num = round(Int64, nfeatures / 2)
        positives = partialsortperm(x, 1:num, rev=true)
        negatives = partialsortperm(x, 1:num, rev=false)
        println(io, names(loadings,2)[dim])
        println(io, "Positive: ", names(loadings,1)[positives])
        println(io, "Negative: ", names(loadings,1)[negatives])
    end
end

function show(io::IO, C::CenteredMatrix)
    print(io, "CenteredMatrix(A=$(C.A), mu=$(C.mu))")
end

function show(io::IO, ::MIME"text/plain", C::CenteredMatrix)
    println(io, "CenteredMatrix:")
    ioc = IOContext(io, :compact=>true, :limit=>true)
    println(ioc, "  A = ", C.A)
    print(ioc, "  mu = ", C.mu)
end

