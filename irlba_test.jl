using SparseArrays
using CSV
using Test

include("scaledata.jl")
include("irlba.jl")

#A = Matrix(CSV.read("tests/A1.csv", header=false))

#U,s,V = irlba(A, 2)

#U2,s2,V2 = irlba(A, 3, U, s, V, init)
#S = svd(A)
#@test s ≈ S.S[1:nu] atol=1e-5

#S_recon = S.U[:,1:nu] * Diagonal(S.S[1:nu]) * S.V[:,1:nu]'
#A_recon = U*Diagonal(s)*V'
#@test S_recon ≈ A_recon atol=2e-6

A = sprand(10000, 500, 0.02)
@time U,s,V = irlba(A, 10)
@time S = svd(Matrix(A))
@time UU,ss,VV = irlba(A'A, 10)
