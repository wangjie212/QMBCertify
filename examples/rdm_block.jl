using QMBCertify
using Graphs

d = 4
tsupp = [UInt16[]]
for i = 0:3, j = 0:3, k = 0:3, l = 0:3
  ind = [i,j,k,l]
  if all(x->iseven(sum(ind .== x)), 1:3)
    inx = ind .!= 0
    bi = QMBCertify.reduce4(UInt16.(3*(Vector(1:d)[inx] .- 1) + ind[inx]), 100) # reduce a monomial to the normal form
    push!(tsupp, bi)
  end
end
unique!(tsupp)
sort!(tsupp)
y = ceil.(Int, rand(length(tsupp))*10^6) # use random integers to represent coefficients
A = zeros(Int, 2^d, 2^d)
Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
for i = 0:3, j = 0:3, k = 0:3, l = 0:3
  ind = [i,j,k,l]
  if all(x->iseven(sum(ind .== x)), 1:3)
    inx = ind .!= 0
    bi = QMBCertify.reduce4(UInt16.(3*(Vector(1:d)[inx] .- 1) + ind[inx]), 100)
    Locb = QMBCertify.bfind(tsupp, length(tsupp), bi)
    A += y[Locb]*real(kron(Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1]))
  end
end

G = SimpleGraph(2^d)
for i = 1:2^d, j = i+1:2^d
    if abs(A[i,j]) > 1e-8
        add_edge!(G, i, j)
    end
end
blocks = connected_components(G)