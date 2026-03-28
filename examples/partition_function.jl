using QMBCertify

# 1d Heisenberg model
supp = [[1;4]]
coe = [3/4]
N = 10 # number of spins
beta = 0.1
@time opt = PFB(supp, coe, beta, N, 3, QUIET=false)
println(2^N*opt)



using LinearAlgebra

Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
H = zeros(Int, 2^N, 2^N)
for i = 1:N-1, j = 2:4
    ind = ones(Int, N)
    ind[i] = ind[i+1] = j
    H += real(kron(Pauli[ind]...))
end
for j = 2:4
    ind = ones(Int, N)
    ind[1] = ind[N] = j
    H += real(kron(Pauli[ind]...))
end
v = eigvals(H/4)
println(sum(exp.(-beta*v)))