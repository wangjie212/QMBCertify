using QMBCertify


# 1d Heisenberg model
supp = Vector{UInt16}[[1;4]]
coe = [3/4]
N = 10 # number of spins
@time opt,cor0,cor1,cor2 = GSB(supp, coe, N, 4, QUIET=true, positivity=0, pso=2, extra=4, correlation=false)

N = 20 # number of spins
@time opt,cor0,cor1,cor2 = GSB(supp, coe, N, 4, QUIET=false, positivity=8, pso=0, extra=9, correlation=false)


# 1d J1-J2 Heisenberg model
N = 10 # number of sites
J2 = 0.3
supp = Vector{UInt16}[[1;4], [1;7]]
coe = [3/4, 3/4*J2]
r = 1
tt = [1;1]
@time opt,cor0,cor1,cor2 = GSB(supp, coe, N, 4, QUIET=true, positivity=0, pso=2, extra=r-1, three_type=tt, correlation=false)


# 2d L×L Heisenberg model
L = 4
supp = [UInt16[1;4]]
coe = [3/2]
@time opt,cor,_,_ = GSB(supp, coe, L, 4, lattice="square", positivity=8, pso=2, extra=0, QUIET=false, correlation=false)


# 2d L×L J1-J2 Heisenberg model
L = 4
supp = Vector{UInt16}[[1;4], [1;7]]
coe = [3/2, 3/2*J2]
@time opt,cor,_,_ = GSB(supp, coe, L, 4, lattice="square", positivity=0, pso=2, extra=0, QUIET=false, correlation=false)


# Ground state computation using DMRG
using ITensors
# 1d Heisenberg model
N = 10 # number of sites
sites = siteinds("S=1/2", N)
os = OpSum()
for j = 1:N-1
  os += "Sx",j,"Sx",j+1
  os += "Sy",j,"Sy",j+1
  os += "Sz",j,"Sz",j+1
end
os += "Sx",1,"Sx",N
os += "Sy",1,"Sy",N
os += "Sz",1,"Sz",N
H = MPO(os, sites)
nsweeps = 5 # number of sweeps is 5
maxdim = [10, 20, 100, 100, 200] # gradually increase states kept
cutoff = [1E-12] # desired truncation error
psi0 = randomMPS(sites, 2)
energy,psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

mo = OpSum()
mo += "Sy",1,"Sz",3,"Sy",4,"Sz",5
val = inner(psi', MPO(mo, sites), psi)
@show val



tsupp = [UInt16[]]
for i = 0:3, j = 0:3, k = 0:3, l = 0:3
  ind = [i,j,k,l]
  if all(x->iseven(sum(ind .== x)), 1:3)
    inx = ind .!= 0
    bi = QMBCertify.reduce4(UInt16.(3*(Vector(1:4)[inx] .- 1) + ind[inx]), 100)
    push!(tsupp, bi)
  end
end
unique!(tsupp)
sort!(tsupp)
y = ceil.(Int, rand(length(tsupp))*1000)
A = zeros(Int, 2^4, 2^4)
Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
for i = 0:3, j = 0:3, k = 0:3, l = 0:3
  ind = [i,j,k,l]
  if all(x->iseven(sum(ind .== x)), 1:3)
    inx = ind .!= 0
    bi = QMBCertify.reduce4(UInt16.(3*(Vector(1:4)[inx] .- 1) + ind[inx]), 100)
    Locb = QMBCertify.bfind(tsupp, length(tsupp), bi)
    A += y[Locb]*real(kron(Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1]))
  end
end
