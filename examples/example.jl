using QMBCertify

# 1d Heisenberg model
supp = Vector{UInt16}[[1;4]]
coe = [3/4]
N = 8 # number of sites
d = 4 # relaxation order
@time begin
opt,cor0,cor1,cor2 = GSB(supp, coe, N, d, QUIET=true, posepsd=false, extra=0, correlation=true)
end

# 1d J1-J2 Heisenberg model
N = 10 # number of sites
J2 = 0.3
supp = Vector{UInt16}[[1;4], [1;7]]
coe = [3/4;3/4*J2]
d,r = 4,5
tt = [1;1]
@time begin
opt,cor0,cor1,cor2 = GSB(supp, coe, N, d, QUIET=true, posepsd=true, extra=r-1, three_type=tt, correlation=false)
end

# 2d L×L Heisenberg model
L = 4 
supp = [UInt16[1;4]]
coe = [3/2]
@time begin
opt,cor,_,_ = GSB(supp, coe, L, 4, lattice="square", extra=true, QUIET=true, correlation=false)
end

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
