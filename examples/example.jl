using QMBCertify

# 1d Heisenberg model
supp = [[1;4]]
coe = [3/4]
N = 10 # number of spins
r = 5
@time opt,data = GSB(supp, coe, N, 4, QUIET=false, rdm=8, pso=0, lso=0, extra=r-1)


# 1d J1-J2 Heisenberg model
N = 20 # number of sites
J2 = 1.0
supp = [[1;4], [1;7]]
coe = [3/4; 3/4*J2]
r = 10
tt = [2;2]
@time opt,data = GSB(supp, coe, N, 4, QUIET=false, rdm=0, pso=2, extra=r-1, three_type=tt)
# @time opt,data = GSB(supp, coe, N, 2, QUIET=false, rdm=0, lso=0, pso=2, extra=r-1, three_type=tt, writetofile="D:/Programs/QMBCertify/data/1dL4j1j2_0.3-2.dat-s")


# 2d L×L Heisenberg model
L = 4
supp = [[1;4]]
coe = [3/2]
@time opt,data = GSB(supp, coe, L, 4, lattice="square", rdm=0, pso=0, lso=0, extra=0, QUIET=false)


# 2d L×L J1-J2 Heisenberg model
L = 4
supp = [[1;4], [1;7]]
J2 = 0.3
coe = [3/2, 3/2*J2]
@time opt,data = GSB(supp, coe, L, 4, lattice="square", rdm=0, pso=0, lso=0, extra=0, QUIET=false)


# Ground state computation using DMRG
using ITensors, ITensorMPS
# 1d Heisenberg model
N = 40 # number of sites
# Js = [0.1, 0.2, 0.241167, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0]
# io = open("D:/Programs/QMBCertify/data/H1D2_100_2.txt", "w")
# for J2 in Js
# println("J2 = $J2")
sites = siteinds("S=1/2", N)
os = OpSum()
for j = 1:N-1
  os += "Sx",j,"Sx",j+1
  os += "Sy",j,"Sy",j+1
  os += "Sz",j,"Sz",j+1
end
for j = 1:N-2
  os += J2,"Sx",j,"Sx",j+2
  os += J2,"Sy",j,"Sy",j+2
  os += J2,"Sz",j,"Sz",j+2
end
os += "Sx",1,"Sx",N
os += "Sy",1,"Sy",N
os += "Sz",1,"Sz",N
os += J2,"Sx",1,"Sx",N-1
os += J2,"Sy",1,"Sy",N-1
os += J2,"Sz",1,"Sz",N-1
os += J2,"Sx",2,"Sx",N
os += J2,"Sy",2,"Sy",N
os += J2,"Sz",2,"Sz",N
H = MPO(os, sites)
nsweeps = 7 # number of sweeps
maxdim = [10, 20, 100, 200, 400, 800, 1600] # gradually increase states kept
cutoff = [1E-12] # desired truncation error
psi0 = randomMPS(sites, 2)
energy,psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
energy = energy/N
mo = OpSum()
mo += "Sx",1,"Sx",2
c1 = real(inner(psi', MPO(mo, sites), psi))
mo = OpSum()
mo += "Sx",1,"Sx",3
c2 = real(inner(psi', MPO(mo, sites), psi))
# write(io, "J2 = $J2, energy = $energy, C1 = $c1, C2 = $c2\n")
# write(io, "psi = $psi\n")
# write(io, "------------------------------------\n")
# end
# close(io)

a = 0
for i = 1:N-1
  mo = OpSum()
  mo += "Sx",1,"Sx",2,"Sx",i,"Sx",i+1
  a += 3*(-1)^(1+i)*inner(psi', MPO(mo, sites), psi)
  mo = OpSum()
  mo += "Sx",1,"Sx",2,"Sy",i,"Sy",i+1
  a += 6*(-1)^(1+i)*inner(psi', MPO(mo, sites), psi)
end
mo = OpSum()
mo += "Sx",1,"Sx",3
a += -3/4*inner(psi', MPO(mo, sites), psi)
