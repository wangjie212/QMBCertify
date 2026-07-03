using QMBCertify

# 1d Heisenberg model
supp = [[1;4]]
coe = [3/4]
lb = [-0.4515446, -0.4452196, -0.4440668, -0.4436649, -0.4434798, -0.4433804, -0.4432808, -0.4432378]


## Bound the partition function and the free energy
N = 50 # number of spins
beta = 1
d = 4

@time LB = PFB(supp, coe, beta, N, d, obj="min", QUIET=false)
a = N*log10(2) + log10(LB)
fa = floor(a)
println([round(10^(a-fa), digits=3), fa])
free_energy = -1/beta*(N*log(2) + log(LB))
println(round(free_energy, digits=2))

@time UB = PFB(supp, coe, beta, N, 3, obj="max", QUIET=false)
a = N*log10(2) + log10(UB)
fa = floor(a)
println([round(10^(a-fa), digits=3), fa])
free_energy = -1/beta*(N*log(2) + log(UB))
println(round(free_energy, digits=2))


## Bound the energy expectation in the thermal state
N = 10 # number of spins
beta = 0.1
d = 3

@time LB = PFB(supp, coe, beta, N, d, obj="min", QUIET=true)
@time UB = PFB(supp, coe, beta, N, d, obj="max", QUIET=true)
@time HLB = PFB(supp, coe, beta, N, d, bound_energy=true, obj="min", QUIET=true)
@time HUB = PFB(supp, coe, beta, N, d, bound_energy=true, obj="max", QUIET=true)
println("$(round(HLB, digits=4)) $(round(HUB, digits=4)) $(round(HLB/UB, digits=4)) $(round(HUB/LB, digits=4))")


## Compute the exact partition function
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



for i = 1:10
beta = i
println([i, exp(beta*0.75)])
end



using Plots

N = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
F = [[-69.35, -69.61], [-138.4, -140.6], [-206.4, -212.9], [-275.0, -282.2], [-343.3, -353.3],
[-412.1, -421.6], [-480.4, -491.7], [-549.3, -561.5], [-618.2, -632.1], [-686.1, -701.1]]
p = plot(N, [F[i][1] for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 0.1"], xlabel="number of sites", ylabel="bounds on free energies")
plot!(p, N, [F[i][2] for i = 1:10], dpi=600, shape=:star, label=["lower bound when β = 0.1"])
savefig("d:/Programs/QMBCertify/data/partition_function/free_energy_beta_0.1.png")


N = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
F = [[-14.01, -14.53], [-27.76, -28.92], [-41.59, -43.08], [-55.43, -57.11], [-69.28, -71.12], 
[-83.14, -85.13], [-96.99, -98.98], [-110.8, -112.9], [-124.7, -126.8], [-138.6, -140.7]]
p = plot(N, [F[i][1] for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 0.5"], xlabel="number of sites", ylabel="bounds on free energies")
plot!(p, N, [F[i][2] for i = 1:10], dpi=600, shape=:star, label=["lower bound when β = 0.5"])
savefig("d:/Programs/QMBCertify/data/partition_function/free_energy_beta_0.5.png")


N = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
F = [[-6.980, -8.690], [-13.85, -16.02], [-20.76, -22.90], [-27.68, -29.88], [-34.60, -36.97], 
[-41.51, -43.76], [-48.44, -50.76], [-55.35, -57.84], [-62.29, -64.46], [-69.20, -71.40]]
p = plot(N, [F[i][1] for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 1.0"], xlabel="number of sites", ylabel="bounds on free energies")
plot!(p, N, [F[i][2] for i = 1:10], dpi=600, shape=:star, label=["lower bound when β = 1.0"])
savefig("d:/Programs/QMBCertify/data/partition_function/free_energy_beta_1.0.png")


N = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
P = [[1.028e3, 1.054e3], [1.022e6, 1.273e6], [9.158e8, 1.754e9], [8.774e11, 1.797e12], [8.104e14, 2.203e15], 
[7.854e17, 2.042e18], [7.336e20, 2.263e21], [7.153e23, 2.423e24], [7.035e26, 2.819e27], [6.284e29, 2.796e30]]
p = plot(N, [log(P[i][2]) for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 0.1"], xlabel="number of sites", ylabel="log(bound) on partition functions")
plot!(p, N, [log(P[i][1]) for i = 1:10], dpi=600, shape=:star, label=["lower bound when β = 0.1"])
savefig("d:/Programs/QMBCertify/data/partition_function/partition_function_beta_0.1.png")


N = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
P = [[1.103e3, 1.426e3], [1.067e6, 1.901e6], [1.075e9, 2.264e9], [1.090e12, 2.514e12], [1.108e15, 2.783e15],
[1.130e18, 3.054e18], [1.152e21, 3.106e21], [1.173e24, 3.333e24], [1.198e27, 3.401e27], [1.224e30, 3.525e30]]
p = plot(N, [log(P[i][2]) for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 0.5"], xlabel="number of sites", ylabel="log(bound) on partition functions")
plot!(p, N, [log(P[i][1]) for i = 1:10], dpi=600, shape=:star, label=["lower bound when β = 0.5"])
savefig("d:/Programs/QMBCertify/data/partition_function/partition_function_beta_0.5.png")


N = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
P = [[1.078e3, 5.915e3], [1.034e6, 9.110e6], [1.040e9, 8.833e9], [1.045e12, 9.494e12], [1.059e15, 1.136e16],
[1.069e18, 1.013e19], [1.085e21, 1.109e22], [1.096e24, 1.322e25], [1.123e27, 9.836e27], [1.135e30, 1.025e31]]
p = plot(N, [log(P[i][2]) for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 1.0"], xlabel="number of sites", ylabel="log(bound) on partition functions")
plot!(p, N, [log(P[i][1]) for i = 1:10], dpi=600, shape=:star, label=["lower bound when β = 1.0"])
savefig("d:/Programs/QMBCertify/data/partition_function/partition_function_beta_1.0.png")


beta = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
P = [[1.028e3, 1.054e3], [1.069e3, 1.106e3], [1.091e3, 1.132e3], [1.101e3, 1.273e3], [1.117e3, 1.377e3],   
[1.099e3, 1.615e3], [1.095e3, 1.803e3], [1.088e3, 4.845e3], [1.081e3, 5.500e3], [1.078e3, 5.915e3]]
p = plot(beta, [P[i][2] for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when N = 10"], xlabel="β", ylabel="bounds on partition functions")
plot!(p, beta, [P[i][1] for i = 1:10], dpi=600, shape=:star, label=["lower bound when N = 10"])
savefig("d:/Programs/QMBCertify/data/partition_function/partition_function_N_10.png")


beta = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
F = [[-69.35, -69.61], [-34.87, -35.04], [-23.32, -23.44], [-17.51, -17.87], [-14.04, -14.45],
[-11.67, -12.31], [-10.00, -10.71], [-8.740, -10.61], [-7.762, -9.569], [-6.983, -8.685]]
p = plot(beta, [F[i][1] for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when N = 10"], legend=:bottomright, xlabel="β", ylabel="bounds on free energies")
plot!(p, beta, [F[i][2] for i = 1:10], dpi=600, shape=:star, label=["lower bound when N = 10"])
savefig("d:/Programs/QMBCertify/data/partition_function/free_energy_N_10.png")
  

beta = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
P = [[8.104e14, 2.203e15], [1.132e15, 3.705e15], [1.125e15, 1.947e15],  [1.118e15, 2.423e15], [1.108e15, 2.783e15], 
[1.100e15, 3.039e15], [1.089e15, 3.323e15], [1.078e15, 9.768e15], [1.070e15, 1.029e16], [1.059e15, 1.136e16]]
p = plot(beta, [P[i][2] for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when N = 50"], xlabel="β", ylabel="bounds on partition functions")
plot!(p, beta, [P[i][1] for i = 1:10], dpi=600, shape=:star, label=["lower bound when N = 50"])
savefig("d:/Programs/QMBCertify/data/partition_function/partition_function_N_50.png")


beta = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
F = [[-343.3, -353.3], [-173.3, -179.2], [-115.5, -117.3], [-86.63, -88.56], [-69.28, -71.12], 
[-57.72, -59.42], [-49.46, -51.06], [-43.27, -46.02], [-38.45, -40.97], [-34.60, -36.97]]
p = plot(beta, [F[i][1] for i = 1:10], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when N = 50"], legend=:bottomright, xlabel="β", ylabel="bounds on free energies")
plot!(p, beta, [F[i][2] for i = 1:10], dpi=600, shape=:star, label=["lower bound when N = 50"])
savefig("d:/Programs/QMBCertify/data/partition_function/free_energy_N_50.png")


N = [10, 20, 30, 40, 50]
E = [[-0.3978, -0.1046], [-3.0665, 0.6132], [-7.484, 3.156], [-6.979, 4.702], [-10.16, 6.545]]
p = plot(N, [E[i][2] for i = 1:5], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 0.1"], xlabel="number of sites", ylabel="bounds on finite sized entropies")
plot!(p, N, [E[i][1] for i = 1:5], dpi=600, shape=:star, label=["lower bound when β = 0.1"])
savefig("d:/Programs/QMBCertify/data/partition_function/finite_sized_entropy_0.1.png")


N = [10, 20, 30, 40, 50]
E = [[-1.592, -0.3758], [-2.100, -0.1561], [-1.617, -0.1957], [-0.8623, -0.8068], [-0.8822, -0.8815]]
p = plot(N, [E[i][2] for i = 1:5], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 0.5"], legend=:bottomright, xlabel="number of sites", ylabel="bounds on finite sized entropies")
plot!(p, N, [E[i][1] for i = 1:5], dpi=600, shape=:star, label=["lower bound when β = 0.5"])
savefig("d:/Programs/QMBCertify/data/partition_function/finite_sized_entropy_0.5.png")


N = [10, 20, 30, 40, 50]
E = [[-1.3569, -0.1847], [-1.571, -0.0750], [-1.8909, -0.0437], [-2.1301, -0.0269], [-2.4454, -0.0195]]
p = plot(N, [E[i][2] for i = 1:5], dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 1.0"], xlabel="number of sites", ylabel="bounds on finite sized entropies")
plot!(p, N, [E[i][1] for i = 1:5], dpi=600, shape=:star, label=["lower bound when β = 1.0"])
savefig("d:/Programs/QMBCertify/data/partition_function/finite_sized_entropy_1.0.png")


P1 = [[1.103e3, 1.426e3], [1.067e6, 1.901e6], [1.075e9, 2.264e9], [1.090e12, 2.514e12], [1.108e15, 2.783e15], 
[1.130e18, 3.054e18], [1.152e21, 3.106e21], [1.173e24, 3.333e24], [1.198e27, 3.401e27], [1.224e30, 3.525e30]]
P2 = [[1.078e3, 5.915e3], [1.034e6, 9.110e6], [1.040e9, 8.833e9], [1.045e12, 9.494e12], [1.059e15, 1.136e16],
[1.069e18, 1.013e19], [1.085e21, 1.109e22], [1.096e24, 1.322e25], [1.123e27, 9.836e27], [1.135e30, 1.025e31]]
R_ub = [P2[i][2]/P1[i][1]^2 for i = 1:10]
R_lb = [P2[i][1]/P1[i][2]^2 for i = 1:10]

N = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
p = plot(N, log.(R_ub), dpi=600, shape=:circle, legendfontsize=10, label=["upper bound when β = 0.5"], xlabel="number of sites", ylabel="log(bound) on purity of thermal states")
plot!(p, N, log.(R_lb), dpi=600, shape=:star, label=["lower bound when β = 0.5"])
savefig("d:/Programs/QMBCertify/data/partition_function/purity_thermal_state_beta_0.5.png")


N = 10
beta = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
UB = [1034, 1064, 1118, 1200, 1315, 1472, 1683, 1965, 2342, 2847]
T = [log(2^N + 2) - log(UB[i]) - beta[i]*0.25*N for i = 1:10]

p = plot(beta, T, dpi=600, shape=:circle, legendfontsize=10, label=["test value when N = 10"], xlabel="β", ylabel="test on absolutely separable")
savefig("d:/Programs/QMBCertify/data/partition_function/absolutely_separable_N_10.png")


N = 50
beta = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
UB = [2.203e15, 3.705e15, 1.947e15, 2.423e15, 2.783e15, 3.039e15, 3.323e15, 9.768e15, 1.029e16, 1.136e16]
T = [log(2^N + 2) - log(UB[i]) - beta[i]*0.249815*N for i = 1:10]

p = plot(beta, T, dpi=600, shape=:circle, legendfontsize=10, label=["test value when N = 50"], xlabel="β", ylabel="test on absolutely separable")
savefig("d:/Programs/QMBCertify/data/partition_function/absolutely_separable_N_50.png")
