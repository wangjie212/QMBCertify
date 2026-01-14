using QMBCertify

J2   = 0.2                      
supp = [[1;4], [1;7]]             
coe  = [3/4; 3/4*J2]               
tt   = [1;1]                    

Ns = [4,6,8,10,12,14,16,18, 20 ,22, 24, 26, 28, 30]

oldbounds  = Float64[]            
newbounds  = Float64[]
shifts     = Float64[]
mineigs_list = []
dims   = []
times = []

for N in Ns
    println("=== N = $N ===")

    time = @elapsed begin

        opt, data = GSB(supp, coe, N, 2;
            QUIET=true, rdm=0, lol=N, extra=1, pso=0, lso=0, three_type=tt, Gram=true)

        result = certify_qmb(data, N, coe[1], opt; tol_gram=1e-15, tol_dft=1e-12, snn=true, J2=J2)

    end

    println("Total time for N=$N: $time seconds")
    push!(times, time)
    push!(oldbounds,  Float64(result.oldbound))
    push!(newbounds,  Float64(result.newbound))
    push!(shifts,     Float64(result.shift))
    push!(mineigs_list, result.mineigs)
    push!(dims,   result.dims)
end