using Pkg
cd("/home/jwang/QMBCertify")
Pkg.activate(".")
using QMBCertify

supp = Vector{UInt16}[[1;4]]
coe = [3/4]
GSB(supp, coe, 10, 4, QUIET=true, positivity=false, pso=3, extra=4, correlation=true)

supp = [UInt16[1;4]]
coe = [3/2]
io = open("/home/jwang/H2D_L6_3.txt", "w")
L = 16
d = 2
write(io, "L = $L, d = $d, pso = 2, extra = 0\n")
t = @elapsed begin
opt,cor,_,_ = GSB(supp, coe, L, d, lattice="square", positivity=10, pso=2, extra=0, QUIET=false, correlation=true)
end
write(io, "E/N = $opt\n")
write(io, "time = $t\n")
write(io, "correlation = $cor\n")
write(io, "------------------------------------\n")
close(io)

GSB(supp, coe, 100, 4, lattice="chain", positivity=10, pso=3, extra=0, QUIET=false, correlation=true)
