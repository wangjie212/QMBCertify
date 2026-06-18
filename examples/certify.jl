using QMBCertify

# 1d Heisenberg model
supp = [[1;4]]
coe = [3/4]
N = 4 # number of spins
opt,data = GSB(supp, coe, N, 3, QUIET=false, rdm=0, pso=2, lso=0, Gram=true)

RHS = zeros(length(data.tsupp))
s = Int.(length.(data.basis)/N)
P = [cos(2*pi*(i-1)*(j-1)/N) + sin(2*pi*(i-1)*(j-1)/N)*im for i = 1:N, j = 1:N]
RHS[1] += data.GramMat[1][1][1,1]
for k = 1:s[1]
    w,c = reduce!(data.basis[1][N*(k-1)+1], L=N)
    if c != 0
        Locb = QMBCertify.bfind(data.tsupp, w)
        RHS[Locb] += 2*sqrt(N)*data.GramMat[1][1][1,k+1]
    end
end
for l = 1:2, i = 1:Int(N/2)+1, j = 1:s[l], k = j:s[l], v = 1:N
    w,c = reduce!([data.basis[l][N*(j-1)+1]; data.basis[l][N*(k-1)+v]], L=N, realify=true)
    if c != 0
        Locb = QMBCertify.bfind(data.tsupp, w)
        if j == k
            if l == 1 && i == 1
                RHS[Locb] += c*data.GramMat[1][1][j+1,k+1]
            else
                RHS[Locb] += c*real(data.GramMat[l][i][j,k]*P[i,v])
            end
        else
            if l == 1 && i == 1
                RHS[Locb] += 2*c*data.GramMat[1][1][j+1,k+1]
            else
                RHS[Locb] += 2*c*real(data.GramMat[l][i][j,k]*P[i,v])
            end
        end
    end
end
println(RHS) # The coefficient vector of RHS w.r.t. the monomial vector stored in data.tsupp



RHS = zeros(length(data.tsupp))
s = Int.(length.(data.basis)/N)
P = [1/sqrt(N)*(cos(2*pi*(i-1)*(j-1)/N) - sin(2*pi*(i-1)*(j-1)/N)*im) for i = 1:N, j = 1:N]
RHS[1] += data.GramMat[1][1][1,1]
for k = 1:s[1]
    w,c = reduce!(data.basis[1][N*(k-1)+1], L=N)
    if c != 0
        Locb = QMBCertify.bfind(data.tsupp, w)
        RHS[Locb] += 2*sqrt(N)*data.GramMat[1][1][1,k+1]
    end
end
for l = 1:2, i = 1:Int(N/2)+1, j = 1:s[l], k = j:s[l], u = 1:N, v = 1:N
    w,c = reduce!([data.basis[l][N*(j-1)+u]; data.basis[l][N*(k-1)+v]], L=N, realify=true)
    if c != 0
        Locb = QMBCertify.bfind(data.tsupp, w)
        if j == k
            if l == 1 && i == 1
                RHS[Locb] += c*data.GramMat[1][1][j+1,k+1]*P[1,u]*conj(P[1,v])
            else
                RHS[Locb] += c*real(data.GramMat[l][i][j,k]*P[i,u]*conj(P[i,v]))
            end
        else
            if l == 1 && i == 1
                RHS[Locb] += 2*c*data.GramMat[1][1][j+1,k+1]*P[1,u]*conj(P[1,v])
            else
                RHS[Locb] += 2*c*real(data.GramMat[l][i][j,k]*P[i,u]*conj(P[i,v]))
            end
        end
    end
end
println(RHS) # The coefficient vector of RHS w.r.t. the monomial vector stored in data.tsupp



N = 10
J2 = 1
supp = [[1;4], [1;7]]
coe = [3/4; 3/4*J2]
opt,data = GSB(supp, coe, N, 4, QUIET=false, rdm=8, pso=3, lso=0, Gram=true)


RHS = zeros(length(data.tsupp))
s = Int.(length.(data.basis)/N)
P = [cos(2*pi*(i-1)*(j-1)/N) + sin(2*pi*(i-1)*(j-1)/N)*im for i = 1:N, j = 1:N]
RHS[1] += data.GramMat[1][1][1,1]
for k = 1:s[1]
    w,c = reduce!(data.basis[1][N*(k-1)+1], L=N)
    if c != 0
        Locb = QMBCertify.bfind(data.tsupp, w)
        RHS[Locb] += 2*sqrt(N)*data.GramMat[1][1][1,k+1]
    end
end
for l = 1:2, i = 1:Int(N/2)+1, j = 1:s[l], k = j:s[l], v = 1:N
    w,c = reduce!([data.basis[l][N*(j-1)+1]; data.basis[l][N*(k-1)+v]], L=N, realify=true)
    if c != 0
        Locb = QMBCertify.bfind(data.tsupp, w)
        if j == k
            if l == 1 && i == 1
                RHS[Locb] += c*data.GramMat[1][1][j+1,k+1]
            else
                RHS[Locb] += c*real(data.GramMat[l][i][j,k]*P[i,v])
            end
        else
            if l == 1 && i == 1
                RHS[Locb] += 2*c*data.GramMat[1][1][j+1,k+1]
            else
                RHS[Locb] += 2*c*real(data.GramMat[l][i][j,k]*P[i,v])
            end
        end
    end
end
s = Int.(length.(data.sbasis)/N)
h = [ceil(Int, item[2]/3) - 1 for item in supp]
for l = 1:2, i = 1:Int(N/2)+1, j = 1:s[l], k = j:s[l], v = 1:N
    ww,cc = QMBCertify.PSDstate_entry(data.sbasis[l][N*(j-1)+1], data.sbasis[l][N*(k-1)+v], h, coe, N)
    for (t, w) in enumerate(ww)
        Locb = QMBCertify.bfind(data.tsupp, w)
        if j == k
            RHS[Locb] += cc[t]*real(data.sGramMat[l][i][j,k]*P[i,v])
        else
            RHS[Locb] += 2*cc[t]*real(data.sGramMat[l][i][j,k]*P[i,v])
        end
    end
end
println(RHS)
