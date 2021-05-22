function GSE1(supp::Vector{Vector{UInt16}}, coe::Vector{Float64}, L::Int, d::Int; energy=[], QUIET=false, lattice="chain",
    solver="Mosek", posepsd=false, extra=0, three_type=[1;1], J2=0, totalspin=false, sector=0, correlation=false)
    basis = Vector{Vector{Vector{UInt16}}}(undef, 4)
    tsupp = Vector{UInt16}[]
    for i = 0:3
        basis[i+1] = split_basis(L, d, i, lattice=lattice, extra=extra, three_type=three_type)
        for j = 1:length(basis[i+1]), k = j:length(basis[i+1])
            @inbounds bi = [basis[i+1][j]; basis[i+1][k]]
            bi,coef = reduce!(bi, L=L, lattice=lattice)
            if coef != 0
                push!(tsupp, bi)
            end
        end
    end
    if posepsd == true && lattice == "chain"
        for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3
            sites = [[1;3;4;5;6;7], [1;2;4;5;6;7], [1;2;3;5;6;7]]
            for r = 1:3
                mon = mono([i,j,k,l,s,t], sites=sites[r])
                if !iszero(mon)
                    push!(tsupp, reduce4(mon, L))
                end
            end
        end
        # for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3
        #     mon = mono([i,j,k,l,s,t,u,v])
        #     if !iszero(mon)
        #         push!(tsupp, reduce4(mon, L))
        #     end
        # end
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp = length(tsupp)
    if solver == "COSMO"
        model = Model(optimizer_with_attributes(COSMO.Optimizer))
        set_optimizer_attributes(model, "eps_abs" => 1e-4, "eps_rel" => 1e-4, "max_iter" => 100000)
    else
        model = Model(optimizer_with_attributes(Mosek.Optimizer))
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    mvar = @variable(model, [1:ltsupp])
    for i = 0:3
        if i==0
            k=Int((length(basis[1])-1)/L)
        else
            k=Int(length(basis[i+1])/L)
        end
        reig=Vector{Matrix{GenericAffExpr{Float64,VariableRef}}}(undef, L)
        ieig=Vector{Matrix{GenericAffExpr{Float64,VariableRef}}}(undef, L)
        for j=1:L
            if i==0&&j==1
                reig[1]=[AffExpr(0) for j1=1:k+1, j2=1:k+1]
                ieig[1]=[AffExpr(0) for j1=1:k+1, j2=1:k+1]
                reig[1][1,1]=AffExpr(1)
                for l=1:k
                    bi,coef=reduce!(basis[1][L*(l-1)+2], L=L, lattice=lattice)
                    if coef!=0
                        Locb=bfind(tsupp, ltsupp, bi)
                        reig[1][1,l+1]=sqrt(L)*mvar[Locb]
                    end
                end
            else
                reig[j]=[AffExpr(0) for j1=1:k, j2=1:k]
                ieig[j]=[AffExpr(0) for j1=1:k, j2=1:k]
            end
        end
        for j1=1:k, j2=j1:k
            if j1==j2
                coef=Vector{ComplexF64}(undef, Int(L/2))
                Locb=Vector{Int}(undef, Int(L/2))
                for r=1:Int(L/2)
                    if i==0
                        bi=[basis[1][L*(j1-1)+2]; basis[1][L*(j2-1)+r+2]]
                    else
                        bi=[basis[i+1][L*(j1-1)+1]; basis[i+1][L*(j2-1)+r+1]]
                    end
                    bi,coef[r]=reduce!(bi, L=L, lattice=lattice)
                    if coef[r]!=0
                        Locb[r]=bfind(tsupp, ltsupp, bi)
                    else
                        Locb[r]=0
                    end
                end
                for l=0:L-1
                    r1,r2=j1,j2
                    if i==0&&l==0
                        r1,r2=j1+1,j2+1
                    end
                    reig[l+1][r1, r2]=1
                    if coef[end]!=0
                        @inbounds add_to_expression!(reig[l+1][r1, r2], (-1)^l, mvar[Locb[end]])
                    end
                    for r=1:Int(L/2)-1
                        if coef[r]^2==1&&abs(cos(2*pi*r*l/L))>=1e-8
                            @inbounds add_to_expression!(reig[l+1][r1, r2], 2*coef[r]*cos(2*pi*r*l/L), mvar[Locb[r]])
                        elseif coef[r]==im&&abs(sin(2*pi*r*l/L))>=1e-8
                            @inbounds add_to_expression!(reig[l+1][r1, r2], -2*sin(2*pi*r*l/L), mvar[Locb[r]])
                        elseif coef[r]==-im&&abs(sin(2*pi*r*l/L))>=1e-8
                            @inbounds add_to_expression!(reig[l+1][r1, r2], 2*sin(2*pi*r*l/L), mvar[Locb[r]])
                        end
                    end
                end
            else
                coef=Vector{ComplexF64}(undef, L)
                Locb=Vector{Int}(undef, L)
                for r=0:L-1
                    if i==0
                        bi=[basis[1][L*(j1-1)+2]; basis[1][L*(j2-1)+r+2]]
                    else
                        bi=[basis[i+1][L*(j1-1)+1]; basis[i+1][L*(j2-1)+r+1]]
                    end
                    bi,coef[r+1]=reduce!(bi, L=L, lattice=lattice)
                    if coef[r+1]!=0
                        Locb[r+1]=bfind(tsupp, ltsupp, bi)
                    else
                        Locb[r+1]=0
                    end
                end
                for l=0:L-1
                    r1,r2=j1,j2
                    if i==0&&l==0
                        r1,r2=j1+1,j2+1
                    end
                    for r=0:L-1
                        if coef[r+1]^2==1
                            if abs(cos(2*pi*r*l/L))>=1e-8
                                @inbounds add_to_expression!(reig[l+1][r1, r2], coef[r+1]*cos(2*pi*r*l/L), mvar[Locb[r+1]])
                            end
                            if abs(sin(2*pi*r*l/L))>=1e-8
                                @inbounds add_to_expression!(ieig[l+1][r1, r2], coef[r+1]*sin(2*pi*r*l/L), mvar[Locb[r+1]])
                            end
                        elseif coef[r+1]==im
                            if abs(sin(2*pi*r*l/L))>=1e-8
                                @inbounds add_to_expression!(reig[l+1][r1, r2], -sin(2*pi*r*l/L), mvar[Locb[r+1]])
                            end
                            if abs(cos(2*pi*r*l/L))>=1e-8
                                @inbounds add_to_expression!(ieig[l+1][r1, r2], cos(2*pi*r*l/L), mvar[Locb[r+1]])
                            end
                        elseif coef[r+1]==-im
                            if abs(sin(2*pi*r*l/L))>=1e-8
                                @inbounds add_to_expression!(reig[l+1][r1, r2], sin(2*pi*r*l/L), mvar[Locb[r+1]])
                            end
                            if abs(cos(2*pi*r*l/L))>=1e-8
                                @inbounds add_to_expression!(ieig[l+1][r1, r2], -cos(2*pi*r*l/L), mvar[Locb[r+1]])
                            end
                        end
                    end
                end
            end
        end
        for l=1:L
            r=k
            if i==0&&l==1
                r=k+1
            end
            pos=[AffExpr(0) for j1=1:2*r, j2=1:2*r]
            pos[1:r,1:r]=reig[l]
            pos[r+1:2r,r+1:2r]=reig[l]
            pos[1:r,r+1:2r]=ieig[l]'-ieig[l]
            @constraint(model, Symmetric(pos) in PSDCone())
        end
    end
    if totalspin == true
        J1 = AffExpr(0)
        for i=1:L, j=1:L
            temp=UInt16[3*(i-1)+1;3*(j-1)+1]
            bi=reduce!(temp, L=L, lattice=lattice)[1]
            Locb=bfind(tsupp,ltsupp,bi)
            @inbounds add_to_expression!(J1, 3, mvar[Locb])
        end
        @constraint(model, J1==4*sector*(sector+1))
        # J2=AffExpr(0)
        # for j1=1:L, j2=1:L, k1=1:L, k2=1:L
        #     temp=UInt16[3*(j1-1)+1;3*(k1-1)+1;3*(j2-1)+1;3*(k2-1)+1]
        #     bi=reduce!(temp, L=L, lattice=lattice)[1]
        #     Locb=bfind(tsupp,ltsupp,bi)
        #     J2+=3*mvar[Locb]
        #     temp=UInt16[3*(j1-1)+1;3*(k1-1)+1;3*(j2-1)+2;3*(k2-1)+2]
        #     bi,coef=reduce!(temp, L=L, lattice=lattice)
        #     Locb=bfind(tsupp,ltsupp,bi)
        #     if coef^2==1
        #         J2+=6*coef*mvar[Locb]
        #     end
        # end
        # @constraint(model, J2==16*sector^2*(sector+1)^2)
    end
    obj = AffExpr(0)
    for i = 1:length(supp)
        Locb = bfind(tsupp,ltsupp,supp[i])
        @inbounds add_to_expression!(obj, coe[i], mvar[Locb])
    end
    @constraint(model, mvar[1]==1)
    if posepsd == true && lattice == "chain"
        posepsd6!(model, mvar, tsupp, L)
        posepsd6!(model, mvar, tsupp, L, sites=[1;3;4;5;6;7])
        posepsd6!(model, mvar, tsupp, L, sites=[1;2;4;5;6;7])
        posepsd6!(model, mvar, tsupp, L, sites=[1;2;3;5;6;7])
        # posepsd8!(model, mvar, tsupp, L)
    end
    if energy != []
        gsen = AffExpr(0)
        Locb = bfind(tsupp, ltsupp, [1;4])
        if lattice == "chain"
            gsen += 3/4*mvar[Locb]
        else
            gsen += 3/2*mvar[Locb]
        end
        if J2 != 0
            Locb = bfind(tsupp, ltsupp, [1;7])
            gsen += 3/4*J2*mvar[Locb]
            if lattice == "square"
                Locb = bfind(tsupp, ltsupp, reduce!(UInt16[1;3*(slabel(L, 2, L=L)-1)+1], L=L, lattice="square")[1])
                gsen += 3/4*J2*mvar[Locb]
            end
        end
        @constraint(model, gsen>=energy[1])
        @constraint(model, gsen<=energy[2])
    end
    @objective(model, Min, obj)
    optimize!(model)
    status = termination_status(model)
    objv = objective_value(model)
    if status != MOI.OPTIMAL
       println("termination status: $status")
       status = primal_status(model)
       println("solution status: $status")
    end
    println("optimum = $objv")
    if correlation == true
        if lattice=="chain"
            cor0=zeros(Int(L/2))
            for i=1:Int(L/2)
                word=UInt16[1; 3*i+1]
                Locb=bfind(tsupp, ltsupp, word)
                cor0[i]=value(mvar[Locb])
            end
            cor1=zeros(Int(L/2-2))
            for i=3:Int(L/2)
                word=UInt16[1; 4; 3*(i-1)+1; 3*i+1]
                word=reduce!(word, L=L, lattice=lattice)[1]
                Locb=bfind(tsupp, ltsupp, word)
                cor1[i-2]=value(mvar[Locb])
            end
            cor2=zeros(Int(L/2-2))
            for i=3:Int(L/2)
                word=UInt16[1; 4; 3*i; 3*i+3]
                word=reduce!(word, L=L, lattice=lattice)[1]
                Locb=bfind(tsupp, ltsupp, word)
                cor2[i-2]=value(mvar[Locb])
            end
        else
            cor0=zeros(L, L)
            for i=1:L, j=1:L
                word=UInt16[1; 3*(slabel(i, j, L=L)-1)+1]
                word=reduce!(word, L=L, lattice=lattice)[1]
                Locb=bfind(tsupp, ltsupp, word)
                cor0[i,j]=value(mvar[Locb])
            end
            cor1,cor2=nothing,nothing
        end
    else
        cor0,cor1,cor2 = nothing,nothing,nothing
    end
    return objv,cor0,cor1,cor2
end

function bfind(A, l, a)
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if A[mid]==a
           return mid
        elseif A[mid]<a
            low=mid+1
        else
            high=mid-1
        end
    end
    return 0
end

function reduce1!(a::Vector{UInt16})
    la = length(a)
    flag = 1
    while flag == 1
        ind = findfirst(x->ceil(Int, a[x]/3) > ceil(Int, a[x+1]/3), 1:la-1)
        if ind != nothing
            temp = a[ind+1]
            a[ind+1] = a[ind]
            a[ind] = temp
            flag = 1
        else
            flag = 0
        end
    end
    return a
end

function reduce2!(a::Vector{UInt16})
    la = length(a)
    flag = 1
    coef = 1
    while flag == 1
        ind = findfirst(x->a[x]!=a[x+1]&&ceil(Int, a[x]/3)==ceil(Int, a[x+1]/3), 1:la-1)
        if ind != nothing
            if mod(a[ind], 3)==1&&mod(a[ind+1], 3)==2
                a[ind] += UInt16(2)
                coef *= im
            elseif mod(a[ind], 3)==2&&mod(a[ind+1], 3)==1
                a[ind] += UInt16(1)
                coef *= -im
            elseif mod(a[ind], 3)==1&&mod(a[ind+1], 3)==0
                a[ind] += UInt16(1)
                coef *= -im
            elseif mod(a[ind], 3)==0&&mod(a[ind+1], 3)==1
                a[ind] -= UInt16(1)
                coef *= im
            elseif mod(a[ind], 3)==2&&mod(a[ind+1], 3)==0
                a[ind] -= UInt16(1)
                coef *= im
            else
                a[ind] -= UInt16(2)
                coef *= -im
            end
            deleteat!(a, ind+1)
            la -= 1
            flag = 1
        else
            flag = 0
        end
    end
    return a,coef
end

function reduce3!(a::Vector{UInt16})
    i = 1
    while i < length(a)
        if a[i] == a[i+1]
            deleteat!(a, i)
            deleteat!(a, i)
        else
            i += 1
        end
    end
    return a
end

function iszero(a::Vector{UInt16})
    return any(i->isodd(count(isequal(i), mod.(a,3))), 0:2)
end

function reduce4(a::Vector{UInt16}, L; lattice="chain")
    l = length(a)
    if l > 0
        pa = Vector{UInt16}[]
        if lattice == "chain"
            for i = 1:l
                ta = [a[i:end];a[1:i-1].+3*L].-3*(ceil(Int, a[i]/3)-1)
                append!(pa, perm(ta))
                M = ceil(Int, ta[l]/3)
                ma = [3*(M-ceil(Int,ta[l+1-j]/3))+smod(ta[l+1-j],3) for j=1:l]
                append!(pa, perm(ma))
            end
        else
            loc = [location(ceil(Int, a[i]/3)) for i=1:l]
            for i = 1:l
                temp = zeros(UInt16, l)
                factor = [[1;1], [-1;1], [1;-1], [-1;-1]]
                for k = 1:4
                    for j = 1:l
                        p = slabel(factor[k][1]*(loc[j][1]-loc[i][1])+1, factor[k][2]*(loc[j][2]-loc[i][2])+1, L=L)
                        temp[j] = 3*p+a[j]-3*ceil(Int, a[j]/3)
                    end
                    append!(pa, perm(sort(temp)))
                end
            end
        end
        return findmin(pa)[1]
    else
        return a
    end
end

function reduce5(a::Vector{UInt16}, L)
    l = length(a)
    ra = zeros(UInt16, l)
    for j = 1:l
        loc = location(ceil(Int, a[j]/3))
        ra[j] = 3*slabel(loc[2], loc[1], L=L)+a[j]-3*ceil(Int, a[j]/3)
    end
    return min(a, ra)
end

function perm(a)
    ra = smod.(a, 3)
    sym = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2], [3;2;1]]
    return [UInt16.(3*(ceil.(Int, a./3).-1) .+ sym[i][ra]) for i=1:6]
end

function reduce!(a::Vector{UInt16}; L=0, lattice="chain", symmetry=true)
    reduce1!(a)
    reduce3!(a)
    a,coef = reduce2!(a)
    reduce3!(a)
    if iszero(a)
        coef = 0
    elseif symmetry == true
        a = reduce4(a, L, lattice=lattice)
    end
    if lattice == "square"
        a = reduce5(a, L)
    end
    return a,coef
end

function slabel(i, j; L=0)
    i = mod(i, L)==0 ? L : mod(i, L)
    j = mod(j, L)==0 ? L : mod(j, L)
    r = max(i,j)
    return r == i ? (r-1)^2+j : r^2+1-i
end

function location(p)
    r = ceil(Int, sqrt(p))
    if p-(r-1)^2 <= r
        return r, p-(r-1)^2
    else
        return r^2+1-p, r
    end
end

function rot(label)
    if label == 1
        return 2,3
    elseif label == 2
        return 3,1
    else
        return 1,2
    end
end

function split_basis(L, d, label; lattice="chain", extra=0, three_type=[1;1])
    if label>0
        basis=Vector{UInt16}[]
        if lattice=="chain"
            for i=1:L
                push!(basis, [3*(i-1)+label])
            end
        else
            for i=1:L, j=1:L
                push!(basis, [3*(slabel(j, i+j-1, L=L)-1)+label])
            end
        end
        if d>1
            a=[[rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1]]]
            for k=1:2
                if lattice=="chain"
                    for s=0:extra, i=1:L
                        push!(basis, [3*(i-1)+a[k][1];smod(3*(i+s)+a[k][2], 3*L)])
                    end
                else
                    tb = [[1;0], [0;1], [1;1], [1;-1], [2;0], [0;2], [2;1], [1;2], [1;-2], [2;-1], [2;2]]
                    for s=1:11, i=1:L, j=1:L
                        push!(basis, [3*(slabel(j, i+j-1, L=L)-1)+a[k][1];3*(slabel(j+tb[s][1], i+j-1+tb[s][2], L=L)-1)+a[k][2]])
                    end
                end
            end
        end
        if d>2
            if lattice=="chain"
                for i=1:L
                    push!(basis, [3*(i-1)+label;smod(3*(i-1+three_type[1])+label, 3*L);smod(3*(i-1+sum(three_type))+label, 3*L)])
                end
                for k=1:3, l=1:2
                    a=rot(label)[l]*ones(Int, 3)
                    a[k]=label
                    for i=1:L
                        push!(basis, [3*(i-1)+a[1];smod(3*(i-1+three_type[1])+a[2], 3*L);smod(3*(i-1+sum(three_type))+a[3], 3*L)])
                    end
                end
            else
                tb = [[0;1;1;1], [0;1;-1;1], [1;0;1;1], [-1;0;-1;1]]
                for s=1:4, i=1:L, j=1:L
                    push!(basis, [3*(slabel(j, i+j-1, L=L)-1)+label;3*(slabel(j+tb[s][1], i+j-1+tb[s][2], L=L)-1)+label;3*(slabel(j+tb[s][3], i+j-1+tb[s][4], L=L)-1)+label])
                end
                for k=1:3, l=1:2
                    a=rot(label)[l]*ones(Int, 3)
                    a[k]=label
                    for s=1:4, i=1:L, j=1:L
                        push!(basis, [3*(slabel(j, i+j-1, L=L)-1)+a[1];3*(slabel(j+tb[s][1], i+j-1+tb[s][2], L=L)-1)+a[2];3*(slabel(j+tb[s][3], i+j-1+tb[s][4], L=L)-1)+a[3]])
                    end
                end
            end
        end
        if d>3
            a=[[label;label;rot(label)[1];rot(label)[2]], [label;rot(label)[1];label;rot(label)[2]], [label;rot(label)[1];rot(label)[2];label], [rot(label)[1];label;label;rot(label)[2]], [rot(label)[1];label;rot(label)[2];label], [rot(label)[1];rot(label)[2];label;label],
            [label;label;rot(label)[2];rot(label)[1]], [label;rot(label)[2];label;rot(label)[1]], [label;rot(label)[2];rot(label)[1];label], [rot(label)[2];label;label;rot(label)[1]], [rot(label)[2];label;rot(label)[1];label], [rot(label)[2];rot(label)[1];label;label],
            [rot(label)[1];rot(label)[1];rot(label)[1];rot(label)[2]], [rot(label)[1];rot(label)[1];rot(label)[2];rot(label)[1]], [rot(label)[1];rot(label)[2];rot(label)[1];rot(label)[1]], [rot(label)[2];rot(label)[1];rot(label)[1];rot(label)[1]],
            [rot(label)[2];rot(label)[2];rot(label)[2];rot(label)[1]], [rot(label)[2];rot(label)[2];rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1];rot(label)[2];rot(label)[2]], [rot(label)[1];rot(label)[2];rot(label)[2];rot(label)[2]]]
            if lattice=="chain"
                for k=1:20, i=1:L
                    push!(basis, [3*(i-1)+a[k][1];smod(3*i+a[k][2], 3*L);smod(3*(i+1)+a[k][3], 3*L);smod(3*(i+2)+a[k][4], 3*L)])
                end
            else
                for k=1:20, i=1:L, j=1:L
                    push!(basis, [3*(slabel(j, i+j-1, L=L)-1)+a[k][1];3*(slabel(j, i+j, L=L)-1)+a[k][2];3*(slabel(j+1, i+j, L=L)-1)+a[k][3];3*(slabel(j+1, i+j-1, L=L)-1)+a[k][4]])
                end
            end
        end
    else
        basis=[UInt16[]]
        if d>1
            for k=1:3
                if lattice=="chain"
                    for s=0:extra, i=1:L
                        push!(basis, UInt16[3*(i-1)+k;smod(3*(i+s)+k, 3*L)])
                    end
                else
                    tb = [[1;0], [0;1], [1;1], [1;-1], [2;0], [0;2], [2;1], [1;2], [1;-2], [2;-1], [2;2]]
                    for s=1:11, i=1:L, j=1:L
                        push!(basis, UInt16[3*(slabel(j, i+j-1, L=L)-1)+k;3*(slabel(j+tb[s][1], i+j-1+tb[s][2], L=L)-1)+k])
                    end
                end
            end
        end
        if d>2
            a=[[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2], [3;2;1]]
            if lattice=="chain"
                for k=1:6, i=1:L
                    push!(basis, UInt16[3*(i-1)+a[k][1];smod(3*(i-1+three_type[1])+a[k][2], 3*L);smod(3*(i-1+sum(three_type))+a[k][3], 3*L)])
                end
            else
                tb = [[0;1;1;1], [0;1;-1;1], [1;0;1;1], [-1;0;-1;1]]
                for s=1:4, k=1:6, i=1:L, j=1:L
                    push!(basis, [3*(slabel(j, i+j-1, L=L)-1)+a[k][1];3*(slabel(j+tb[s][1], i+j-1+tb[s][2], L=L)-1)+a[k][2];3*(slabel(j+tb[s][3], i+j-1+tb[s][4], L=L)-1)+a[k][3]])
                end
            end
        end
        if d>3
            a=[[1;1;1;1], [2;2;2;2], [3;3;3;3], [1;1;2;2], [1;2;1;2], [1;2;2;1], [2;1;1;2], [2;1;2;1], [2;2;1;1],
            [1;1;3;3], [1;3;1;3], [1;3;3;1], [3;1;1;3], [3;1;3;1], [3;3;1;1], [3;3;2;2], [3;2;3;2], [3;2;2;3], [2;3;3;2], [2;3;2;3], [2;2;3;3]]
            if lattice=="chain"
                for k=1:21, i=1:L
                    push!(basis, [3*(i-1)+a[k][1];smod(3*i+a[k][2], 3*L);smod(3*(i+1)+a[k][3], 3*L);smod(3*(i+2)+a[k][4], 3*L)])
                end
            else
                for k=1:21, i=1:L, j=1:L
                    push!(basis, [3*(slabel(j, i+j-1, L=L)-1)+a[k][1];3*(slabel(j, i+j, L=L)-1)+a[k][2];3*(slabel(j+1, i+j, L=L)-1)+a[k][3];3*(slabel(j+1, i+j-1, L=L)-1)+a[k][4]])
                end
            end
        end
    end
    return basis
end

function smod(i, s)
    r = mod(i, s)
    return r == 0 ? s : r
end

function get_ncbasis(n, d; ind=UInt16[i for i=1:n])
    basis=[UInt16[]]
    for i=1:d
        append!(basis, _get_ncbasis_deg(n, i, ind=ind))
    end
    return basis
end

function _get_ncbasis_deg(n, d; ind=UInt16[i for i=1:n])
    if d>0
        basis=Vector{UInt16}[]
        for i=1:n
            temp=_get_ncbasis_deg(n, d-1, ind=ind)
            push!.(temp, ind[i])
            append!(basis, temp)
        end
        return basis
    else
        return [UInt16[]]
    end
end

function GSE2(supp::Vector{Vector{UInt16}}, coe::Vector{Float64}, L::Int, d::Int; lattice="square", A=[], clique=3, solver="COSMO", CS=true, TS="block", merge=false, QUIET=false)
    # basis=get_ncbasis(3*L^2, d)
    # basis=basis[is_basis.(basis, L, lattice=lattice)]
    # tsupp=copy(supp)
    # sort!(tsupp)
    # blocks,cl,blocksize=get_ncblocks(tsupp,basis,TS=TS,QUIET=QUIET,merge=merge)
    # if CS==true
        if lattice=="square"
            cliques=Vector{Vector{UInt16}}(undef, Int(L^2/2))
            l=1
            for i=1:2:L-1, j=1:2:L-1
                cliques[l]=UInt16[]
                for k=1:3
                    push!(cliques[l], 3*(slabel(i-1,j,L=L)-1)+k, 3*(slabel(i,j-1,L=L)-1)+k, 3*(slabel(i,j,L=L)-1)+k, 3*(slabel(i+1,j,L=L)-1)+k, 3*(slabel(i,j+1,L=L)-1)+k)
                end
                sort!(cliques[l])
                l+=1
            end
            for i=2:2:L, j=2:2:L
                cliques[l]=UInt16[]
                for k=1:3
                    push!(cliques[l], 3*(slabel(i-1,j,L=L)-1)+k, 3*(slabel(i,j-1,L=L)-1)+k, 3*(slabel(i,j,L=L)-1)+k, 3*(slabel(i+1,j,L=L)-1)+k, 3*(slabel(i,j+1,L=L)-1)+k)
                end
                sort!(cliques[l])
                l+=1
            end
        elseif lattice=="chain"
            cliques=Vector{Vector{UInt16}}(undef, L)
            for i=1:L
                cliques[i]=UInt16[]
                for j=0:clique-1
                    push!(cliques[i], 3*smod(i+j, L)-2, 3*smod(i+j, L)-1, 3*smod(i+j, L))
                end
            end
            sort!(cliques[i])
        else
            cliques=get_cliques(A, CS=CS)
            # cliques=Vector{Vector{UInt16}}(undef, L^2)
            # for i=1:L, j=1:L
                # cliques[(i-1)*L+j]=UInt16[]
                # cliques[(i-1)*L+j]=get_clique(i, j, L)
                # for l=1:5
                #     push!(cliques[(i-1)*L+j], 3*clique[l]-2, 3*clique[l]-1, 3*clique[l])
                # end
                # sort!(cliques[(i-1)*L+j])
            # end
        end
    # else
    #     cliques=[UInt16[i for i=1:3*L^2]]
    # end
    blocks,cl,blocksize,basis=get_blocks_mix(L, d, supp, cliques, A=A, lattice=lattice, TS=TS)
    opt=blockpop(supp, coe, basis, blocks, cl, blocksize, solver=solver, QUIET=QUIET)
    return opt
end

function get_blocks_mix(L, d, supp, cliques; A=[], lattice="chain", TS="block")
    cql=length(cliques)
    blocks=Vector{Vector{Vector{UInt16}}}(undef,cql)
    cl=Vector{UInt16}(undef,cql)
    blocksize=Vector{Vector{Int}}(undef,cql)
    basis=Vector{Vector{Vector{UInt16}}}(undef,cql)
    for i=1:cql
        if lattice=="chain"||lattice=="square"
            ind=[issubset(supp[j], cliques[i]) for j=1:length(supp)]
            basis[i]=get_ncbasis(length(cliques[i]), d, ind=cliques[i])
            basis[i]=basis[i][is_basis.(basis[i], L, lattice=lattice)]
        else
            varclique=UInt16[]
            for l=1:length(cliques[i])
                push!(varclique, 3*cliques[i][l]-2, 3*cliques[i][l]-1, 3*cliques[i][l])
            end
            sort!(varclique)
            ind=[issubset(supp[j], varclique) for j=1:length(supp)]
            if cql==1
                cliques0=get_cliques(A)
                basis[1]=Vector{UInt16}[]
                for j=1:length(cliques0)
                    append!(basis[1], Kagome_basis(cliques0[j], d))
                end
                unique!(basis[1])
            else
                basis[i]=Kagome_basis(cliques[i], d)
            end
        end
        tsupp=copy(supp[ind])
        sort!(tsupp)
        blocks[i],cl[i],blocksize[i]=get_ncblocks(tsupp,basis[i],TS=TS,QUIET=true)
    end
    return blocks,cl,blocksize,basis
end

function get_ncgraph(tsupp, basis)
    lb=length(basis)
    G=SimpleGraph(lb)
    ltsupp=length(tsupp)
    for i = 1:lb, j = i+1:lb
        bi = [basis[i][end:-1:1]; basis[j]]
        bi=reduce!(bi, symmetrty=false)[1]
        if bfind(tsupp, ltsupp, bi)!=0
           add_edge!(G, i, j)
        end
    end
    return G
end

function get_ncblocks(tsupp, basis; TS="block", QUIET=true, merge=false)
    G=get_ncgraph(tsupp, basis)
    if TS=="block"
        blocks=connected_components(G)
        blocksize=length.(blocks)
        cl=length(blocksize)
    else
        blocks,cl,blocksize=chordal_cliques!(G, method=TS, minimize=false)
        if merge==true
            blocks,cl,blocksize=clique_merge!(blocks, cl, QUIET=true)
        end
    end
    ub=unique(blocksize)
    sizes=[sum(blocksize.== i) for i in ub]
    if QUIET==false
        println("------------------------------------------------------")
        println("The sizes of blocks:\n$ub\n$sizes")
        println("------------------------------------------------------")
    end
    return blocks,cl,blocksize
end

function is_basis(a::Vector{UInt16}, L::Int; lattice="chain")
    for i=1:length(a)-1
        if ceil(Int, a[i]/3)>=ceil(Int, a[i+1]/3)
            return false
        end
    end
    if lattice=="chain"
        if length(a)==2
            if ceil(Int, a[2]/3) - ceil(Int, a[1]/3)!=1&&ceil(Int, a[2]/3) - ceil(Int, a[1]/3)!=L-1
                return false
            end
        elseif length(a)==3
            if (ceil(Int, a[1]/3)==1&&ceil(Int, a[2]/3)==L-1&&ceil(Int, a[3]/3)==L)||(ceil(Int, a[1]/3)==1&&ceil(Int, a[2]/3)==2&&ceil(Int, a[3]/3)==L)
                return true
            elseif ceil(Int, a[3]/3) - ceil(Int, a[2]/3) > 1||ceil(Int, a[2]/3) - ceil(Int, a[1]/3) > 1
                return false
            end
        end
    elseif lattice=="square"
        if length(a)==2
            i1,j1=location(ceil(Int, a[1]/3))
            i2,j2=location(ceil(Int, a[2]/3))
            if (i1!=i2&&j1!=j2)||(i1==i2&&abs(j1-j2)!=1&&abs(j1-j2)!=L-1)||(j1==j2&&abs(i1-i2)!=1&&abs(i1-i2)!=L-1)
                return false
            end
        end
        if length(a)==3
            i1,j1=location(ceil(Int, a[1]/3))
            i2,j2=location(ceil(Int, a[2]/3))
            i3,j3=location(ceil(Int, a[3]/3))
            if i1==i2&&i2==i3
                temp=sort([j1;j2;j3])
                if temp==[1;L-1;L]||temp==[1;2;L]
                    return true
                elseif temp[2]-temp[1]>1||temp[3]-temp[2]>1
                    return false
                end
            elseif j1==j2&&j2==j3
                temp=sort([i1;i2;i3])
                if temp==[1;L-1;L]||temp==[1;2;L]
                    return true
                elseif temp[2]-temp[1]>1||temp[3]-temp[2]>1
                    return false
                end
            else
                return false
            end
        end
    else
        if length(a)==2
            ind1=ceil(Int, a[1]/3)
            ind2=ceil(Int, a[2]/3)
            i1,j1=location(ceil(Int, ind1/3))
            i2,j2=location(ceil(Int, ind2/3))
            if i1==i2&&j1>j2&&(mod(ind1, 3)!=1||mod(ind2, 3)!=0)
                return false
            elseif i1==i2&&j1<j2&&(mod(ind1, 3)!=0||mod(ind2, 3)!=1)
                return false
            elseif j1==j2&&i1>i2&&(mod(ind1, 3)!=1||mod(ind2, 3)!=2)
                return false
            elseif j1==j2&&i1<i2&&(mod(ind1, 3)!=2||mod(ind2, 3)!=1)
                return false
            elseif i1==i2+1&&j1==j2-1&&(mod(ind1, 3)!=0||mod(ind2, 3)!=2)
                return false
            elseif i2==i1+1&&j2==j1-1&&(mod(ind1, 3)!=2||mod(ind2, 3)!=0)
                return false
            elseif i1!=i2&&j1!=j2&&(i1!=i2+1||j1!=j2-1)&&(i2!=i1+1||j2!=j1-1)
                return false
            end
        end
    end
    return true
end

function blockpop(supp, coe, basis, blocks, cl, blocksize; solver="Mosek", QUIET=false)
    tsupp=Vector{UInt16}[]
    cql=length(cl)
    for k=1:cql, i=1:cl[k], j=1:blocksize[k][i], r=j:blocksize[k][i]
        @inbounds bi=[basis[k][blocks[k][i][j]]; basis[k][blocks[k][i][r]]]
        bi=reduce!(bi, symmetrty=false)[1]
        push!(tsupp, bi)
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp=length(tsupp)
    if solver=="COSMO"
        model=Model(optimizer_with_attributes(COSMO.Optimizer))
    else
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons=[AffExpr(0) for i=1:ltsupp]
    for k=1:cql, i=1:cl[k]
        bs=blocksize[k][i]
        pos=@variable(model, [1:2*bs, 1:2*bs], PSD)
        for j=1:bs, r=j:bs
            @constraint(model, pos[j,r]==pos[j+bs,r+bs])
            @constraint(model, pos[r,j+bs]+pos[j,r+bs]==0)
            @inbounds bi = [basis[k][blocks[k][i][j]]; basis[k][blocks[k][i][r]]]
            bi,coef=reduce!(bi, symmetrty=false)
            Locb=bfind(tsupp, ltsupp, bi)
            if r==j
                @inbounds add_to_expression!(cons[Locb], pos[j,r])
            else
                if coef==1
                    @inbounds add_to_expression!(cons[Locb], 2,  pos[j,r])
                elseif coef==im
                    @inbounds add_to_expression!(cons[Locb], -2,  pos[j,r+bs])
                elseif coef==-1
                    @inbounds add_to_expression!(cons[Locb], -2,  pos[j,r])
                else
                    @inbounds add_to_expression!(cons[Locb], 2, pos[j,r+bs])
                end
            end
        end
    end
    bc=zeros(ltsupp)
    for i=1:length(supp)
        Locb=bfind(tsupp, ltsupp, supp[i])
        if Locb==0
           @error "The monomial basis is not enough!"
           return nothing
        else
           bc[Locb]=coe[i]
        end
    end
    @variable(model, lower)
    cons[1]+=lower
    @constraint(model, con[i=1:ltsupp], cons[i]==bc[i])
    @objective(model, Max, lower)
    optimize!(model)
    status=termination_status(model)
    objv = objective_value(model)
    if status!=MOI.OPTIMAL
       println("termination status: $status")
       status=primal_status(model)
       println("solution status: $status")
    end
    println("optimum = $objv")
    return objv
end

function get_clique(i, j, L)
    clique=UInt16[]
    ind0=slabel(i,j,L=L)
    ind1=slabel(i-1,j+1,L=L)
    ind2=slabel(i,j+1,L=L)
    push!(clique, 3*ind0-2, 3*ind0-1, 3*ind0, 3*ind1-1, 3*ind2-2)
    return clique
end

function get_points(A)
    points=Vector{Int}[]
    for i=minimum(A[1,:]):maximum(A[1,:]), j=minimum(A[2,:]):maximum(A[2,:])
        if !isodd(i)||!iseven(j)
            μ=A\[i;j]
            if μ[1]>=0&&μ[1]<1&&μ[2]>=0&&μ[2]<1
                push!(points, [i;j])
            end
        end
    end
    sort!(points)
    ind=[iseven(points[i][2]) for i=1:length(points)]
    basepoints=points[ind]
    return points,basepoints
end

function get_cliques(A; CS=true)
    points,basepoints=get_points(A)
    lp=length(points)
    if CS==true
        cql=length(basepoints)
        cliques=Vector{Vector{UInt16}}(undef, cql)
        for i=1:cql
            mem1=bfind(points, lp, normform!(A, [basepoints[i][1];basepoints[i][2]-1]))
            mem2=bfind(points, lp, normform!(A, [basepoints[i][1]+1;basepoints[i][2]-1]))
            mem3=bfind(points, lp, basepoints[i])
            mem4=bfind(points, lp, normform!(A, [basepoints[i][1];basepoints[i][2]+1]))
            mem5=bfind(points, lp, normform!(A, [basepoints[i][1]-1;basepoints[i][2]+1]))
            cliques[i]=[mem1;mem2;mem3;mem4;mem5]
        end
    else
        cliques=[UInt16[i for i=1:lp]]
    end
    return cliques
end

function normform!(A, point)
    μ=A\point
    for i=1:2
        if μ[i]<0
            point+=A[:,i]
        elseif μ[i]>=1
            point-=A[:,i]
        end
    end
    return point
end

function clique_decomp(L::Int,supp::Vector{Vector{UInt16}};alg="MF",minimize=true)
    G=SimpleGraph(L^2)
    for i=1:L, j=1:L
        add_edge!(G, slabel(i,j,L=L), slabel(i+1,j,L=L))
        add_edge!(G, slabel(i,j,L=L), slabel(i,j+1,L=L))
    end
    cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("------------------------------------------------------")
    return cliques,cql,cliquesize
end

function Kagome_basis(clique, d)
    basis=Vector{UInt16}[[]]
    for i=1:length(clique), j=1:3
        push!(basis, [3*(clique[i]-1)+j])
    end
    if d>1
        for i=1:3, j=1:3
            push!(basis, sort([3*(clique[1]-1)+i;3*(clique[2]-1)+j]), sort([3*(clique[1]-1)+i;3*(clique[3]-1)+j]), sort([3*(clique[2]-1)+i;3*(clique[3]-1)+j]))
            push!(basis, sort([3*(clique[3]-1)+i;3*(clique[4]-1)+j]), sort([3*(clique[3]-1)+i;3*(clique[5]-1)+j]), sort([3*(clique[4]-1)+i;3*(clique[5]-1)+j]))
        end
    end
    if d>2
        for i=1:3, j=1:3, k=1:3
            push!(basis, sort([3*(clique[1]-1)+i;3*(clique[2]-1)+j;3*(clique[3]-1)+k]))
            push!(basis, sort([3*(clique[3]-1)+i;3*(clique[4]-1)+j;3*(clique[5]-1)+k]))
        end
    end
    return basis
end

# function GSE0(supp::Vector{Vector{UInt16}}, coe::Vector{Float64}, L::Int, d::Int; QUIET=false, lattice="chain", totalspin=false, sector=0, correlation=false)
#     basis=Vector{Vector{Vector{UInt16}}}(undef, 4)
#     tsupp=Vector{UInt16}[]
#     for i=0:3
#         basis[i+1]=split_basis(L, d, i, lattice=lattice)
#         for j=1:length(basis[i+1]), k=j:length(basis[i+1])
#             @inbounds bi=[basis[i+1][j]; basis[i+1][k]]
#             bi=reduce!(bi, L=L, lattice=lattice)[1]
#             push!(tsupp, bi)
#         end
#     end
#     sort!(tsupp)
#     unique!(tsupp)
#     ltsupp=length(tsupp)
#     model=Model(optimizer_with_attributes(Mosek.Optimizer))
#     set_optimizer_attribute(model, MOI.Silent(), QUIET)
#     cons=[AffExpr(0) for i=1:ltsupp]
#     for i=1:4
#         bs=length(basis[i])
#         pos=@variable(model, [1:2*bs, 1:2*bs], PSD)
#         for j=1:bs, r=j:bs
#             @constraint(model, pos[j,r]==pos[j+bs,r+bs])
#             @constraint(model, pos[r,j+bs]+pos[j,r+bs]==0)
#             @inbounds bi = [basis[i][j]; basis[i][r]]
#             bi,coef=reduce!(bi, L=L, lattice=lattice)
#             Locb=bfind(tsupp, ltsupp, bi)
#             if r==j
#                 @inbounds cons[Locb]+=pos[j,r]
#             else
#                 if coef==1
#                     @inbounds cons[Locb]+=2*pos[j,r]
#                 elseif coef==im
#                     @inbounds cons[Locb]-=2*pos[j,r+bs]
#                 elseif coef==-1
#                     @inbounds cons[Locb]-=2*pos[j,r]
#                 else
#                     @inbounds cons[Locb]+=2*pos[j,r+bs]
#                 end
#             end
#         end
#     end
#     bc=zeros(ltsupp)
#     for i=1:length(supp)
#         Locb=bfind(tsupp, ltsupp, supp[i])
#         if Locb==0
#            @error "The monomial basis is not enough!"
#            return nothing,nothing
#         else
#            bc[Locb]=coe[i]
#         end
#     end
#     @variable(model, lower)
#     cons[1]+=lower
#     @constraint(model, con[i=1:ltsupp], cons[i]==bc[i])
#     @objective(model, Max, lower)
#     optimize!(model)
#     status=termination_status(model)
#     objv = objective_value(model)
#     if status!=MOI.OPTIMAL
#        println("termination status: $status")
#        status=primal_status(model)
#        println("solution status: $status")
#     end
#     println("optimum = $objv")
#     if correlation==true
#         maxl=min(Int(n/3), ceil(Int, L/2))
#         dual_var=-dual.(con)
#         cor=zeros(3, maxl-1)
#         for i=1:maxl-1, j=1:3
#             word=UInt16[j; 3*i+j]
#             Locb=bfind(tsupp, ltsupp, word)
#             cor[j, i]=dual_var[Locb]
#         end
#     else
#         cor=nothing
#     end
#     return objv,cor
# end
