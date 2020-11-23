function GSE(supp::Vector{Vector{UInt16}}, coe::Vector{Float64}, L::Int, d::Int; QUIET=false, dim=1, totalspin=false, correlation=true)
    basis=Vector{Vector{Vector{UInt16}}}(undef, 4)
    tsupp=Vector{UInt16}[]
    for i=0:3
        basis[i+1]=split_basis(L, d, i)
        for j=1:length(basis[i+1]), k=j:length(basis[i+1])
            @inbounds bi=[basis[i+1][j]; basis[i+1][k]]
            bi,coef=reduce!(bi, L=L, dim=dim)
            if coef!=0
                push!(tsupp, bi)
            end
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp=length(tsupp)
    model=Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    mvar=@variable(model, [1:ltsupp])
    for i=0:3
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
                    bi,coef=reduce!(basis[1][L*(l-1)+2], L=L, dim=dim)
                    if coef!=0
                        Locb=ncbfind(tsupp, ltsupp, bi)
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
                    bi,coef[r]=reduce!(bi, L=L, dim=dim)
                    if coef[r]!=0
                        Locb[r]=ncbfind(tsupp, ltsupp, bi)
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
                        reig[l+1][r1, r2]+=(-1)^l*mvar[Locb[end]]
                    end
                    for r=1:Int(L/2)-1
                        if coef[r]^2==1&&abs(cos(2*pi*r*l/L))>=1e-8
                            reig[l+1][r1, r2]+=2*coef[r]*cos(2*pi*r*l/L)*mvar[Locb[r]]
                        elseif coef[r]==im&&abs(sin(2*pi*r*l/L))>=1e-8
                            reig[l+1][r1, r2]-=2*sin(2*pi*r*l/L)*mvar[Locb[r]]
                        elseif coef[r]==-im&&abs(sin(2*pi*r*l/L))>=1e-8
                            reig[l+1][r1, r2]+=2*sin(2*pi*r*l/L)*mvar[Locb[r]]
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
                    bi,coef[r+1]=reduce!(bi, L=L, dim=dim)
                    if coef[r+1]!=0
                        Locb[r+1]=ncbfind(tsupp, ltsupp, bi)
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
                                reig[l+1][r1, r2]+=coef[r+1]*cos(2*pi*r*l/L)*mvar[Locb[r+1]]
                            end
                            if abs(sin(2*pi*r*l/L))>=1e-8
                                ieig[l+1][r1, r2]+=coef[r+1]*sin(2*pi*r*l/L)*mvar[Locb[r+1]]
                            end
                        elseif coef[r+1]==im
                            if abs(sin(2*pi*r*l/L))>=1e-8
                                reig[l+1][r1, r2]-=sin(2*pi*r*l/L)*mvar[Locb[r+1]]
                            end
                            if abs(cos(2*pi*r*l/L))>=1e-8
                                ieig[l+1][r1, r2]+=cos(2*pi*r*l/L)*mvar[Locb[r+1]]
                            end
                        elseif coef[r+1]==-im
                            if abs(sin(2*pi*r*l/L))>=1e-8
                                reig[l+1][r1, r2]+=sin(2*pi*r*l/L)*mvar[Locb[r+1]]
                            end
                            if abs(cos(2*pi*r*l/L))>=1e-8
                                ieig[l+1][r1, r2]-=cos(2*pi*r*l/L)*mvar[Locb[r+1]]
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
    if totalspin==true
        J1=AffExpr(3*L)
        for i=1:Int(L/2)-1
            word=UInt16[1; 3*i+1]
            Locb=ncbfind(tsupp, ltsupp, word)
            J1+=6*L*mvar[Locb]
        end
        Locb=ncbfind(tsupp, ltsupp, UInt16[1; 3*Int(L/2)+1])
        J1+=3*L*mvar[Locb]
        @constraint(model, J1==0)
        # J2=AffExpr(0)
        # for j1=1:L, j2=1:L, k1=1:L, k2=1:L
        #     temp=UInt16[3*(j1-1)+1;3*(k1-1)+1;3*(j2-1)+1;3*(k2-1)+1]
        #     bi=reduce!(temp, L=L, dim=dim)[1]
        #     Locb=ncbfind(tsupp,ltsupp,bi)
        #     J2+=3*mvar[Locb]
        #     temp=UInt16[3*(j1-1)+1;3*(k1-1)+1;3*(j2-1)+2;3*(k2-1)+2]
        #     bi,coef=reduce!(temp, L=L, dim=dim)
        #     Locb=ncbfind(tsupp,ltsupp,bi)
        #     if coef^2==1
        #         J2+=6*coef*mvar[Locb]
        #     end
        # end
        # @constraint(model, J2==0)
    end
    obj=AffExpr(0)
    for i=1:length(supp)
        Locb=ncbfind(tsupp,ltsupp,supp[i])
        obj+=coe[i]*mvar[Locb]
    end
    @constraint(model, mvar[1]==1)
    @objective(model, Min, obj)
    optimize!(model)
    status=termination_status(model)
    objv = objective_value(model)
    if status!=MOI.OPTIMAL
       println("termination status: $status")
       status=primal_status(model)
       println("solution status: $status")
    end
    println("optimum = $objv")
    if correlation==true
        cor=zeros(3, Int(L/2))
        for i=1:Int(L/2), j=1:3
            word=UInt16[j; 3*i+j]
            Locb=ncbfind(tsupp, ltsupp, word)
            cor[j, i]=value(mvar[Locb])
        end
    else
        cor=nothing
    end
    return objv,cor
end

function ncbfind(A, l, a)
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if isequal(A[mid], a)
           return mid
        elseif isless(A[mid], a)
            low=mid+1
        else
            high=mid-1
        end
    end
    return 0
end

function reduce1!(a::Vector{UInt16})
    la=length(a)
    flag=1
    while flag==1
        ind=findfirst(x->a[x]>a[x+1]&&ceil(Int, a[x]/3)!=ceil(Int, a[x+1]/3), 1:la-1)
        if ind!=nothing
            temp=a[ind+1]
            a[ind+1]=a[ind]
            a[ind]=temp
            flag=1
        else
            flag=0
        end
    end
    return a
end

function reduce2!(a::Vector{UInt16})
    la=length(a)
    flag=1
    coef=1
    while flag==1
        ind=findfirst(x->a[x]!=a[x+1]&&ceil(Int, a[x]/3)==ceil(Int, a[x+1]/3), 1:la-1)
        if ind!=nothing
            if mod(a[ind], 3)==1&&mod(a[ind+1], 3)==2
                a[ind]+=UInt16(2)
                coef*=im
            elseif mod(a[ind], 3)==2&&mod(a[ind+1], 3)==1
                a[ind]+=UInt16(1)
                coef*=-im
            elseif mod(a[ind], 3)==1&&mod(a[ind+1], 3)==0
                a[ind]+=UInt16(1)
                coef*=-im
            elseif mod(a[ind], 3)==0&&mod(a[ind+1], 3)==1
                a[ind]-=UInt16(1)
                coef*=im
            elseif mod(a[ind], 3)==2&&mod(a[ind+1], 3)==0
                a[ind]-=UInt16(1)
                coef*=im
            else
                a[ind]-=UInt16(2)
                coef*=-im
            end
            deleteat!(a, ind+1)
            la-=1
            flag=1
        else
            flag=0
        end
    end
    return a,coef
end

function reduce3!(a::Vector{UInt16})
    i=1
    while i<length(a)
        if a[i]==a[i+1]
            deleteat!(a, i)
            deleteat!(a, i)
        else
            i+=1
        end
    end
    return a
end

function check(a::Vector{UInt16})
    b=mod.(a, 3)
    if any(i->isodd(count(isequal(i), b)), 0:2)
        return false
    end
    return true
end

function reduce4!(a::Vector{UInt16}; L=0, dim=1)
    l=length(a)
    if dim==1&&l>0
        loc=UInt16[ceil(UInt16, a[i]/3) for i=1:l]
        b=[[loc[i:end];loc[1:i-1].+UInt16(L)].-loc[i] for i=1:l]
        ind=findmini(b)
        a=[a[ind:end];a[1:ind-1].+UInt16(3*L)].-UInt16(3)*(loc[ind]-UInt16(1))
    else
        if l==1&&a[1]>3
            a.-=3*(ceil(UInt16, a[1]/3)-1)
        elseif l>1
            loc=[location(ceil(Int, a[i]/3)) for i=1:l]
            origin=findmin(loc)[1]
            if origin!=(1,1)
                for i=1:l
                    p=label(mod(loc[i][1]-origin[1]+1, L), mod(loc[i][2]-origin[2]+1, L), L=L)
                    a[i]=3*(p-1)+a[i]-3*(ceil(UInt16, a[i]/3)-1)
                end
            end
        end
    end
    return a
end

function findmini(b::Vector{Vector{UInt16}})
    ind = sum.(b).==minimum(sum.(b))
    if count(ind)==1
        return findfirst(ind)
    else
        return findmin(b[ind])[2]
    end
end

function reduce!(a::Vector{UInt16}; L=0, dim=1)
    reduce1!(a)
    reduce3!(a)
    a,coef=reduce2!(a)
    reduce3!(a)
    if check(a)
        a=reduce4!(a, L=L, dim=dim)
    else
        coef=0
    end
    return a,coef
end

function label(i, j; L=0)
    i=i==0 ? L : i
    j=j==0 ? L : j
    r=max(i,j)
    return r==i ? (r-1)^2+j : r^2+1-i
end

function location(p)
    r=ceil(Int, sqrt(p))
    if p-(r-1)^2<=r
        return r, p-(r-1)^2
    else
        return r^2+1-p, r
    end
end

function rot(label)
    if label==1
        return 2,3
    elseif label==2
        return 3,1
    else
        return 1,2
    end
end

function split_basis(L, d, label)
    if label>0
        basis=Vector{UInt16}[]
        for i=1:L
            push!(basis, [3*(i-1)+label])
        end
        if d>1
            a=[[rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1]]]
            for k=1:2
                for i=1:L
                    push!(basis, [gmod(3*(i-1)+a[k][1], L);gmod(3*i+a[k][2], L)])
                end
            end
        end
        if d>2
            for i=1:L-2
                push!(basis, [3*(i-1)+label;3*i+label;3*(i+1)+label])
            end
            push!(basis, [label;3*(L-2)+label;3*(L-1)+label], [label;3+label;3*(L-1)+label])
            for k=1:3, l=1:2
                a=rot(label)[l]*ones(Int, 3)
                a[k]=label
                for i=1:L-2
                    push!(basis, [3*(i-1)+a[1];3*i+a[2];3*(i+1)+a[3]])
                end
                push!(basis, [a[3];3*(L-2)+a[1];3*(L-1)+a[2]], [a[2];3+a[3];3*(L-1)+a[1]])
            end
        end
    else
        basis=[UInt16[]]
        if d>1
            for k=1:3
                for i=1:L-1
                    push!(basis, UInt16[3*(i-1)+k;3*i+k])
                end
                push!(basis, UInt16[k;3*(L-1)+k])
            end
        end
        if d>2
            a=[[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2], [3;2;1]]
            for k=1:6
                for i=1:L-2
                    push!(basis, UInt16[3*(i-1)+a[k][1];3*i+a[k][2];3*(i+1)+a[k][3]])
                end
                push!(basis, UInt16[a[k][3];3*(L-2)+a[k][1];3*(L-1)+a[k][2]], UInt16[a[k][2];3+a[k][3];3*(L-1)+a[k][1]])
            end
        end
    end
    return basis
end

function gmod(i, L)
    return i>3*L ? i-3*L : i
end

# function GSE(supp::Vector{Vector{UInt16}}, coe::Vector{Float64}, L::Int, n::Int, d::Int; dim=1, conse=true, TS="block", merge=false, QUIET=false, correlation=false)
#     basis=get_ncbasis(n, d)
#     basis=basis[is_basis.(basis, L=L, dim=dim, conse=conse)]
#     tsupp=copy(supp)
#     sort!(tsupp)
#     blocks,cl,blocksize,ub,sizes,_=get_ncblocks(tsupp,basis,TS=TS,QUIET=QUIET,merge=merge,field="complex",L=L,dim=dim,PBC=PBC)
#     opt,tsupp,cor=ncblockupop_com(n,supp,coe,basis,blocks,cl,blocksize,QUIET=QUIET,L=L,dim=dim,PBC=PBC,correlation=correlation)
#     data=data_type(supp,basis,coe,"eigen",tsupp,ub,sizes)
#     return opt,cor,data
# end
#
# function is_basis(a::Vector{UInt16}; L=0, dim=1, conse=true)
#     for i=1:length(a)-1
#         if (a[i]>a[i+1]&&ceil(Int, a[i]/3)!=ceil(Int, a[i+1]/3))||(a[i]!=a[i+1]&&ceil(Int, a[i]/3)==ceil(Int, a[i+1]/3))||a[i]==a[i+1]
#             return false
#         end
#     end
#     if dim==1&&conse==true
#         if length(a)==2
#             dis=ceil(Int, a[2]/3)-ceil(Int, a[1]/3)
#             if min(dis, L-dis) > 1
#             # if ceil(Int, a[2]/3) - ceil(Int, a[1]/3) > 1
#                 return false
#             end
#         elseif length(a)==3
#             if ceil(Int, a[3]/3) - ceil(Int, a[2]/3) > 1||ceil(Int, a[2]/3) - ceil(Int, a[1]/3) > 1
#                 return false
#             end
#         end
#     elseif dim==2&&conse==true
#         if length(a)==2
#             i1,j1=location(ceil(Int, a[1]/3))
#             i2,j2=location(ceil(Int, a[2]/3))
#             if abs(mod(i1, L)-mod(i2, L))+abs(mod(j1, L)-mod(j2, L))>1
#                 return false
#             end
#         end
#     end
#     return true
# end
#
# function ncblockupop_com(n,supp,coe,basis,blocks,cl,blocksize;QUIET=true,L=0,dim=1,PBC=false,correlation=false)
#     tsupp=Vector{Vector{UInt16}}(undef, Int(sum(blocksize.^2+blocksize)/2))
#     k=1
#     for i=1:cl, j=1:blocksize[i], r=j:blocksize[i]
#         @inbounds bi=[basis[blocks[i][j]]; basis[blocks[i][r]]]
#         @inbounds tsupp[k]=reduce!(bi, L=L, dim=dim, PBC=PBC)[1]
#         k+=1
#     end
#     sort!(tsupp)
#     unique!(tsupp)
#     ltsupp=length(tsupp)
#     model=Model(optimizer_with_attributes(Mosek.Optimizer))
#     set_optimizer_attribute(model, MOI.Silent(), QUIET)
#     cons=[AffExpr(0) for i=1:ltsupp]
#     for i=1:cl
#         bs=blocksize[i]
#         pos=@variable(model, [1:2*bs, 1:2*bs], PSD)
#         for j=1:bs, r=j:bs
#             @constraint(model, pos[j,r]==pos[j+bs,r+bs])
#             @constraint(model, pos[r,j+bs]+pos[j,r+bs]==0)
#             @inbounds bi = [basis[blocks[i][j]]; basis[blocks[i][r]]]
#             bi,coef=reduce!(bi, L=L, dim=dim, PBC=PBC)
#             Locb=ncbfind(tsupp,ltsupp,bi)
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
#         Locb=ncbfind(tsupp,ltsupp,supp[i])
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
#         if PBC==true
#             maxl=min(Int(n/3), ceil(Int, L/2))
#         else
#             maxl=Int(n/3)
#         end
#         dual_var=-dual.(con)
#         cor=Vector{Matrix{Float64}}(undef, cl)
#         for i=1:cl
#             cor[i]=zeros(blocksize[i], blocksize[i])
#             for j=1:blocksize[i],k=1:blocksize[i]
#                 bi = [basis[blocks[i][j]]; basis[blocks[i][k]]]
#                 bi,_=reduce!(bi, L=L, dim=dim, PBC=PBC)
#                 Locb=ncbfind(tsupp, ltsupp, bi)
#                 cor[i][j,k]=dual_var[Locb]
#             end
#         end
#         # cor=zeros(3, maxl-1)
#         # for i=1:maxl-1, j=1:3
#         #     word=UInt16[j; 3*i+j]
#         #     Locb=ncbfind(tsupp, ltsupp, word)
#         #     cor[j, i]=dual_var[Locb]
#         # end
#     else
#         cor=nothing
#     end
#     return objv,tsupp,cor
# end
