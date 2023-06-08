function bfind(A, l, a)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if A[mid] == a
           return mid
        elseif A[mid] < a
            low = mid + 1
        else
            high = mid - 1
        end
    end
    return nothing
end

function reduce1!(a::Vector{UInt16})
    la = length(a)
    flag = 1
    while flag == 1
        ind = findfirst(x->ceil(Int, a[x]/3) > ceil(Int, a[x+1]/3), 1:la-1)
        if ind !== nothing
            a[ind],a[ind+1] = a[ind+1],a[ind]
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
        if ind !== nothing
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

function isz(a::Vector{UInt16})
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
    if isz(a)
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
                    if L > 4
                        push!(tb, [2;-2])
                    end
                    if extra == true
                        if L == 6
                            push!(tb, [0;3], [1;3], [2;3], [3;3], [3;2], [3;1], [3;0], [3;-1], [3;-2])
                        end
                    end
                    for s=1:length(tb), i=1:L, j=1:L
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
                tb = [[0;1;1;1], [0;1;-1;1], [1;0;1;1], [-1;0;-1;1], [1;0;2;0], [0;1;0;2]]
                for s=1:length(tb), i=1:L, j=1:L
                    push!(basis, [3*(slabel(j, i+j-1, L=L)-1)+label;3*(slabel(j+tb[s][1], i+j-1+tb[s][2], L=L)-1)+label;3*(slabel(j+tb[s][3], i+j-1+tb[s][4], L=L)-1)+label])
                end
                for k=1:3, l=1:2
                    a=rot(label)[l]*ones(Int, 3)
                    a[k]=label
                    for s=1:length(tb), i=1:L, j=1:L
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
                    if L > 4
                        push!(tb, [2;-2])
                    end
                    if extra == true
                        if L == 6
                            push!(tb, [0;3], [1;3], [2;3], [3;3], [3;2], [3;1], [3;0], [3;-1], [3;-2])
                        end
                    end
                    for s=1:length(tb), i=1:L, j=1:L
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
                tb = [[0;1;1;1], [0;1;-1;1], [1;0;1;1], [-1;0;-1;1], [1;0;2;0], [0;1;0;2]]
                for s=1:length(tb), k=1:6, i=1:L, j=1:L
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
