mutable struct mosek_para
    tol_pfeas::Float64
    tol_dfeas::Float64
    tol_relgap::Float64
end

mosek_para() = mosek_para(1e-8, 1e-8, 1e-8)

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
        ind = findfirst(x -> a[x] != a[x+1] && ceil(Int, a[x]/3) == ceil(Int, a[x+1]/3), 1:la-1)
        if ind !== nothing
            s = mod.(a[ind:ind+1], 3)
            if s == [1, 2]
                a[ind] += UInt16(2)
                coef *= im
            elseif s == [2, 1]
                a[ind] += UInt16(1)
                coef *= -im
            elseif s == [1, 0]
                a[ind] += UInt16(1)
                coef *= -im
            elseif s == [0, 1]
                a[ind] -= UInt16(1)
                coef *= im
            elseif s == [2, 0]
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
                    for j = 1:l
                        p = slabel(factor[k][1]*(loc[j][2]-loc[i][2])+1, factor[k][2]*(loc[j][1]-loc[i][1])+1, L=L)
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

function perm(a)
    ra = smod.(a, 3)
    sym = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2], [3;2;1]]
    return [UInt16.(3*(ceil.(Int, a./3).-1) .+ sym[i][ra]) for i=1:6]
end

function reduce!(a::Vector{UInt16}; L=0, lattice="chain")
    reduce1!(a)
    reduce3!(a)
    a,coef = reduce2!(a)
    reduce3!(a)
    if isz(a)
        coef = 0
    else
        a = reduce4(a, L, lattice=lattice)
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
    if label > 0
        basis = Vector{UInt16}[]
        if lattice == "chain"
            for i = 1:L
                push!(basis, [3*(i-1)+label])
            end
        else
            for i = 1:L, j = 1:L
                push!(basis, [3*(slabel(i, j, L=L)-1)+label])
            end
        end
        if d > 1
            a = [[rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1]]]
            for k = 1:2
                if lattice == "chain"
                    for s = 0:extra, i = 1:L
                        push!(basis, sort([3*(i-1)+a[k][1];smod(3*(i+s)+a[k][2], 3*L)]))
                    end
                else
                    tb = [[1;0], [0;1], [1;1], [1;-1], [2;0], [0;2], [2;1], [1;2], [2;2], [2;-1]]
                    if L > 4
                        push!(tb, [1;-2], [2;-2])
                    end
                    if extra == true
                        if L >= 6
                            push!(tb, [0;3], [1;3], [2;3], [3;3], [3;2], [3;1], [3;0], [3;-1], [3;-2])
                        end
                        if L >= 8
                            push!(tb, [3;-3], [2;-3], [1;-3])
                        end
                    end
                    for s in tb, i = 1:L, j = 1:L
                        push!(basis, sort([3*(slabel(i, j, L=L)-1)+a[k][1];3*(slabel(i+s[1], j+s[2], L=L)-1)+a[k][2]]))
                    end
                end
            end
        end
        if d > 2
            if lattice == "chain"
                for i = 1:L
                    push!(basis, sort([3*(i-1)+label;smod(3*(i-1+three_type[1])+label, 3*L);smod(3*(i-1+sum(three_type))+label, 3*L)]))
                end
                for k = 1:3, l = 1:2
                    a = rot(label)[l]*ones(Int, 3)
                    a[k] = label
                    for i = 1:L
                        push!(basis, sort([3*(i-1)+a[1];smod(3*(i-1+three_type[1])+a[2], 3*L);smod(3*(i-1+sum(three_type))+a[3], 3*L)]))
                    end
                end
            else
                tb = [[0;1;1;1], [0;1;-1;1], [1;0;1;1], [-1;0;-1;1], [1;0;2;0], [0;1;0;2]]
                for s in tb, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+label;3*(slabel(i+s[1], j+s[2], L=L)-1)+label;3*(slabel(i+s[3], j+s[4], L=L)-1)+label]))
                end
                for k = 1:3, l = 1:2
                    a = rot(label)[l]*ones(Int, 3)
                    a[k] = label
                    for s in tb, i = 1:L, j = 1:L
                        push!(basis, sort([3*(slabel(i, j, L=L)-1)+a[1];3*(slabel(i+s[1], j+s[2], L=L)-1)+a[2];3*(slabel(i+s[3], j+s[4], L=L)-1)+a[3]]))
                    end
                end
            end
        end
        if d > 3
            a = [[label;label;rot(label)[1];rot(label)[2]], [label;rot(label)[1];label;rot(label)[2]], [label;rot(label)[1];rot(label)[2];label], [rot(label)[1];label;label;rot(label)[2]], [rot(label)[1];label;rot(label)[2];label], [rot(label)[1];rot(label)[2];label;label],
            [label;label;rot(label)[2];rot(label)[1]], [label;rot(label)[2];label;rot(label)[1]], [label;rot(label)[2];rot(label)[1];label], [rot(label)[2];label;label;rot(label)[1]], [rot(label)[2];label;rot(label)[1];label], [rot(label)[2];rot(label)[1];label;label],
            [rot(label)[1];rot(label)[1];rot(label)[1];rot(label)[2]], [rot(label)[1];rot(label)[1];rot(label)[2];rot(label)[1]], [rot(label)[1];rot(label)[2];rot(label)[1];rot(label)[1]], [rot(label)[2];rot(label)[1];rot(label)[1];rot(label)[1]],
            [rot(label)[2];rot(label)[2];rot(label)[2];rot(label)[1]], [rot(label)[2];rot(label)[2];rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1];rot(label)[2];rot(label)[2]], [rot(label)[1];rot(label)[2];rot(label)[2];rot(label)[2]]]
            if lattice == "chain"
                for k = 1:20, i = 1:L
                    push!(basis, sort([3*(i-1)+a[k][1];smod(3*i+a[k][2], 3*L);smod(3*(i+1)+a[k][3], 3*L);smod(3*(i+2)+a[k][4], 3*L)]))
                end
            else
                for k = 1:20, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+a[k][1];3*(slabel(i+1, j, L=L)-1)+a[k][2];3*(slabel(i, j+1, L=L)-1)+a[k][3];3*(slabel(i+1, j+1, L=L)-1)+a[k][4]]))
                end
            end
        end
    else
        basis = Vector{UInt16}[]
        if d > 1
            for k = 1:3
                if lattice == "chain"
                    for s = 0:extra, i = 1:L
                        push!(basis, sort([3*(i-1)+k;smod(3*(i+s)+k, 3*L)]))
                    end
                else
                    tb = [[1;0], [0;1], [1;1], [1;-1], [2;0], [0;2], [2;1], [1;2], [2;2], [2;-1]]
                    if L > 4
                        push!(tb, [1;-2], [2;-2])
                    end
                    if extra == true
                        if L >= 6
                            push!(tb, [0;3], [1;3], [2;3], [3;3], [3;2], [3;1], [3;0], [3;-1], [3;-2])
                        end
                        if L >= 8
                            push!(tb, [3;-3], [2;-3], [1;-3])
                        end
                    end
                    for s in tb, i = 1:L, j = 1:L
                        push!(basis, sort([3*(slabel(i, j, L=L)-1)+k;3*(slabel(i+s[1], j+s[2], L=L)-1)+k]))
                    end
                end
            end
        end
        if d > 2
            a = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2], [3;2;1]]
            if lattice == "chain"
                for k = 1:6, i = 1:L
                    push!(basis, sort([3*(i-1)+a[k][1];smod(3*(i-1+three_type[1])+a[k][2], 3*L);smod(3*(i-1+sum(three_type))+a[k][3], 3*L)]))
                end
            else
                tb = [[0;1;1;1], [0;1;-1;1], [1;0;1;1], [-1;0;-1;1], [1;0;2;0], [0;1;0;2]]
                for s in tb, k = 1:6, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+a[k][1];3*(slabel(i+s[1], j+s[2], L=L)-1)+a[k][2];3*(slabel(i+s[3], j+s[4], L=L)-1)+a[k][3]]))
                end
            end
        end
        if d > 3
            a = [[1;1;1;1], [2;2;2;2], [3;3;3;3], [1;1;2;2], [1;2;1;2], [1;2;2;1], [2;1;1;2], [2;1;2;1], [2;2;1;1],
            [1;1;3;3], [1;3;1;3], [1;3;3;1], [3;1;1;3], [3;1;3;1], [3;3;1;1], [3;3;2;2], [3;2;3;2], [3;2;2;3], [2;3;3;2], [2;3;2;3], [2;2;3;3]]
            if lattice == "chain"
                for k = 1:21, i = 1:L
                    push!(basis, sort([3*(i-1)+a[k][1];smod(3*i+a[k][2], 3*L);smod(3*(i+1)+a[k][3], 3*L);smod(3*(i+2)+a[k][4], 3*L)]))
                end
            else
                for k = 1:21, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+a[k][1];3*(slabel(i+1, j, L=L)-1)+a[k][2];3*(slabel(i, j+1, L=L)-1)+a[k][3];3*(slabel(i+1, j+1, L=L)-1)+a[k][4]]))
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

function eigen_circmat(supp, coe, L; symmetry=false)
    seig = [Vector{UInt16}[] for i = 1:L]
    if symmetry == false
        ceig = [ComplexF64[] for i = 1:L]
        for i = 1:L
            for j = 1:L, (s,c) in enumerate(coe[j])
                if c != 0
                    push!(seig[i], supp[j][s])
                    push!(ceig[i], real(c)*cos(2*pi*(i-1)*(j-1)/L)-imag(c)*sin(2*pi*(i-1)*(j-1)/L)+(real(c)*sin(2*pi*(i-1)*(j-1)/L)+imag(c)*cos(2*pi*(i-1)*(j-1)/L))*im)
                end
            end
            seig[i], ceig[i] = resort(seig[i], ceig[i])
        end
    else
        ceig = [Float64[] for i = 1:L]
        for i = 1:L 
            for (s,c) in enumerate(coe[1])
                if c != 0
                    push!(seig[i], supp[1][s])
                    push!(ceig[i], c)
                end
            end
            for j = 2:Int(L/2), (s,c) in enumerate(coe[j])
                if c != 0
                    push!(seig[i], supp[j][s])
                    push!(ceig[i], 2*c*cos(2*pi*(i-1)*(j-1)/L))
                end
            end
            for (s,c) in enumerate(coe[end])
                if c != 0
                    push!(seig[i], supp[end][s])
                    push!(ceig[i], c*(-1)^(i-1))
                end
            end
            seig[i], ceig[i] = resort(seig[i], ceig[i])
        end
    end
    return seig, ceig
end
