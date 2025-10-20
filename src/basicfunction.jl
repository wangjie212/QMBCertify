mutable struct mosek_para
    tol_pfeas::Float64
    tol_dfeas::Float64
    tol_relgap::Float64
    num_threads::Int64
end

mosek_para() = mosek_para(1e-8, 1e-8, 1e-8, 0)

function get_basis(L, label, d; lattice="chain", extra=0, three_type=[1;1])
    basis = Vector{UInt16}[]
    if lattice == "square"
        tb2 = [[1;0], [0;1], [1;1], [1;-1], [2;0], [0;2], [2;1], [1;2], [2;2]]
        lb2 = 9
        if L >= 6
            lb2 = 12
            push!(tb2, [2;-1], [1;-2], [2;-2], [0;3], [1;3], [2;3], [3;3], [3;2], [3;1], [3;0]) # 7
        end
        if L >= 8
            push!(tb2, [3;-1], [3;-2], [3;-3], [2;-3], [1;-3], [0;4], [1;4], [2;4], [3;4], [4;4], [4;3], [4;2], [4;1], [4;0]) # 21
        end
        if L >= 10
            push!(tb2, [4;-1], [4;-2], [4;-3], [4;-4], [3;-4], [2;-4], [1;-4], [0;5], [1;5], [2;5], [3;5], [4;5], [5;5], [5;4], [5;3], [5;2], [5;1], [5;0]) # 39
        end
        tb3 = [[0;1;1;1], [0;1;-1;1], [1;0;1;1], [-1;0;-1;1], [1;0;2;0], [0;1;0;2]]
    end
    if label > 0
        a1 = [[rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1]]]
        a2 = [[label;label;rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1];label;label], [label;rot(label)[1];label;rot(label)[2]], [rot(label)[2];label;rot(label)[1];label], [label;rot(label)[1];rot(label)[2];label], [label;rot(label)[2];rot(label)[1];label], 
        [rot(label)[1];label;label;rot(label)[2]], [rot(label)[2];label;label;rot(label)[1]], [rot(label)[1];label;rot(label)[2];label], [label;rot(label)[2];label;rot(label)[1]], [rot(label)[1];rot(label)[2];label;label], [label;label;rot(label)[2];rot(label)[1]],  
        [rot(label)[1];rot(label)[1];rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1];rot(label)[1];rot(label)[1]], [rot(label)[1];rot(label)[1];rot(label)[2];rot(label)[1]], [rot(label)[1];rot(label)[2];rot(label)[1];rot(label)[1]], 
        [rot(label)[2];rot(label)[2];rot(label)[2];rot(label)[1]], [rot(label)[1];rot(label)[2];rot(label)[2];rot(label)[2]], [rot(label)[2];rot(label)[2];rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1];rot(label)[2];rot(label)[2]]]
        if lattice == "chain"
            for i = 1:L
                push!(basis, [3*(i-1)+label])
            end
            if d > 2
                for i = 1:L
                    push!(basis, sort([3*(i-1)+label;smod(3*(i-1+three_type[1])+label, 3*L);smod(3*(i-1+sum(three_type))+label, 3*L)]))
                end
                for l = 1:2, k = 1:3
                    ind = rot(label)[l]*ones(Int, 3)
                    ind[k] = label
                    for i = 1:L
                        push!(basis, sort([3*(i-1)+ind[1];smod(3*(i-1+three_type[1])+ind[2], 3*L);smod(3*(i-1+sum(three_type))+ind[3], 3*L)]))
                    end
                end
            end
            if d > 1
                for s = 0:extra, k in a1, i = 1:L
                    push!(basis, sort([3*(i-1)+k[1];smod(3*(i+s)+k[2], 3*L)]))
                end
            end
            if d > 3
                for k in a2, i = 1:L
                    push!(basis, sort([3*(i-1)+k[1];smod(3*i+k[2], 3*L);smod(3*(i+1)+k[3], 3*L);smod(3*(i+2)+k[4], 3*L)]))
                end
            end
        else
            for i = 1:L, j = 1:L
                push!(basis, [3*(slabel(i, j, L=L)-1)+label])
            end
            if d > 2
                for s in tb3, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+label;3*(slabel(i+s[1], j+s[2], L=L)-1)+label;3*(slabel(i+s[3], j+s[4], L=L)-1)+label]))
                end
                for k = 1:3, l = 1:2
                    ind = rot(label)[l]*ones(Int, 3)
                    ind[k] = label
                    for s in tb3, i = 1:L, j = 1:L
                        push!(basis, sort([3*(slabel(i, j, L=L)-1)+ind[1];3*(slabel(i+s[1], j+s[2], L=L)-1)+ind[2];3*(slabel(i+s[3], j+s[4], L=L)-1)+ind[3]]))
                    end
                end
            end
            if d > 1
                for k in a1, s in tb2[1:lb2+extra], i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+k[1];3*(slabel(i+s[1], j+s[2], L=L)-1)+k[2]]))
                end
            end
            if d > 3
                for k in a2, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+k[1];3*(slabel(i+1, j, L=L)-1)+k[2];3*(slabel(i, j+1, L=L)-1)+k[3];3*(slabel(i+1, j+1, L=L)-1)+k[4]]))
                end
            end
        end  
    else
        a1 = [[1;2;3], [3;2;1], [1;3;2], [2;3;1], [2;1;3], [3;1;2]]
        a2 = [[1;1;1;1], [2;2;2;2], [3;3;3;3], [1;2;2;1], [2;1;1;2], [1;3;3;1], [3;1;1;3], [3;2;2;3], [2;3;3;2], [1;1;2;2], [2;2;1;1], [1;2;1;2], [2;1;2;1], [1;1;3;3], [3;3;1;1], [1;3;1;3], [3;1;3;1], [3;3;2;2], [2;2;3;3], [3;2;3;2], [2;3;2;3]]
        if lattice == "chain"
            for s = 0:extra, k = 1:3, i = 1:L
                push!(basis, sort([3*(i-1)+k;smod(3*(i+s)+k, 3*L)]))
            end
            if d > 3
                for k in a2, i = 1:L
                    push!(basis, sort([3*(i-1)+k[1];smod(3*i+k[2], 3*L);smod(3*(i+1)+k[3], 3*L);smod(3*(i+2)+k[4], 3*L)]))
                end
            end
            if d > 2
                for k in a1, i = 1:L
                    push!(basis, sort([3*(i-1)+k[1];smod(3*(i-1+three_type[1])+k[2], 3*L);smod(3*(i-1+sum(three_type))+k[3], 3*L)]))
                end
            end
        else
            for s in tb2[1:lb2+extra], k = 1:3, i = 1:L, j = 1:L
                push!(basis, sort([3*(slabel(i, j, L=L)-1)+k;3*(slabel(i+s[1], j+s[2], L=L)-1)+k]))
            end
            if d > 3
                for k in a2, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+k[1];3*(slabel(i+1, j, L=L)-1)+k[2];3*(slabel(i, j+1, L=L)-1)+k[3];3*(slabel(i+1, j+1, L=L)-1)+k[4]]))
                end
            end
            if d > 2
                for s in tb3, k in a1, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+k[1];3*(slabel(i+s[1], j+s[2], L=L)-1)+k[2];3*(slabel(i+s[3], j+s[4], L=L)-1)+k[3]]))
                end
            end
        end
    end
    return basis
end

# binary search in a sorted sequence
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

# reduction to the normal form
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

# reduction to the normal form
function reduce2!(a::Vector{UInt16}; realify=false)
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
            elseif s == [0, 2]
                a[ind] -= UInt16(2)
                coef *= -im
            elseif s == [2, 1] || s == [1, 0]
                a[ind] += UInt16(1)
                coef *= -im
            else
                a[ind] -= UInt16(1)
                coef *= im  
            end
            deleteat!(a, ind+1)
            la -= 1
            flag = 1
        else
            flag = 0
        end
    end
    if realify == true && !isreal(coef)
        coef = imag(coef)
    end
    return a,coef
end

# reduction to the normal form
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

# identify zeros by sign symmetry
function isz(a::Vector{UInt16})
    return any(i->isodd(count(isequal(i), mod.(a,3))), 0:2)
end

# reduction w.r.t symmetries
function reduce4(a::Vector{UInt16}, L; lattice="chain")
    l = length(a)
    if l > 0
        pa = Vector{UInt16}[]
        if lattice == "chain"
            for i = 1:l
                ta = [a[i:end]; a[1:i-1] .+ 3*L] .- 3*(ceil(UInt16, a[i]/3) - 1)
                append!(pa, perm(ta))
                rta = reverse(ta)
                ma = 3*(ceil(UInt16, ta[end]/3) .- ceil.(UInt16, rta/3)) + smod.(rta, 3)
                append!(pa, perm(ma))
            end
        else
            factor = [[1;1], [-1;1], [1;-1], [-1;-1]]
            loc = location.(ceil.(UInt16, a/3))
            for i = 1:l
                temp = zeros(UInt16, l)
                for k = 1:4
                    for j = 1:l
                        p = slabel(factor[k][1]*(loc[j][1]-loc[i][1])+1, factor[k][2]*(loc[j][2]-loc[i][2])+1, L=L)
                        temp[j] = 3*p + a[j] - 3*ceil(UInt16, a[j]/3)
                    end
                    append!(pa, perm(sort(temp)))
                    for j = 1:l
                        p = slabel(factor[k][1]*(loc[j][2]-loc[i][2])+1, factor[k][2]*(loc[j][1]-loc[i][1])+1, L=L)
                        temp[j] = 3*p + a[j] - 3*ceil(UInt16, a[j]/3)
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

# implement all reductions
function reduce!(a::Vector{UInt16}; L=0, lattice="chain", realify=false)
    reduce1!(a)
    reduce3!(a)
    a,coef = reduce2!(a, realify=realify)
    reduce3!(a)
    if isz(a)
        coef = 0
    else
        a = reduce4(a, L, lattice=lattice)
    end
    return a,coef
end

function perm(a)
    ra = smod.(a, 3)
    sym = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2], [3;2;1]]
    return [UInt16.(3*(ceil.(Int, a./3).-1) .+ sym[i][ra]) for i=1:6]
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

function smod(i, s)
    r = mod(i, s)
    return r == 0 ? s : r
end

# compute the eigenvalues of a (symmetry, real) circulant matrix
function eigen_circmat(supp, coe, L; symmetry=false, real_matrix=false)
    ne = real_matrix == true ? Int(L/2)+1 : L
    seig = [Vector{UInt16}[] for i = 1:ne]
    if symmetry == false
        ceig = [ComplexF64[] for i = 1:ne]
        for i = 1:ne
            for j = 1:L, (s,c) in enumerate(coe[j])
                if c != 0
                    push!(seig[i], supp[j][s])
                    push!(ceig[i], c*(cos(2*pi*(i-1)*(j-1)/L) + sin(2*pi*(i-1)*(j-1)/L)*im))
                end
            end
            if !isempty(ceig[i])
                seig[i], ceig[i] = resort(seig[i], ceig[i])
            end
        end
    else
        ceig = [Float64[] for i = 1:ne]
        for i = 1:ne
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
            if !isempty(ceig[i])
                seig[i], ceig[i] = resort(seig[i], ceig[i])
            end
        end
    end
    return seig, ceig
end

function add_SU2_equality!(model, tsupp, ltsupp, cons; L=0, lattice="chain")
    ind = findall(item->length(item) == 4 && all(smod.(item, 3) .== 1), tsupp)
    for item in tsupp[ind]
        fr = @variable(model)
        Locb = bfind(tsupp, ltsupp, item)
        add_to_expression!(cons[Locb], fr)
        for i = 2:4
            a = copy(item)
            a[1] += 1
            a[i] += 1
            a = reduce!(a, L=L, lattice=lattice)[1]
            Locb = bfind(tsupp, ltsupp, a)
            add_to_expression!(cons[Locb], -1, fr)
        end
    end
    ind = findall(item->length(item) == 6 && sum(smod.(item, 3) .== 1) == 4, tsupp)
    for item in tsupp[ind]
        ino = Vector(1:6)[smod.(item, 3) .== 1]
        fr = @variable(model)
        Locb = bfind(tsupp, ltsupp, item)
        add_to_expression!(cons[Locb], fr)
        for i = 2:4
            a = copy(item)
            a[ino[1]] += 2
            a[ino[i]] += 2
            a = reduce!(a, L=L, lattice=lattice)[1]
            Locb = bfind(tsupp, ltsupp, a)
            add_to_expression!(cons[Locb], -1, fr)
        end
    end
    ind = findall(item->length(item) == 6 && all(smod.(item, 3) .== 1), tsupp)
    for item in tsupp[ind]
        fr = @variable(model)
        Locb = bfind(tsupp, ltsupp, item)
        add_to_expression!(cons[Locb], fr)
        for i = 2:6
            ino = [Vector(2:i-1); Vector(i+1:6)]
            for j = 2:4
                ine = [Vector(2:j-1); Vector(j+1:4)]
                a = copy(item)
                a[ino[1]] += 1
                a[ino[j]] += 1
                a[ino[ine[1]]] += 2
                a[ino[ine[2]]] += 2
                a = reduce!(a, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, ltsupp, a)
                add_to_expression!(cons[Locb], -1, fr)
            end
        end
    end
    ind = findall(item->length(item) == 8 && sum(smod.(item, 3) .== 1) == 6, tsupp)
    for item in tsupp[ind]
        ino = Vector(1:8)[smod.(item, 3) .== 1]
        for i = 1:6
            fr = @variable(model)
            Locb = bfind(tsupp, ltsupp, item)
            add_to_expression!(cons[Locb], fr)
            for j in [Vector(1:i-1); Vector(i+1:6)]
                a = copy(item)
                a[ino[i]] += 2
                a[ino[j]] += 2
                a = reduce!(a, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, ltsupp, a)
                add_to_expression!(cons[Locb], -1, fr)
            end
        end
    end
    ind = findall(item->length(item) == 8 && sum(smod.(item, 3) .== 1) == 4 && sum(smod.(item, 3) .== 2) == 4, tsupp)
    for item in tsupp[ind]
        ino = Vector(1:8)[smod.(item, 3) .== 2]
        fr = @variable(model)
        Locb = bfind(tsupp, ltsupp, item)
        add_to_expression!(cons[Locb], fr)
        for i in 2:4
            a = copy(item)
            a[ino[1]] += 1
            a[ino[i]] += 1
            a = reduce!(a, L=L, lattice=lattice)[1]
            Locb = bfind(tsupp, ltsupp, a)
            add_to_expression!(cons[Locb], -1, fr)
        end
    end
end
