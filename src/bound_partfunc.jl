function PFB(supp::Vector{Vector{Int}}, coe::Vector{Float64}, beta, L::Int, d::Int; QUIET=false)
    println("*********************************** QMBCertify ***********************************")
    println("QMBCertify is launching...")
    if QUIET == false
        println("Generating the monomial basis...")
    end
    time = @elapsed begin
    basis = [get_pfbasis(L, i, d) for i = 1:4]
    basis_loct = [get_pfbasis_loct(L, i) for i = 1:2]
    basis_loc1 = [get_pfbasis_loc1(L, i) for i = 1:2]
    coe1 = Vector{Vector{Vector{Int8}}}(undef, 4)
    bi1 = Vector{Vector{Vector{Vector{UInt16}}}}(undef, 4)
    coe2 = Vector{Vector{Vector{Int8}}}(undef, 4)
    bi2 = Vector{Vector{Vector{Vector{UInt16}}}}(undef, 4)
    tsupp = Vector{UInt16}[]
    for i = 1:4
        p = length(basis[i][1])
        k = Int(length(basis[i][2])/L)
        for j = 1:p, s = j:p
            bi = sadd(involution(basis[i][1][j]), basis[i][1][s], realify=true)[1]
            push!(tsupp, reduce_pf(bi, L))
        end
        for j = 1:p, s = 1:k
            bi = sadd(involution(basis[i][1][j]), basis[i][2][L*(s-1)+1], realify=true)[1]
            if !iszero(bi)
                push!(tsupp, reduce_pf(bi, L))
            end
        end
        # processing diagonal blocks
        coe1[i] = Vector{Vector{Int8}}(undef, k)
        bi1[i] = Vector{Vector{Vector{UInt16}}}(undef, k)
        for j = 1:k
            coe1[i][j] = Vector{Int8}(undef, Int(L/2)+1)
            bi1[i][j] = Vector{Vector{UInt16}}(undef, Int(L/2)+1)
            for r = 1:Int(L/2)+1
                bi,coe1[i][j][r] = sadd(involution(basis[i][2][L*(j-1)+1]), basis[i][2][L*(j-1)+r], realify=true)
                if iszero(bi)
                    coe1[i][j][r] = 0
                else
                    bi1[i][j][r] = reduce_pf(bi, L)
                    push!(tsupp, bi1[i][j][r])
                end
            end
        end
        # processing non-diagonal blocks
        coe2[i] = Vector{Vector{Int8}}(undef, Int(k*(k-1)/2))
        bi2[i] = Vector{Vector{Vector{UInt16}}}(undef, Int(k*(k-1)/2))
        for j1 = 1:k-1, j2 = j1+1:k
            j = Int((2*k-j1)*(j1-1)/2) + j2 - j1
            coe2[i][j] = Vector{Int8}(undef, L)
            bi2[i][j] = Vector{Vector{UInt16}}(undef, L)
            for r = 1:L
                bi,coe2[i][j][r] = sadd(involution(basis[i][2][L*(j1-1)+1]), basis[i][2][L*(j2-1)+r], realify=true)
                if iszero(bi)
                    coe2[i][j][r] = 0
                else
                    bi2[i][j][r] = reduce_pf(bi, L)
                    push!(tsupp, bi2[i][j][r])
                end
            end
        end
    end
    # basis0 = [get_pfbasis0(L, i) for i = 1:4]
    # for i = 1:4, j = 1:length(basis0[i]), k = j:length(basis0[i])
    #     bi = sadd(involution(basis0[i][j]), basis0[i][k])[1]
    #     if !iszero(bi)
    #         push!(tsupp, reduce_pf(bi, L))
    #     end
    # end
    sort!(tsupp)
    unique!(tsupp)
    end
    if QUIET == false
        println("Obtained the monomial basis in $time seconds.")
        println("Starting block-diagonalization...")
    end
    time = @elapsed begin
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:length(tsupp)]
    mb = 0
    for i = 1:4
        p = length(basis[i][1])
        k = Int(length(basis[i][2])/L)
        mb = maximum([2k, k+p, mb])
        for l = 1:Int(L/2) + 1
            if l == 1
                pos = @variable(model, [1:k+p, 1:k+p], PSD)
                for j = 1:p, s = j:p
                    bi = sadd(involution(basis[i][1][j]), basis[i][1][s], realify=true)[1]
                    Locb = bfind(tsupp, reduce_pf(bi, L))
                    if j == s
                        @inbounds add_to_expression!(cons[Locb], pos[j,s])
                    else
                        @inbounds add_to_expression!(cons[Locb], 2, pos[j,s])
                    end
                end
                for j = 1:p, s = 1:k
                    bi,coe = sadd(involution(basis[i][1][j]), basis[i][2][L*(s-1)+1], realify=true)
                    if !iszero(bi)
                        Locb = bfind(tsupp, reduce_pf(bi, L))
                        @inbounds add_to_expression!(cons[Locb], 2*coe*sqrt(L), pos[j,s+p])
                    end
                end
            elseif l == Int(L/2) + 1
                pos = @variable(model, [1:k, 1:k], PSD)
            else
                pos = @variable(model, [1:2k, 1:2k], PSD)
            end
            for s = 1:k
                if l == 1
                    pp = pos[s+p, s+p]
                elseif l == Int(L/2) + 1
                    pp = pos[s, s]
                else
                    pp = pos[s, s] + pos[s+k, s+k]
                end
                if coe1[i][s][1] != 0
                    Locb = bfind(tsupp, bi1[i][s][1])
                    @inbounds add_to_expression!(cons[Locb], coe1[i][s][1], pp)
                end
                if coe1[i][s][end] != 0
                    Locb = bfind(tsupp, bi1[i][s][end])
                    @inbounds add_to_expression!(cons[Locb], coe1[i][s][end]*(-1)^(l-1), pp)
                end
                for r = 2:Int(L/2)
                    if coe1[i][s][r] != 0
                        Locb = bfind(tsupp, bi1[i][s][r])
                        @inbounds add_to_expression!(cons[Locb], 2*coe1[i][s][r]*cos(2*pi*(r-1)*(l-1)/L), pp)
                    end
                end
            end
            for j1 = 1:k-1, j2 = j1+1:k
                if l == 1
                    pp1 = pos[j1+p, j2+p]
                elseif l == Int(L/2) + 1
                    pp1 = pos[j1, j2]
                else
                    pp1 = pos[j1, j2] + pos[j1+k, j2+k]
                    pp2 = pos[j1+k, j2] - pos[j2+k, j1]
                end
                j = Int((2*k-j1)*(j1-1)/2) + j2 - j1
                for r = 1:L
                    if coe2[i][j][r] != 0
                        Locb = bfind(tsupp, bi2[i][j][r])
                        @inbounds add_to_expression!(cons[Locb], 2*coe2[i][j][r]*cos(2*pi*(r-1)*(l-1)/L), pp1)
                        if l != 1 && l != Int(L/2) + 1
                            @inbounds add_to_expression!(cons[Locb], -2*coe2[i][j][r]*sin(2*pi*(r-1)*(l-1)/L), pp2)
                        end
                    end
                end
            end
        end
    end
    # for i = 1:4
    #     mb = max(length(basis0[i]), mb)
    #     pos = @variable(model, [1:length(basis0[i]), 1:length(basis0[i])], PSD)
    #     for j = 1:length(basis0[i]), k = j:length(basis0[i])
    #         bi,coef = sadd(involution(basis0[i][j]), basis0[i][k], realify=true)
    #         if !iszero(bi)
    #             Locb = bfind(tsupp, reduce_pf(bi, L))
    #             if j == k
    #                 @inbounds add_to_expression!(cons[Locb], coef, pos[j,k])
    #             else
    #                 @inbounds add_to_expression!(cons[Locb], 2*coef, pos[j,k])
    #             end
    #         end
    #     end
    # end
    pos = @variable(model, [1:length(basis_loct[1]), 1:length(basis_loct[1])], PSD)
    for i = 1:length(basis_loct[1]), j = i:length(basis_loct[1])
        bi,coef = sadd(basis_loct[1][i], basis_loct[1][j], realify=true)
        Locb1 = bfind(tsupp, reduce_pf(UInt16[1;bi], L))
        Locb2 = bfind(tsupp, reduce_pf(UInt16[1;1;bi], L))
        if i == j
            @inbounds add_to_expression!(cons[Locb1], coef, pos[i,j])
            @inbounds add_to_expression!(cons[Locb2], -coef, pos[i,j])
        else
            @inbounds add_to_expression!(cons[Locb1], 2*coef, pos[i,j])
            @inbounds add_to_expression!(cons[Locb2], -2*coef, pos[i,j])
        end
    end
    # pos = @variable(model, [1:length(basis_loct[1]), 1:length(basis_loct[1])], PSD)
    # for i = 1:length(basis_loct[1]), j = i:length(basis_loct[1])
    #     bi,coef = sadd(basis_loct[1][i], UInt16[2;basis_loct[1][j]], realify=true)
    #     Locb = bfind(tsupp, reduce_pf(bi, L))
    #     if i == j
    #         @inbounds add_to_expression!(cons[Locb], coef, pos[i,j])
    #     else
    #         @inbounds add_to_expression!(cons[Locb], 2*coef, pos[i,j])
    #     end
    # end
    pos = @variable(model, [1:length(basis_loct[1]), 1:length(basis_loct[1])], PSD)
    for i = 1:length(basis_loct[1]), j = i:length(basis_loct[1])
        bi,coef1 = sadd(basis_loct[1][i], UInt16[2;basis_loct[1][j]], realify=true)
        Locb1 = bfind(tsupp, reduce_pf(bi, L))
        bi,coef2 = sadd(basis_loct[1][i], UInt16[2;2;basis_loct[1][j]], realify=true)
        Locb2 = bfind(tsupp, reduce_pf(bi, L))
        if i == j
            @inbounds add_to_expression!(cons[Locb1], 10*coef1, pos[i,j])
            @inbounds add_to_expression!(cons[Locb2], -coef2, pos[i,j])
        else
            @inbounds add_to_expression!(cons[Locb1], 20*coef1, pos[i,j])
            @inbounds add_to_expression!(cons[Locb2], -2*coef2, pos[i,j])
        end
    end
    pos = @variable(model, [1:Int(L/2) + 1], lower_bound=0)
    for l = 1:Int(L/2) + 1, r = 1:L
        bi,coef = sadd(basis_loct[2][1], UInt16[2;basis_loct[2][r]], realify=true)
        Locb = bfind(tsupp, reduce_pf(bi, L))
        @inbounds add_to_expression!(cons[Locb], coef, pos[l]*cos(2*pi*(r-1)*(l-1)/L))
    end
    pos = @variable(model, [1:length(basis_loc1[1]), 1:length(basis_loc1[1])], PSD)
    for i = 1:length(basis_loc1[1]), j = i:length(basis_loc1[1])
        bi,coef = sadd(basis_loc1[1][i], UInt16[3;basis_loc1[1][j]], realify=true)
        Locb = bfind(tsupp, reduce_pf(bi, L))
        if i == j
            @inbounds add_to_expression!(cons[Locb], coef, pos[i,j])
        else
            @inbounds add_to_expression!(cons[Locb], 2*coef, pos[i,j])
        end
    end
    # pos = @variable(model, [1:length(basis_loc1[1]), 1:length(basis_loc1[1])], PSD)
    # for i = 1:length(basis_loc1[1]), j = i:length(basis_loc1[1])
    #     bi,coef1 = sadd(basis_loc1[1][i], basis_loc1[1][j], realify=true)
    #     Locb1 = bfind(tsupp, reduce_pf(bi, L))
    #     bi,coef2 = sadd(basis_loc1[1][i], UInt16[3;basis_loc1[1][j]], realify=true)
    #     Locb2 = bfind(tsupp, reduce_pf(bi, L))
    #     if i == j
    #         @inbounds add_to_expression!(cons[Locb1], 10*coef1, pos[i,j])
    #         @inbounds add_to_expression!(cons[Locb2], -coef2, pos[i,j])
    #     else
    #         @inbounds add_to_expression!(cons[Locb1], 20*coef1, pos[i,j])
    #         @inbounds add_to_expression!(cons[Locb2], -2*coef2, pos[i,j])
    #     end
    # end
    pos = @variable(model, [1:Int(L/2) + 1], lower_bound=0)
    for l = 1:Int(L/2) + 1, r = 1:L
        bi,coef = sadd(basis_loc1[2][1], UInt16[3;basis_loc1[2][r]], realify=true)
        Locb = bfind(tsupp, reduce_pf(bi, L))
        @inbounds add_to_expression!(cons[Locb], coef, pos[l]*cos(2*pi*(r-1)*(l-1)/L))
    end
    free = @variable(model, [1:2d])
    for i = 1:2d
        Locb = bfind(tsupp, ones(UInt16, i))
        @inbounds add_to_expression!(cons[Locb], free[i])
        @inbounds add_to_expression!(cons[1], -1/(i+1), free[i])
    end
    free = @variable(model, [1:3])
    Locb = bfind(tsupp, [2;4;7])
    @inbounds add_to_expression!(cons[Locb], -3L*beta, free[1])
    Locb = bfind(tsupp, [3])
    @inbounds add_to_expression!(cons[Locb], -1, free[1])
    @inbounds add_to_expression!(cons[1], 1, free[1])
    Locb = bfind(tsupp, [2;2;4;7])
    @inbounds add_to_expression!(cons[Locb], -6L*beta, free[2])
    Locb = bfind(tsupp, [3;3])
    @inbounds add_to_expression!(cons[Locb], -1, free[2])
    @inbounds add_to_expression!(cons[1], 1, free[2])
    Locb = bfind(tsupp, [2])
    @inbounds add_to_expression!(cons[Locb], 1, free[3])
    Locb = bfind(tsupp, [1;2;4;7])
    @inbounds add_to_expression!(cons[Locb], -3L*beta, free[3])
    Locb = bfind(tsupp, [3])
    @inbounds add_to_expression!(cons[Locb], -1, free[3])
    end
    if QUIET == false
        println("Finished block-diagonalization in $time seconds.")
        println("SDP size: n = $mb, m = $(length(tsupp))")
    end
    cons[bfind(tsupp, [3])] -= 1
    @variable(model, lower)
    @objective(model, Max, lower)
    cons[1] += lower
    @constraint(model, cons .== 0)
    if QUIET == false
        println("Solving the SDP...")
    end
    time = @elapsed begin
    optimize!(model)
    end
    if QUIET == false
        println("SDP solving time: $time seconds.")
    end
    status = termination_status(model)
    objv = objective_value(model)
    if status != MOI.OPTIMAL
       println("termination status: $status")
       status = primal_status(model)
       println("solution status: $status")
    end
    println("optimum = $objv")
    return objv,tsupp
end

function involution(a::Vector{UInt16})
    ind = findall(x -> x == 2 || x == 3, a)
    if !isempty(ind)
        b = copy(a)
        b[ind] = reverse(b[ind])
        return b
    end
    return a
end

function sadd(a::Vector{UInt16}, b::Vector{UInt16}; realify=false)
    s = count(x -> x == 1, a) + count(x -> x == 1, b)
    a = a[findall(x -> x != 1, a)]
    b = b[findall(x -> x != 1, b)]
    word,coef = [a; b],1
    if !isempty(a) && !isempty(b)
        ind1 = findlast(x -> x == 2 || x == 3, a)
        ind1 = ind1 === nothing ? 0 : ind1
        if ind1 < length(a)
            ind2 = findfirst(x -> x == 2 || x == 3, b)
            ind2 = ind2 === nothing ? length(b)+1 : ind2
            if ind2 > 1
                temp = [a[ind1+1:end]; b[1:ind2-1]]
                temp,coef = reduce2!(reduce3!(reduce1!(temp)))
                reduce3!(temp)
                word = [a[1:ind1]; temp; b[ind2:end]]
            end
        end
        if length(word) > 1 && word[1] > 3 && word[end] > 3 && ceil(Int, word[1]/3) == ceil(Int, word[end]/3)
            if word[1] == word[end]
                word = word[2:end-1]
            else
                w,c = reduce2!([word[end];word[1]])
                coef *= c
                word = [word[2:end-1];w]
            end
        end
        if realify == true && !isreal(coef)
            coef = imag(coef)
        end
    end
    ind = findfirst(x -> x == 2 || x == 3, word)
    if ind !== nothing
        word = [word[ind:end]; word[1:ind-1]]
    end
    return [ones(UInt16, s); word],coef
end

function iszero(a::Vector{UInt16})
    return !isempty(a) && !all(x -> x == 1, a) && !any(x -> x == 2 || x == 3, a)
end

function reduce_pf(a::Vector{UInt16}, L)
    s = count(x -> x == 1, a)
    ind = findall(x -> x != 1, a)
    pa = unique(vcat(cyclic_extension.(reverse_extension(a[ind]))...))
    if findfirst(x -> x > 3, a) !== nothing
        c = Vector{UInt16}[]
        for b in pa
            ind = findall(x -> x > 3, b)
            pb = symmetry_extension(b[ind], L)
            for d in pb
                temp = copy(b)
                temp[ind] = d
                push!(c, [ones(UInt16, s); temp])
            end
        end
        return findmin(c)[1]
    else
        return [ones(UInt16, s); findmin(pa)[1]]
    end
end

function cyclic_extension(a::Vector{UInt16})
    ind = findall(x -> x == 2 || x == 3, a)
    if length(ind) > 1
        return unique([[a[i:end]; a[1:i-1]] for i in ind])
    else
        return [a]
    end
end

function reverse_extension(a::Vector{UInt16})
    return unique([a, involution(a)])
end

function permutation_extension(a::Vector{UInt16})
    ind = findall(x -> x > 3, a)
    if !isempty(ind)
        b = unique(perm(a[ind]))
        c = [copy(a) for i = 1:length(b)]
        for i = 1:length(b)
            c[i][ind] = b[i]
        end
        return c        
    end
    return [a]
end

function translation_extension(a::Vector{UInt16}, L)
    ind = findall(x -> x > 3, a)
    if !isempty(ind)
        b = a[ind]
        sites = unique(ceil.(UInt16, b/3))
        c = Vector{Vector{UInt16}}(undef, length(sites))
        for (i, s) in enumerate(sites)
            c[i] = copy(a)
            temp = b .- 3*(s - 2)
            for j = 1:length(temp)
                if temp[j] <= 3
                    temp[j] += 3L
                end
            end
            c[i][ind] = sort(temp)
        end
        return unique(c)
    end
    return [a]
end

function mirror_extension(a::Vector{UInt16})
    ind = findall(x -> x > 3, a)
    if !isempty(ind)
        b = copy(a)
        c = reverse(b[ind])
        sites = ceil.(UInt16, c/3)
        b[ind] = 3*(maximum(sites) + 1 .- sites) + smod.(c, 3)
        return unique([a, b])        
    end
    return [a]
end

function symmetry_extension(a::Vector{UInt16}, L)
    return unique(vcat(permutation_extension.(vcat(translation_extension.(mirror_extension(a), L)...))...))
end

function get_pfbasis(L, label, d)
    basis = Vector{Vector{Vector{UInt16}}}(undef, 2)
    if label == 1
        basis[1] = Vector{UInt16}[[], [1], [2]]
        basis[2] = Vector{UInt16}[]
        if d > 1
            push!(basis[1], [1;1], [1;2], [2;2])
            append!(basis[2], [sort([3i+1;smod(3i+1, 3L)+3]) for i = 1:L])
        end
        if d > 2
            push!(basis[1], [1;1;1], [1;1;2], [1;2;2], [2;2;2])
        end
    elseif label == 2
        basis[1] = Vector{UInt16}[]
        basis[2] = Vector{UInt16}[[3i+1] for i = 1:L]
        if d > 1
            append!(basis[2], [[[1;3i+1] for i = 1:L]; [[3i+1;2] for i = 1:L]; [[2;3i+1] for i = 1:L]])
        end
        if d > 2
            append!(basis[2], [[[1;2;3i+1] for i = 1:L]; [[1;1;3i+1] for i = 1:L]; [[2;2;3i+1] for i = 1:L]; [[2;3i+1;2] for i = 1:L]; [[3i+1;2;2] for i = 1:L]; [[1;3i+1;2] for i = 1:L]; 
            [[1;3i+2;3i+6] for i = 1:L]; [[1;3i+3;3i+5] for i = 1:L]; [[2;3i+2;3i+6] for i = 1:L]; [[3i+2;2;3i+6] for i = 1:L]; [[3i+2;3i+6;2] for i = 1:L]; [[2;3i+3;3i+5] for i = 1:L]; 
            [[3i+3;2;3i+5] for i = 1:L]; [[3i+3;3i+5;2] for i = 1:L]])
        end
        # if d > 1
        #     append!(basis[2], [[[3i+2;3i+6] for i = 1:L]; [[3i+3;3i+5] for i = 1:L]])
        # end
    elseif label == 3
        basis[1] = Vector{UInt16}[[], [3]]
        basis[2] = Vector{UInt16}[]
        if d > 1
            push!(basis[1], [3;3])
            append!(basis[2], [sort([3i+1;smod(3i+1, 3L)+3]) for i = 1:L])
        end
        if d > 2
            push!(basis[1], [3;3;3])
        end
    else
        basis[1] = Vector{UInt16}[]
        basis[2] = Vector{UInt16}[[3i+1] for i = 1:L]
        if d > 1
            append!(basis[2], [[[3i+1;3] for i = 1:L]; [[3;3i+1] for i = 1:L]])
        end
        if d > 2
            append!(basis[2], [[[3;3;3i+1] for i = 1:L]; [[3;3i+1;3] for i = 1:L]; [[3i+1;3;3] for i = 1:L]; [[3;3i+2;3i+6] for i = 1:L]; [[3i+2;3;3i+6] for i = 1:L]; 
            [[3i+2;3i+6;3] for i = 1:L]; [[3;3i+3;3i+5] for i = 1:L]; [[3i+3;3;3i+5] for i = 1:L]; [[3i+3;3i+5;3] for i = 1:L]])
        end
        # if d > 1
        #     append!(basis[2], [[[3i+2;3i+6] for i = 1:L]; [[3i+3;3i+5] for i = 1:L]])
        # end
    end
    return basis
end

function get_pfbasis_loct(L, label)
    if label == 1
        basis = Vector{UInt16}[[], [1], [2]]
    else
        basis = Vector{UInt16}[[3i+1] for i = 1:L]
    end
    return basis
end

function get_pfbasis_loc1(L, label)
    if label == 1
        basis = Vector{UInt16}[[], [3]]
    else
        basis = Vector{UInt16}[[3i+1] for i = 1:L]
    end
    return basis
end

function get_pfbasis0(L, label)
    if label == 1
        basis = Vector{UInt16}[[[], [1], [2], [1;1], [1;2], [2;2]]; [sort([3i+1;smod(3i+1, 3L)+3]) for i = 1:L]]
    elseif label == 2
        basis = Vector{UInt16}[[[3i+1] for i = 1:L]; [[1;3i+1] for i = 1:L]; [[3i+1;2] for i = 1:L]; [[2;3i+1] for i = 1:L]]
    elseif label == 3
        basis = Vector{UInt16}[[[], [3], [3;3]]; [sort([3i+1;smod(3i+1, 3L)+3]) for i = 1:L]]
    else
        basis = Vector{UInt16}[[[3i+1] for i = 1:L]; [[3i+1;3] for i = 1:L]; [[3;3i+1] for i = 1:L]]
    end
    return basis
end
