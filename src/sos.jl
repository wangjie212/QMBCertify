function GSE(supp::Vector{Vector{UInt16}}, coe::Vector{Float64}, L::Int, d::Int; energy=[], QUIET=false, lattice="chain",
    posepsd=false, extra=0, three_type=[1;1], J2=0, correlation=false)
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
        # for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3
        #     sites = [[1;3;4;5;6;7], [1;2;4;5;6;7], [1;2;3;5;6;7]]
        #     for r = 1:3
        #         mon = mono([i,j,k,l,s,t], sites=sites[r])
        #         if !isz(mon)
        #             push!(tsupp, reduce4(mon, L))
        #         end
        #     end
        # end
        for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3
            mon = mono([i,j,k,l,s,t,u,v])
            if !isz(mon)
                push!(tsupp, reduce4(mon, L))
            end
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp = length(tsupp)
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:ltsupp]
    mb = 0
    for i = 0:3
        if i == 0
            k = Int((length(basis[1])-1)/L)
        else
            k = Int(length(basis[i+1])/L)
        end
        pos = Vector{Symmetric{VariableRef}}(undef, L)
        for j = 1:L
            if i == 0 && j == 1
                pos[j] = @variable(model, [1:2*(k+1), 1:2*(k+1)], PSD)
                if 2*(k+1) > mb
                    mb = 2*(k+1)
                end
                @inbounds add_to_expression!(cons[1], pos[1][1,1]+pos[1][k+2,k+2])
                for l = 1:k
                    bi,coef = reduce!(basis[1][L*(l-1)+2], L=L, lattice=lattice)
                    if coef != 0
                        Locb = bfind(tsupp, ltsupp, bi)
                        @inbounds add_to_expression!(cons[Locb], 2*sqrt(L), pos[1][1,l+1]+pos[1][k+2,l+k+2])
                    end
                end
            else
                pos[j] = @variable(model, [1:2*k, 1:2*k], PSD)
                if 2*k > mb
                    mb = 2*k
                end
            end
        end
        for j = 1:k
            coef = Vector{ComplexF64}(undef, Int(L/2))
            Locb = Vector{Int}(undef, Int(L/2))
            for r = 1:Int(L/2)
                if i == 0
                    bi = [basis[1][L*(j-1)+2]; basis[1][L*(j-1)+r+2]]
                else
                    bi = [basis[i+1][L*(j-1)+1]; basis[i+1][L*(j-1)+r+1]]
                end
                bi,coef[r] = reduce!(bi, L=L, lattice=lattice)
                if coef[r] != 0
                    Locb[r] = bfind(tsupp, ltsupp, bi)
                else
                    Locb[r] = 0
                end
            end
            for l = 1:L
                if i == 0 && l == 1
                    pp = pos[l][j+1, j+1] + pos[l][j+k+2, j+k+2]
                else
                    pp = pos[l][j, j]+pos[l][j+k, j+k]
                end
                @inbounds add_to_expression!(cons[1], pp)
                if coef[end] != 0
                    @inbounds add_to_expression!(cons[Locb[end]], (-1)^(l-1), pp)
                end
                for r = 1:Int(L/2)-1
                    if coef[r]^2 == 1
                        @inbounds add_to_expression!(cons[Locb[r]], 2*coef[r]*cos(2*pi*r*(l-1)/L), pp)
                    elseif coef[r] == im 
                        @inbounds add_to_expression!(cons[Locb[r]], -2*sin(2*pi*r*(l-1)/L), pp)
                    elseif coef[r] == -im 
                        @inbounds add_to_expression!(cons[Locb[r]], 2*sin(2*pi*r*(l-1)/L), pp)
                    end
                end
            end
        end
        for j1 = 1:k-1, j2 = j1+1:k
            coef = Vector{ComplexF64}(undef, L)
            Locb = Vector{Int}(undef, L)
            for r = 1:L
                if i == 0
                    bi = [basis[1][L*(j1-1)+2]; basis[1][L*(j2-1)+r+1]]
                else
                    bi = [basis[i+1][L*(j1-1)+1]; basis[i+1][L*(j2-1)+r]]
                end
                bi,coef[r] = reduce!(bi, L=L, lattice=lattice)
                if coef[r] != 0
                    Locb[r] = bfind(tsupp, ltsupp, bi)
                else
                    Locb[r] = 0
                end
            end
            for l = 1:L
                if i == 0 && l == 1
                    pp1 = pos[l][j1+1, j2+1] + pos[l][j1+k+2, j2+k+2]
                    pp2 = pos[l][j1+k+2, j2+1] - pos[l][j2+k+2, j1+1]
                else
                    pp1 = pos[l][j1, j2] + pos[l][j1+k, j2+k]
                    pp2 = pos[l][j1+k, j2] - pos[l][j2+k, j1]
                end
                for r = 1:L
                    if coef[r]^2 == 1
                        @inbounds add_to_expression!(cons[Locb[r]], 2*coef[r]*cos(2*pi*(r-1)*(l-1)/L), pp1)
                        @inbounds add_to_expression!(cons[Locb[r]], -2*coef[r]*sin(2*pi*(r-1)*(l-1)/L), pp2)
                    elseif coef[r] == im
                        @inbounds add_to_expression!(cons[Locb[r]], -2*sin(2*pi*(r-1)*(l-1)/L), pp1)
                        @inbounds add_to_expression!(cons[Locb[r]], -2*cos(2*pi*(r-1)*(l-1)/L), pp2)
                    elseif coef[r] == -im
                        @inbounds add_to_expression!(cons[Locb[r]], 2*sin(2*pi*(r-1)*(l-1)/L), pp1)
                        @inbounds add_to_expression!(cons[Locb[r]], 2*cos(2*pi*(r-1)*(l-1)/L), pp2)
                    end
                end
            end
        end
    end
    bc = zeros(ltsupp)
    for i = 1:length(supp)
        Locb = bfind(tsupp, ltsupp, supp[i])
        if Locb == 0
           @error "The monomial basis is not enough!"
           return nothing,nothing,nothing,nothing
        else
           bc[Locb] = coe[i]
        end
    end
    if posepsd == true && lattice == "chain"
        # posepsd6!(model, cons, tsupp, L)
        # posepsd6!(model, cons, tsupp, L, sites=[1;3;4;5;6;7])
        # posepsd6!(model, cons, tsupp, L, sites=[1;2;4;5;6;7])
        # posepsd6!(model, cons, tsupp, L, sites=[1;2;3;5;6;7])
        posepsd8!(model, cons, tsupp, L)
        if mb > 256
            mb = 256
        end
    end
    println("SDP size: n = $mb, m = $ltsupp")
    @variable(model, lower)
    cons[1] += lower
    @constraint(model, con[i=1:ltsupp], cons[i]==bc[i])
    @objective(model, Max, lower)
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
            if lattice == "chain"
                gsen += 3/4*J2*mvar[Locb]
            else
                gsen += 3/2*J2*mvar[Locb]
            end
        end
        @constraint(model, gsen>=energy[1])
        @constraint(model, gsen<=energy[2])
    end
    optimize!(model)
    status = termination_status(model)
    objv = objective_value(model)
    if status != MOI.OPTIMAL
       println("termination status: $status")
       status = primal_status(model)
       println("solution status: $status")
    end
    println("optimum = $objv")
    mvar = -dual.(con)
    if correlation == true
        if lattice == "chain"
            cor0 = zeros(Int(L/2))
            for i = 1:Int(L/2)
                word = UInt16[1; 3*i+1]
                Locb = bfind(tsupp, ltsupp, word)
                cor0[i] = value(mvar[Locb])
            end
            cor1 = zeros(Int(L/2-2))
            for i = 3:Int(L/2)
                word = UInt16[1; 4; 3*(i-1)+1; 3*i+1]
                word = reduce!(word, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, ltsupp, word)
                cor1[i-2] = value(mvar[Locb])
            end
            cor2 = zeros(Int(L/2-2))
            for i = 3:Int(L/2)
                word = UInt16[1; 4; 3*i; 3*i+3]
                word = reduce!(word, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, ltsupp, word)
                cor2[i-2] = value(mvar[Locb])
            end
        else
            cor0 = zeros(L, L)
            for i = 1:L, j = 1:L
                word = UInt16[1; 3*(slabel(i, j, L=L)-1)+1]
                word = reduce!(word, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, ltsupp, word)
                cor0[i,j] = value(mvar[Locb])
            end
            cor1 = cor2 = nothing
        end
    else
        cor0 = cor1 = cor2 = nothing
    end
    return objv,cor0,cor1,cor2
end
