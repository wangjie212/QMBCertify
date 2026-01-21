function bound_partfunc(supp::Vector{Vector{Int}}, coe::Vector{Float64}, L::Int, d::Int; QUIET=false)
    println("*********************************** QMBCertify ***********************************")
    println("QMBCertify is launching...")
    if QUIET == false
        println("Generating the monomial basis...")
    end
    time = @elapsed begin
    basis = Vector{Vector{Vector{UInt16}}}(undef, 2)
    tsupp = [UInt16[]]
    coe1 = Vector{Vector{Vector{Int8}}}(undef, 2)
    bi1 = Vector{Vector{Vector{Vector{UInt16}}}}(undef, 2)
    coe2 = Vector{Vector{Vector{Int8}}}(undef, 2)
    bi2 = Vector{Vector{Vector{Vector{UInt16}}}}(undef, 2)
    for i = 1:2
        basis[i] = get_basis(L, i-1, d, lattice=lattice, extra=extra, three_type=three_type)
        k = Int(length(basis[i])/L)
        # processing diagonal blocks
        coe1[i] = Vector{Vector{Int8}}(undef, k)
        bi1[i] = Vector{Vector{Vector{UInt16}}}(undef, k)
        for j = 1:k
            coe1[i][j] = Vector{Int8}(undef, Int(L/2))
            bi1[i][j] = Vector{Vector{UInt16}}(undef, Int(L/2))
            for r = 1:Int(L/2)
                bi = [basis[i][L*(j-1)+1]; basis[i][L*(j-1)+r+1]]
                bi1[i][j][r], coe1[i][j][r] = reduce!(bi, L=L, lattice=lattice, realify=true)
                if coe1[i][j][r] != 0
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
                bi = [basis[i][L*(j1-1)+1]; basis[i][L*(j2-1)+r]]
                bi2[i][j][r],coe2[i][j][r] = reduce!(bi, L=L, lattice=lattice, realify=true)
                if coe2[i][j][r] != 0
                    push!(tsupp, bi2[i][j][r])
                end
            end
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp = length(tsupp)
    end
    if QUIET == false
        println("Obtained the monomial basis in $time seconds.")
        println("Starting block-diagonalization...")
    end
    time = @elapsed begin
    model = Model(optimizer_with_attributes(Mosek.Optimizer, QUIET))
    cons = [AffExpr(0) for i=1:ltsupp]
    mb = 0
    gram = Vector{Vector{Symmetric{VariableRef}}}(undef, 2)
    for i = 1:2
        k = Int(length(basis[i])/L)
        mb = max(2*k, mb)
        gram[i] = Vector{Symmetric{VariableRef}}(undef, Int(L/2)+1)
        for l = 1:Int(L/2) + 1
            if i == 1 && l == 1
                gram[i][1] = @variable(model, [1:k+1, 1:k+1], PSD)
                @inbounds add_to_expression!(cons[1], gram[i][1][1,1])
                for j = 1:k
                    bi,coef = reduce!(basis[1][L*(j-1)+1], L=L, lattice=lattice)
                    if coef != 0
                        Locb = bfind(tsupp, ltsupp, bi)
                        @inbounds add_to_expression!(cons[Locb], 2*sqrt(L), gram[i][1][1,j+1])
                    end
                end
            elseif l == 1 || l == Int(L/2) + 1
                gram[i][l] = @variable(model, [1:k, 1:k], PSD)
            else
                gram[i][l] = @variable(model, [1:2*k, 1:2*k], PSD)
            end
            for j = 1:k
                if i == 1 && l == 1
                    pp = gram[i][l][j+1, j+1]
                elseif l == 1 || l == Int(L/2) + 1
                    pp = gram[i][l][j, j]
                else
                    pp = gram[i][l][j, j] + gram[i][l][j+k, j+k]
                end
                @inbounds add_to_expression!(cons[1], pp)
                if coe1[i][j][end] != 0
                    Locb = bfind(tsupp, ltsupp, bi1[i][j][end])
                    @inbounds add_to_expression!(cons[Locb], coe1[i][j][end]*(-1)^(l-1), pp)
                end
                for r = 1:Int(L/2)-1
                    if coe1[i][j][r] != 0
                        Locb = bfind(tsupp, ltsupp, bi1[i][j][r])
                        @inbounds add_to_expression!(cons[Locb], 2*coe1[i][j][r]*cos(2*pi*r*(l-1)/L), pp)
                    end
                end
            end
            for j1 = 1:k-1, j2 = j1+1:k
                if i == 1 && l == 1
                    pp1 = gram[i][l][j1+1, j2+1]
                elseif l == 1 || l == Int(L/2) + 1
                    pp1 = gram[i][l][j1, j2]
                else
                    pp1 = gram[i][l][j1, j2] + gram[i][l][j1+k, j2+k]
                    pp2 = gram[i][l][j1+k, j2] - gram[i][l][j2+k, j1]
                end
                j = Int((2*k-j1)*(j1-1)/2) + j2 - j1
                for r = 1:L
                    if coe2[i][j][r] != 0
                        Locb = bfind(tsupp, ltsupp, bi2[i][j][r])
                        @inbounds add_to_expression!(cons[Locb], 2*coe2[i][j][r]*cos(2*pi*(r-1)*(l-1)/L), pp1)
                        if l != 1 && l != Int(L/2) + 1
                            @inbounds add_to_expression!(cons[Locb], -2*coe2[i][j][r]*sin(2*pi*(r-1)*(l-1)/L), pp2)
                        end
                    end
                end                
            end
        end
    end
    end
    if QUIET == false
        println("Finished block-diagonalization in $time seconds.")
        println("SDP size: n = $mb, m = $ltsupp")
    end
    obj = @variable(model, lower)
    @objective(model, Max, obj)
    for i = 1:length(supp)
        Locb = bfind(tsupp, ltsupp, supp[i])
        if Locb === nothing
           @error "The monomial basis is not enough!"
           return nothing
        else
           cons[Locb] -= coe[i]
        end
    end
    cons[1] += lower
    @constraint(model, con, cons==zeros(ltsupp))
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
    return objv
end