function GSB(supp::Vector{Vector{Int}}, coe::Vector{Float64}, L::Int, d::Int; H_supp=supp, H_coe=coe, SU2_symmetry=false, lso=true, lol=L, pso=3, energy=[], QUIET=false, lattice="chain",
    rdm=false, extra=0, three_type=[1;1], J2=0, correlation=false, Gram=false, writetofile=false, mosek_setting=mosek_para())
    println("*********************************** QMBCertify ***********************************")
    println("QMBCertify is launching...")
    if QUIET == false
        println("Generating the monomial basis...")
    end
    time = @elapsed begin
    basis = Vector{Vector{Vector{UInt16}}}(undef, 2)
    if pso > 0
        gb = Vector{Vector{Vector{UInt16}}}(undef, 2)
    end
    tsupp = [UInt16[]]
    if lattice == "chain" 
        coe1 = Vector{Vector{Vector{Int8}}}(undef, 2)
        bi1 = Vector{Vector{Vector{Vector{UInt16}}}}(undef, 2)
        coe2 = Vector{Vector{Vector{Int8}}}(undef, 2)
        bi2 = Vector{Vector{Vector{Vector{UInt16}}}}(undef, 2)
        if pso > 0
            coe3 = Vector{Vector{Vector{Vector{Float16}}}}(undef, 2)
            bi3 = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, 2)
            coe4 = Vector{Vector{Vector{Vector{Float16}}}}(undef, 2)
            bi4 = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, 2)
        end
    else
        veig1 = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, 2)
        ceig1 = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, 2)
        nveig1 = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, 2)
        nceig1 = Vector{Vector{Vector{Union{Vector{Float64},Vector{ComplexF64}}}}}(undef, 2)
        if pso > 0
            veig2 = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, 2)
            ceig2 = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, 2)
            nveig2 = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, 2)
            nceig2 = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, 2)
        end
    end
    for i = 1:2
        basis[i] = get_basis(L, i-1, d, lattice=lattice, extra=extra, three_type=three_type)
        if lattice == "chain"
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
        else
            k = Int(length(basis[i])/L^2)
            # implementing the first block-diagonalization
            sveig = Vector{Vector{Vector{Vector{UInt16}}}}(undef, Int(k*(k+1)/2)*L)
            sceig = Vector{Vector{Union{Vector{Float64},Vector{ComplexF64}}}}(undef, Int(k*(k+1)/2)*L)
            for j = 1:k
                # processing diagonal blocks
                ind = Int(((j-1)*(2k-j+2)*L)/2) + 1
                ctemp = Vector{Vector{Int8}}(undef, Int(L/2)+1)
                vtemp = Vector{Vector{Vector{UInt16}}}(undef, Int(L/2)+1)
                for r = 1:Int(L/2)+1
                    bi = [basis[i][L^2*(j-1)+1]; basis[i][L^2*(j-1)+r]]
                    v1,c1 = reduce!(bi, L=L, lattice=lattice, realify=true)
                    vtemp[r], ctemp[r] = [v1], [c1]
                end
                sveig[ind], sceig[ind] = eigen_circmat(vtemp, ctemp, L, symmetry=true, real_matrix=true)
                # processing non-diagonal blocks
                for q = 1:(k-j+1)*L - 1
                    ctemp = Vector{Vector{Int8}}(undef, L)
                    vtemp = Vector{Vector{Vector{UInt16}}}(undef, L)
                    for r = 1:L
                        bi = [basis[i][L^2*(j-1)+1]; basis[i][L^2*(j-1)+q*L+r]]
                        v1,c1 = reduce!(bi, L=L, lattice=lattice, realify=true)
                        vtemp[r], ctemp[r] = [v1], [c1]
                    end
                    sveig[ind+q], sceig[ind+q] = eigen_circmat(vtemp, ctemp, L, real_matrix=true)
                end
            end
            # implementing the second block-diagonalization
            veig1[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, Int(k*(k+1)/2)*(Int(L/2)-1))
            ceig1[i] = Vector{Vector{Vector{ComplexF64}}}(undef, Int(k*(k+1)/2)*(Int(L/2)-1))
            nveig1[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, k*(k+1))
            nceig1[i] = Vector{Vector{Union{Vector{Float64},Vector{ComplexF64}}}}(undef, k*(k+1))
            for j = 1:k, l = 0:k-j, p = 1:Int(L/2)-1
                ind = Int(((j-1)*(2k-j+2)*L)/2) + l*L
                ind1 = Int(((j-1)*(2k-j+2)*(Int(L/2)-1))/2) + l*(Int(L/2)-1)
                vtemp = [sveig[ind+q][p+1] for q=1:L]
                ctemp = [sceig[ind+q][p+1] for q=1:L]
                veig1[i][ind1+p], ceig1[i][ind1+p] = eigen_circmat(vtemp, ctemp, L)
                append!(tsupp, veig1[i][ind1+p]...)
            end
            for j = 1:k, l = 0:k-j, (s,p) in enumerate([1, Int(L/2)+1])
                ind = Int(((j-1)*(2k-j+2)*L)/2) + l*L
                ind1 = (j-1)*(2k-j+2) + l*2
                vtemp = [sveig[ind+q][p] for q=1:L]
                ctemp = [sceig[ind+q][p] for q=1:L]
                nveig1[i][ind1+s], nceig1[i][ind1+s] = eigen_circmat(vtemp, ctemp, L, real_matrix=true)
                append!(tsupp, nveig1[i][ind1+s]...)
            end
        end
        if pso > 0
            h = [ceil(Int, item[2]/3) - 1 for item in H_supp]
            gb[i] = basis[i][length.(basis[i]) .<= pso]
            if lattice == "chain"
                k = Int(length(gb[i])/L)
                coe3[i] = Vector{Vector{Vector{Float16}}}(undef, k)
                bi3[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, k)
                for j = 1:k
                    coe3[i][j] = Vector{Vector{Float16}}(undef, Int(L/2)+1)
                    bi3[i][j] = Vector{Vector{Vector{UInt16}}}(undef, Int(L/2)+1)
                    for r = 1:Int(L/2)+1
                        bi3[i][j][r], coe3[i][j][r] = PSDstate_entry(gb[i][L*(j-1)+1], gb[i][L*(j-1)+r], h, H_coe, L, lattice=lattice)
                    end
                end
                coe4[i] = Vector{Vector{Vector{Float16}}}(undef, Int(k*(k-1)/2))
                bi4[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, Int(k*(k-1)/2))
                for j1 = 1:k-1, j2 = j1+1:k
                    j = Int((2*k-j1)*(j1-1)/2) + j2 - j1
                    coe4[i][j] = Vector{Vector{Float16}}(undef, L)
                    bi4[i][j] = Vector{Vector{Vector{UInt16}}}(undef, L)
                    for r = 1:L
                        bi4[i][j][r], coe4[i][j][r] = PSDstate_entry(gb[i][L*(j1-1)+1], gb[i][L*(j2-1)+r], h, H_coe, L, lattice=lattice)
                        append!(tsupp, bi4[i][j][r])
                    end        
                end
            else
                k = Int(length(gb[i])/L^2)
                sveig = Vector{Vector{Vector{Vector{UInt16}}}}(undef, Int(k*(k+1)/2)*L)
                sceig = Vector{Vector{Union{Vector{Float64},Vector{ComplexF64}}}}(undef, Int(k*(k+1)/2)*L)
                for j = 1:k
                    ind = Int(((j-1)*(2k-j+2)*L)/2) + 1
                    ctemp = Vector{Vector{Float16}}(undef, Int(L/2)+1)
                    vtemp = Vector{Vector{Vector{UInt16}}}(undef, Int(L/2)+1)
                    for r = 1:Int(L/2)+1
                        vtemp[r], ctemp[r] = PSDstate_entry(gb[i][L^2*(j-1)+1], gb[i][L^2*(j-1)+r], h, H_coe, L, lattice=lattice)
                    end
                    sveig[ind], sceig[ind] = eigen_circmat(vtemp, ctemp, L, symmetry=true, real_matrix=true)
                    for q = 1:(k-j+1)*L - 1
                        ctemp = Vector{Vector{Float16}}(undef, L)
                        vtemp = Vector{Vector{Vector{UInt16}}}(undef, L)
                        for r = 1:L
                            vtemp[r], ctemp[r] = PSDstate_entry(gb[i][L^2*(j-1)+1], gb[i][L^2*(j-1)+q*L+r], h, H_coe, L, lattice=lattice)
                        end
                        sveig[ind+q], sceig[ind+q] = eigen_circmat(vtemp, ctemp, L, real_matrix=true)
                    end
                end
                veig2[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, Int(k*(k+1)/2)*(Int(L/2)-1))
                ceig2[i] = Vector{Vector{Vector{ComplexF64}}}(undef, Int(k*(k+1)/2)*(Int(L/2)-1))
                nveig2[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, k*(k+1))
                nceig2[i] = Vector{Vector{Union{Vector{Float64},Vector{ComplexF64}}}}(undef, k*(k+1))
                for j = 1:k, l = 0:k-j, p = 1:Int(L/2)-1
                    ind = Int(((j-1)*(2k-j+2)*L)/2) + l*L
                    ind1 = Int(((j-1)*(2k-j+2)*(Int(L/2)-1))/2) + l*(Int(L/2)-1)
                    vtemp = [sveig[ind+q][p+1] for q=1:L]
                    ctemp = [sceig[ind+q][p+1] for q=1:L]
                    veig2[i][ind1+p], ceig2[i][ind1+p] = eigen_circmat(vtemp, ctemp, L)
                    append!(tsupp, veig2[i][ind1+p]...)
                end
                for j = 1:k, l = 0:k-j, (s,p) in enumerate([1, Int(L/2)+1])
                    ind = Int(((j-1)*(2k-j+2)*L)/2) + l*L
                    ind1 = (j-1)*(2k-j+2) + l*2
                    vtemp = [sveig[ind+q][p] for q=1:L]
                    ctemp = [sceig[ind+q][p] for q=1:L]
                    nveig2[i][ind1+s], nceig2[i][ind1+s] = eigen_circmat(vtemp, ctemp, L, real_matrix=true)
                    append!(tsupp, nveig2[i][ind1+s]...)
                end
            end
        end
    end
    if rdm != false
        if rdm == 8
            for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3
                ind = [i,j,k,l,s,t,u,v]
                if all(x->iseven(sum(ind .== x)), 1:3)
                    inx = ind .!= 0
                    push!(tsupp, reduce4(UInt16.(3*(Vector(1:8)[inx] .- 1) + ind[inx]), L, lattice=lattice))
                end
            end
        elseif rdm == 9
            for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3, w = 0:3
                ind = [i,j,k,l,s,t,u,v,w]
                if all(x->iseven(sum(ind .== x)), 1:3)
                    inx = ind .!= 0
                    push!(tsupp, reduce4(UInt16.(3*(Vector(1:9)[inx] .- 1) + ind[inx]), L, lattice=lattice))
                end
            end
        elseif rdm == 10
            for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3, w = 0:3, z = 0:3
                ind = [i,j,k,l,s,t,u,v,w,z]
                if all(x->iseven(sum(ind .== x)), 1:3)
                    inx = ind .!= 0
                    push!(tsupp, reduce4(UInt16.(3*(Vector(1:10)[inx] .- 1) + ind[inx]), L, lattice=lattice))
                end
            end
        else
            println("Adding rdm > 10 is not supported!")
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
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:ltsupp]
    mb = 0
    gram = Vector{Vector{Symmetric{VariableRef}}}(undef, 2)
    for i = 1:2
        if lattice == "chain"
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
        else
            k = Int(length(basis[i])/L^2)
            mb = max(2*k, mb)
            pos = Matrix{Symmetric{VariableRef}}(undef, Int(L/2)-1, L)
            npos = Matrix{Symmetric{VariableRef}}(undef, 2, Int(L/2)+1)
            for l = 1:Int(L/2)-1, u = 1:L
                pos[l, u] = @variable(model, [1:2*k, 1:2*k], PSD)
            end
            for l = 1:2, u = 1:Int(L/2)+1
                if i == 1 && l == 1 && u == 1
                    npos[1, 1] = @variable(model, [1:k+1, 1:k+1], PSD)
                    @inbounds add_to_expression!(cons[1], npos[1,1][1,1])
                    for j = 1:k
                        bi,coef = reduce!(basis[1][L^2*(j-1)+1], L=L, lattice=lattice)
                        if coef != 0
                            Locb = bfind(tsupp, ltsupp, bi)
                            @inbounds add_to_expression!(cons[Locb], 2*L, npos[1,1][1,j+1])
                        end
                    end
                elseif u == 1 || u == Int(L/2)+1
                    npos[l, u] = @variable(model, [1:k, 1:k], PSD)
                else
                    npos[l, u] = @variable(model, [1:2k, 1:2k], PSD)
                end
            end
            for j1 = 1:k, j2 = j1:k, l = 1:Int(L/2)-1
                ind = Int(((j1-1)*(2k-j1+2)*(Int(L/2)-1))/2) + (j2-j1)*(Int(L/2)-1)
                for u = 1:L
                    pp1 = pos[l, u][j1, j2] + pos[l, u][j1+k, j2+k]
                    pp2 = pos[l, u][j1+k, j2] - pos[l, u][j2+k, j1]
                    for (s,c) in enumerate(ceig1[i][ind+l][u])
                        Locb = bfind(tsupp, ltsupp, veig1[i][ind+l][u][s])
                        if j1 == j2
                            @inbounds add_to_expression!(cons[Locb], real(c), pp1)
                        else
                            @inbounds add_to_expression!(cons[Locb], 2*real(c), pp1)
                            @inbounds add_to_expression!(cons[Locb], -2*imag(c), pp2)
                        end
                    end
                end
            end
            for j1 = 1:k, j2 = j1:k, l = 1:2
                ind = (j1-1)*(2k-j1+2) + (j2-j1)*2
                for u = 1:Int(L/2)+1
                    if i == 1 && l == 1 && u == 1
                        pp1 = npos[l, u][j1+1, j2+1]
                    elseif u == 1 || u == Int(L/2)+1
                        pp1 = npos[l, u][j1, j2]
                    else
                        pp1 = npos[l, u][j1, j2] + npos[l, u][j1+k, j2+k]
                        pp2 = npos[l, u][j1+k, j2] - npos[l, u][j2+k, j1]
                    end
                    for (s,c) in enumerate(nceig1[i][ind+l][u])
                        Locb = bfind(tsupp, ltsupp, nveig1[i][ind+l][u][s])
                        if j1 == j2
                            @inbounds add_to_expression!(cons[Locb], real(c), pp1)
                        else
                            @inbounds add_to_expression!(cons[Locb], 2*real(c), pp1)
                            if u != 1 && u != Int(L/2)+1
                                @inbounds add_to_expression!(cons[Locb], -2*imag(c), pp2)
                            end
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
    if SU2_symmetry == true
        add_SU2_equality!(model, tsupp, ltsupp, cons, L=L, lattice=lattice)
    end
    if lso == true
        if QUIET == false
            println("Adding linear state optimality constraints...")
        end
        time = @elapsed begin
            h = [ceil(Int, item[2]/3) - 1 for item in H_supp]
            if lattice == "chain"
                mons = generate_mons(L, lol, rdm-1)
                mons = filter_mons(mons, tsupp, h, H_coe, L, lattice="chain")
                for mon in mons
                    fr = @variable(model)
                    for (u, v) in enumerate(h), i = 1:L, j = 1:3
                        word = UInt16[3*(i-1)+j; 3*mod(i-1+v, L)+j; mon]
                        word,coef = reduce!(word, L=L, lattice="chain")
                        if imag(coef) != 0
                            Locb = bfind(tsupp, ltsupp, word)
                            add_to_expression!(cons[Locb], H_coe[u]*imag(coef), fr)
                        end
                    end
                end
            else
                mons = generate_mons(L^2, lol, 0)
                mons = filter_mons(mons, tsupp, h, H_coe, L, lattice="square")
                for mon in mons
                    fr = @variable(model)
                    for (u, v) in enumerate(h), i = 1:L, w = 1:L, j = 1:3
                        word = UInt16[3*(slabel(i, w, L=L)-1)+j; 3*(slabel(i+v, w, L=L)-1)+j; mon]
                        word,coef = reduce!(word, L=L, lattice="square")
                        if imag(coef) != 0
                            Locb = bfind(tsupp, ltsupp, word)
                            add_to_expression!(cons[Locb], H_coe[u]*imag(coef), fr)
                        end
                        word = UInt16[3*(slabel(i, w, L=L)-1)+j; 3*(slabel(i, w+v, L=L)-1)+j; mon]
                        word,coef = reduce!(word, L=L, lattice="square")
                        if imag(coef) != 0
                            Locb = bfind(tsupp, ltsupp, word)
                            add_to_expression!(cons[Locb], H_coe[u]*imag(coef), fr)
                        end
                    end
                end
            end
        end
        if QUIET == false
            println("Added linear state optimality constraints in $time seconds.")
        end
    end
    if pso > 0
        if QUIET == false
            println("Adding PSD state optimality constraints...")
        end
        time = @elapsed begin
            if lattice == "chain"
                for i = 1:2
                    k = Int(length(gb[i])/L)
                    pos = Vector{Symmetric{VariableRef}}(undef, Int(L/2)+1)
                    for l = 1:Int(L/2)+1
                        if l == 1 || l == Int(L/2)+1
                            pos[l] = @variable(model, [1:k, 1:k], PSD)
                        else
                            pos[l] = @variable(model, [1:2*k, 1:2*k], PSD)
                        end
                        for j = 1:k
                            if l == 1 || l == Int(L/2)+1
                                pp = pos[l][j, j]
                            else
                                pp = pos[l][j, j] + pos[l][j+k, j+k]
                            end
                            for s = 1:length(coe3[i][j][1])
                                Locb = bfind(tsupp, ltsupp, bi3[i][j][1][s])
                                @inbounds add_to_expression!(cons[Locb], coe3[i][j][1][s], pp)
                            end
                            for r = 1:Int(L/2)-1, s = 1:length(coe3[i][j][r+1])
                                Locb = bfind(tsupp, ltsupp, bi3[i][j][r+1][s])
                                @inbounds add_to_expression!(cons[Locb], 2*coe3[i][j][r+1][s]*cos(2*pi*r*(l-1)/L), pp)
                            end
                            for s = 1:length(coe3[i][j][end])
                                Locb = bfind(tsupp, ltsupp, bi3[i][j][end][s])
                                @inbounds add_to_expression!(cons[Locb], coe3[i][j][end][s]*(-1)^(l-1), pp)
                            end
                        end
                        for j1 = 1:k-1, j2 = j1+1:k
                            if l == 1 || l == Int(L/2)+1
                                pp1 = pos[l][j1, j2]
                            else
                                pp1 = pos[l][j1, j2] + pos[l][j1+k, j2+k]
                                pp2 = pos[l][j1+k, j2] - pos[l][j2+k, j1]
                            end
                            j = Int((2*k-j1)*(j1-1)/2) + j2 - j1
                            for r = 1:L, s = 1:length(coe4[i][j][r])
                                Locb = bfind(tsupp, ltsupp, bi4[i][j][r][s])
                                if coe4[i][j][r][s] != 0
                                    @inbounds add_to_expression!(cons[Locb], 2*coe4[i][j][r][s]*cos(2*pi*(r-1)*(l-1)/L), pp1)
                                    if l != 1 && l != Int(L/2)+1
                                        @inbounds add_to_expression!(cons[Locb], -2*coe4[i][j][r][s]*sin(2*pi*(r-1)*(l-1)/L), pp2)
                                    end
                                end
                            end
                        end
                    end
                end
            else
                for i = 1:2
                    k = Int(length(gb[i])/L^2)
                    pos = Matrix{Symmetric{VariableRef}}(undef, Int(L/2)-1, L)
                    npos = Matrix{Symmetric{VariableRef}}(undef, 2, Int(L/2)+1)
                    for l = 1:Int(L/2)-1, u = 1:L
                        pos[l, u] = @variable(model, [1:2*k, 1:2*k], PSD)
                    end
                    for l = 1:2, u = 1:Int(L/2)+1
                        if u == 1 || u == Int(L/2)+1
                            npos[l, u] = @variable(model, [1:k, 1:k], PSD)
                        else
                            npos[l, u] = @variable(model, [1:2k, 1:2k], PSD)
                        end
                    end
                    for j1 = 1:k, j2 = j1:k, l = 1:Int(L/2)-1
                        ind = Int(((j1-1)*(2k-j1+2)*(Int(L/2)-1))/2) + (j2-j1)*(Int(L/2)-1)
                        for u = 1:L
                            pp1 = pos[l, u][j1, j2] + pos[l, u][j1+k, j2+k]
                            pp2 = pos[l, u][j1+k, j2] - pos[l, u][j2+k, j1]
                            for (s,c) in enumerate(ceig2[i][ind+l][u])
                                Locb = bfind(tsupp, ltsupp, veig2[i][ind+l][u][s])
                                if j1 == j2
                                    @inbounds add_to_expression!(cons[Locb], real(c), pp1)
                                else
                                    @inbounds add_to_expression!(cons[Locb], 2*real(c), pp1)
                                    @inbounds add_to_expression!(cons[Locb], -2*imag(c), pp2)
                                end
                            end
                        end
                    end
                    for j1 = 1:k, j2 = j1:k, l = 1:2
                        ind = (j1-1)*(2k-j1+2) + (j2-j1)*2
                        for u = 1:Int(L/2)+1
                            if u == 1 || u == Int(L/2)+1
                                pp1 = npos[l, u][j1, j2]
                            else
                                pp1 = npos[l, u][j1, j2] + npos[l, u][j1+k, j2+k]
                                pp2 = npos[l, u][j1+k, j2] - npos[l, u][j2+k, j1]
                            end
                            for (s,c) in enumerate(nceig2[i][ind+l][u])
                                Locb = bfind(tsupp, ltsupp, nveig2[i][ind+l][u][s])
                                if j1 == j2
                                    @inbounds add_to_expression!(cons[Locb], real(c), pp1)
                                else
                                    @inbounds add_to_expression!(cons[Locb], 2*real(c), pp1)
                                    if u != 1 && u != Int(L/2)+1
                                        @inbounds add_to_expression!(cons[Locb], -2*imag(c), pp2)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        if QUIET == false
            println("Added PSD state optimality constraints in $time seconds.")
        end
    end
    if rdm != false
        if QUIET == false
            println("Adding the rdm constraint...")
        end
        time = @elapsed begin
        if rdm == 8
            posepsd8!(model, cons, tsupp, L, lattice=lattice)
        elseif rdm == 9
            posepsd9!(model, cons, tsupp, L, lattice=lattice)
        elseif rdm == 10
            posepsd10!(model, cons, tsupp, L, lattice=lattice)
        end
        end
        if QUIET == false
            println("Added the rdm constraint in $time seconds.")
        end
    end
    obj = @variable(model, lower)
    if energy != []
        mul = @variable(model, [1:2], lower_bound=0)
        Locb = bfind(tsupp, ltsupp, [1;4])
        if lattice == "chain"
            @inbounds add_to_expression!(cons[Locb], 3/4, mul[1]-mul[2])
        else
            @inbounds add_to_expression!(cons[Locb], 3/2, mul[1]-mul[2])
        end
        if J2 != 0
            Locb = bfind(tsupp, ltsupp, [1;7])
            if lattice == "chain"
                @inbounds add_to_expression!(cons[Locb], 3/4*J2, mul[1]-mul[2])
            else
                @inbounds add_to_expression!(cons[Locb], 3/2*J2, mul[1]-mul[2])
            end
        end
        obj += energy[1]*mul[1]
        obj -= energy[2]*mul[2]
    end
    @objective(model, Max, obj)
    for i = 1:length(supp)
        Locb = bfind(tsupp, ltsupp, supp[i])
        if Locb === nothing
           @error "The monomial basis is not enough!"
           return nothing,nothing,nothing,nothing
        else
           cons[Locb] -= coe[i]
        end
    end
    cons[1] += lower
    @constraint(model, con, cons==zeros(ltsupp))
    if writetofile != false
        write_to_file(dualize(model), writetofile)
    end
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
    cor0 = cor1 = cor2 = GramMat = multiplier = moment = nothing
    mvar = -dual(con) # extract moments
    if correlation == true
        if lattice == "chain"
            cor0 = zeros(Int(L/2))
            for i = 1:Int(L/2)
                word = UInt16[1; 3*i+1]
                Locb = bfind(tsupp, ltsupp, word)
                cor0[i] = mvar[Locb]
            end
            cor1 = zeros(Int(L/2-2))
            for i = 3:Int(L/2)
                word = UInt16[1; 4; 3*(i-1)+1; 3*i+1]
                word = reduce!(word, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, ltsupp, word)
                cor1[i-2] = mvar[Locb]
            end
            cor2 = zeros(Int(L/2-2))
            for i = 3:Int(L/2)
                word = UInt16[1; 4; 3*i; 3*i+3]
                word = reduce!(word, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, ltsupp, word)
                cor2[i-2] = mvar[Locb]
            end
        else
            cor0 = zeros(L, L)
            for i = 1:L, j = 1:L
                word = UInt16[1; 3*(slabel(i, j, L=L)-1)+1]
                word = reduce!(word, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, ltsupp, word)
                cor0[i,j] = mvar[Locb]
            end
        end
    end
    if Gram == true
        GramMat = Vector{Vector{Union{Matrix{Float64}, Matrix{ComplexF64}}}}(undef, 2)
        for i = 1:2
            k = Int(length(basis[i])/L)
            GramMat[i] = Vector{Union{Matrix{Float64}, Matrix{ComplexF64}}}(undef, Int(L/2)+1)
            for l = 1:Int(L/2)+1
                if l == 1 || l == Int(L/2) + 1
                    GramMat[i][l] = value.(gram[i][l])
                else
                    A = value.(gram[i][l])
                    GramMat[i][l] = A[1:k,1:k] + A[k+1:2k,k+1:2k] + (A[k+1:2k,1:k] - A[1:k,k+1:2k])*im
                end
            end
        end
        if energy != []
            multiplier = value.(mul)
        end
    end
    # moment = Vector{Vector{Union{Matrix{Float64},Matrix{ComplexF64}}}}(undef, 2)
    # for i = 1:2
    #     moment[i] = Vector{Union{Matrix{Float64},Matrix{ComplexF64}}}(undef, Int(L/2) + 1)
    #     k = Int(length(basis[i])/L)
    #     for l = 1:Int(L/2) + 1
    #         if i == 1 && l == 1
    #             moment[1][1] = zeros(k+1, k+1)
    #             moment[1][1][1,1] += mvar[1]
    #             for j = 1:k
    #                 bi,coef = reduce!(basis[1][L*(j-1)+1], L=L, lattice=lattice)
    #                 if coef != 0
    #                     Locb = bfind(tsupp, ltsupp, bi)
    #                     moment[1][1][1,j+1] += sqrt(L)*mvar[Locb]
    #                 end
    #             end
    #         elseif l == 1 || l == Int(L/2) + 1
    #             moment[i][l] = zeros(k, k)
    #         else
    #             moment[i][l] = zeros(ComplexF64, k, k)
    #         end
    #         for j = 1:k
    #             if i == 1 && l == 1
    #                 moment[1][1][j+1, j+1] += mvar[1]
    #             else
    #                 moment[i][l][j, j] += mvar[1]
    #             end
    #             if coe1[i][j][end] != 0
    #                 Locb = bfind(tsupp, ltsupp, bi1[i][j][end])
    #                 if i == 1 && l == 1
    #                     moment[1][1][j+1, j+1] += coe1[1][j][end]*mvar[Locb]
    #                 else
    #                     moment[i][l][j, j] += coe1[i][j][end]*(-1)^(l-1)*mvar[Locb]
    #                 end
    #             end
    #             for r = 1:Int(L/2)-1
    #                 if coe1[i][j][r] != 0
    #                     Locb = bfind(tsupp, ltsupp, bi1[i][j][r])
    #                     if i == 1 && l == 1
    #                         moment[1][1][j+1, j+1] += 2*coe1[1][j][r]*mvar[Locb]
    #                     else
    #                         moment[i][l][j, j] += 2*coe1[i][j][r]*cos(2*pi*r*(l-1)/L)*mvar[Locb]
    #                     end
    #                 end
    #             end
    #         end
    #         for j1 = 1:k-1, j2 = j1+1:k
    #             j = Int((2*k-j1)*(j1-1)/2) + j2 - j1
    #             for r = 1:L
    #                 if coe2[i][j][r] != 0
    #                     Locb = bfind(tsupp, ltsupp, bi2[i][j][r])
    #                     if i == 1 && l == 1
    #                         moment[1][1][j1+1, j2+1] += coe2[1][j][r]*mvar[Locb]
    #                     elseif l == 1 || l == Int(L/2) + 1
    #                         moment[i][l][j1, j2] += coe2[i][j][r]*cos(2*pi*(r-1)*(l-1)/L)*mvar[Locb]
    #                     else
    #                         moment[i][l][j1, j2] += coe2[i][j][r]*(cos(2*pi*(r-1)*(l-1)/L)+sin(2*pi*(r-1)*(l-1)/L)*im)*mvar[Locb]
    #                     end
    #                 end
    #             end
    #         end
    #         moment[i][l] += moment[i][l]'
    #         moment[i][l] -= diagm(diag(moment[i][l]))/2
    #     end
    # end
    data = qmb_data(cor0,cor1,cor2,basis,tsupp,GramMat,multiplier,moment)
    return objv,data
end

function resort(Locb, coef)
    nLocb = copy(Locb)
    sort!(nLocb)
    unique!(nLocb)
    ncoef = zeros(typeof(coef[1]), length(nLocb))
    for i = 1:length(Locb)
        ind = bfind(nLocb, length(nLocb), Locb[i])
        ncoef[ind] += coef[i]
    end
    ind = ncoef .!= 0
    return nLocb[ind],ncoef[ind]
end

function update!(coef, bis, bi, L; λ=1, lattice="chain")
    bi,co = reduce!(bi, L=L, lattice=lattice, realify=true)
    if co != 0
        push!(coef, λ*co)
        push!(bis, bi)
    end
end

function PSDstate_entry(ba1, ba2, h, c, L; lattice="chain")
    coef = Float16[]
    bis = Vector{UInt16}[]
    if lattice == "chain"
        for (i, k) in enumerate(h), u = 1:L, v = 1:3
            s = count_reduce(u, v, k, ba1, L)
            t = count_reduce(u, v, k, ba2, L)
            if s + t > 0
                bi = UInt16[ba1; 3*(u-1)+v; 3*mod(u-1+k, L)+v; ba2]
                update!(coef, bis, bi, L, λ=c[i]*(s+t), lattice="chain")
            end
        end
    else
        for (i, k) in enumerate(h), u = 1:L, w = 1:L, v = 1:3
            s = count_reduce(u, w, v, k, ba1, L)
            t = count_reduce(u, w, v, k, ba2, L)
            if s[1] + t[1] > 0
                bi = UInt16[ba1; 3*(slabel(u, w, L=L)-1)+v; 3*(slabel(u+k, w, L=L)-1)+v; ba2]
                update!(coef, bis, bi, L, λ=c[i]*(s[1]+t[1]), lattice="square")
            end
            if s[2] + t[2] > 0
                bi = UInt16[ba1; 3*(slabel(u, w, L=L)-1)+v; 3*(slabel(u, w+k, L=L)-1)+v; ba2]
                update!(coef, bis, bi, L, λ=c[i]*(s[2]+t[2]), lattice="square")
            end
        end
    end
    if !isempty(coef)
        bis, coef = resort(bis, coef)
    end
    return bis, coef
end

function generate_mons(N, b1, b2)
    mons = Vector{UInt16}[]
    for i = 2:N-1, j = i+1:N
        push!(mons, [1; 3*(i-1)+2; 3*(j-1)+3])
    end
    for i = 2:b1-1, j = i+1:b1, k = 2:b1-1, l = k+1:b1
        if length(unique([i;j;k;l])) == 4
            for t = 1:3
                push!(mons, sort([1; 3*(i-1)+2; 3*(j-1)+3; 3*(k-1)+t; 3*(l-1)+t]))
            end
        end
    end
    for i = 2:b2-1, j = i+1:b2, k = 2:b2-3, l = k+1:b2-2, u = l+1:b2-1, v = u+1:b2
        if length(unique([i;j;k;l;u;v])) == 6
            for t = 1:3
                push!(mons, sort([1; 3*(i-1)+2; 3*(j-1)+3; 3*(k-1)+t; 3*(l-1)+t; 3*(u-1)+t; 3*(v-1)+t]))
            end
        end
    end
    for i = 2:b2-1, j = i+1:b2, k = 2:b2-1, l = k+1:b2, u = 2:b2-1, v = u+1:b2
        if length(unique([i;j;k;l;u;v])) == 6
            for s = 1:2, t = s+1:3
                push!(mons, sort([1; 3*(i-1)+2; 3*(j-1)+3; 3*(k-1)+s; 3*(l-1)+s; 3*(u-1)+t; 3*(v-1)+t]))
            end
        end
    end
    unique!(mons)
    return mons
end

function filter_mons(mons, tsupp, h, c, L; lattice="chain")
    ltsupp = length(tsupp)
    lc = zeros(Float64, length(mons))
    rd = ceil.(Int, rand(ltsupp)*10^8)
    if lattice == "chain"
        for (k, mon) in enumerate(mons)
            flag = 1
            for (u, v) in enumerate(h), i = 1:L, j = 1:3
                word = UInt16[3*(i-1)+j; 3*mod(i-1+v, L)+j; mon]
                word,coef = reduce!(word, L=L, lattice="chain")
                if imag(coef) != 0
                    Locb = bfind(tsupp, ltsupp, word)
                    if Locb === nothing
                        flag = 0
                        break
                    end
                    lc[k] += c[u]*imag(coef)*rd[Locb]
                end
            end
            if flag == 0
                lc[k] = 0
            end
        end
    else
        for (k, mon) in enumerate(mons)
            flag = 1
            for (u, v) in enumerate(h), i = 1:L, w = 1:L, j = 1:3
                word = UInt16[3*(slabel(i, w, L=L)-1)+j; 3*(slabel(i+v, w, L=L)-1)+j; mon]
                word,coef = reduce!(word, L=L, lattice="square")
                if imag(coef) != 0
                    Locb = bfind(tsupp, ltsupp, word)
                    if Locb === nothing
                        flag = 0
                        break
                    end
                    lc[k] += c[u]*imag(coef)*rd[Locb]
                end
                word = UInt16[3*(slabel(i, w, L=L)-1)+j; 3*(slabel(i, w+v, L=L)-1)+j; mon]
                word,coef = reduce!(word, L=L, lattice="square")
                if imag(coef) != 0
                    Locb = bfind(tsupp, ltsupp, word)
                    if Locb === nothing
                        flag = 0
                        break
                    end
                    lc[k] += c[u]*imag(coef)*rd[Locb]
                end
            end
            if flag == 0
                lc[k] = 0
            end
        end
    end
    ind = unique(i -> abs(lc[i]), eachindex(lc))
    return mons[ind[abs.(lc[ind]) .> 0]]
end

function count_reduce(i, s, t, ba, L)
    loc = ceil.(UInt16, ba/3)
    u = bfind(loc, length(loc), i)
    v = bfind(loc, length(loc), smod(i+t, L))
    return (u !== nothing && mod(ba[u], 3) != mod(s, 3)) + (v !== nothing && mod(ba[v], 3) != mod(s, 3)) == 1
end

function count_reduce(i, j, s, t, ba, L)
    loc = ceil.(UInt16, ba/3)
    u = bfind(loc, length(loc), slabel(i, j, L=L))
    v = bfind(loc, length(loc), slabel(i+t, j, L=L))
    w = bfind(loc, length(loc), slabel(i, j+t, L=L))
    a = u !== nothing && mod(ba[u], 3) != mod(s, 3)
    b = v !== nothing && mod(ba[v], 3) != mod(s, 3)
    c = w !== nothing && mod(ba[w], 3) != mod(s, 3)
    return [a + b == 1, a + c == 1]
end
