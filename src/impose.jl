function posepsd6!(model, mvar, tsupp, L; lattice="chain", sites=[1;2;3;4;5;6])
    ltsupp = length(tsupp)
    Pauli = [[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    pos = zeros(GenericAffExpr{Float64,VariableRef}, 2^6, 2^6)
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3
        mon = mono([i,j,k,l,s,t], sites=sites)
        if !iszero(mon)
            if lattice=="chain"
                Locb = bfind(tsupp, ltsupp, reduce4(mon, L))
            else
                Locb = bfind(tsupp, ltsupp, reduce5(reduce4(mon, L, lattice="square"), L))
            end
            tp = tensor([Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1], Pauli[s+1], Pauli[t+1]])
            pos += mvar[Locb]*tp
        end
    end
    @constraint(model, pos in PSDCone())
end

function posepsd8!(model, mvar, tsupp, L)
    ltsupp = length(tsupp)
    Pauli = [[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    pos = zeros(GenericAffExpr{Float64,VariableRef}, 2^8, 2^8)
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3
        mon = mono([i,j,k,l,s,t,u,v])
        if !iszero(mon)
            Locb = bfind(tsupp, ltsupp, reduce4(mon, L))
            tp = tensor([Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1], Pauli[s+1], Pauli[t+1], Pauli[u+1], Pauli[v+1]])
            pos += mvar[Locb]*tp
        end
    end
    @constraint(model, pos in PSDCone())
end

function tensor(matrices)
    O = matrices[1]
    for i = 2:length(matrices)
        O = kron(O, matrices[i])
    end
    return O
end

function mono(indices; sites=[1;2;3;4;5;6;7;8])
    mon = UInt16[]
    for i = 1:length(indices)
        if indices[i] != 0
            push!(mon, 3*(sites[i]-1)+indices[i])
        end
    end
    return mon
end
