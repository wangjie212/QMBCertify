function posepsd!(model, mvar, tsupp, L)
    nb = 6
    ltsupp = length(tsupp)
    Pauli = Vector{Matrix{ComplexF16}}(undef, 4)
    Pauli[1] = [1 0; 0 1]
    Pauli[2] = [0 1; 1 0]
    Pauli[3] = [0 -im; im 0]
    Pauli[4] = [1 0; 0 -1]
    pos = zeros(GenericAffExpr{Float64,VariableRef}, 2^nb, 2^nb)
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3
        mon = mono([i,j,k,l,s,t])
        if !iszero(mon)
            Locb = bfind(tsupp, ltsupp, reduce4(mon, L))
            tp = tensor([Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1], Pauli[s+1], Pauli[t+1]])
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

function mono(indices)
    mon = UInt16[]
    for i = 1:length(indices)
        if indices[i] != 0
            push!(mon, 3*(i-1)+indices[i])
        end
    end
    return mon
end
