function posepsd!(model, mvar, tsupp, L)
    nb = 8
    ltsupp = length(tsupp)
    Pauli = Vector{Matrix{ComplexF16}}(undef, 4)
    Pauli[1] = [1 0; 0 1]
    Pauli[2] = [0 1; 1 0]
    Pauli[3] = [0 -im; im 0]
    Pauli[4] = [1 0; 0 -1]
    rfbcons = zeros(GenericAffExpr{Float64,VariableRef}, 2^nb, 2^nb)
    ifbcons = zeros(GenericAffExpr{Float64,VariableRef}, 2^nb, 2^nb)
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3
        mon = mono([i,j,k,l,s,t,u,v])
        if !iszero(mon)
            Locb = bfind(tsupp, ltsupp, reduce4(mon, L))
            tp = tensor([Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1], Pauli[s+1], Pauli[t+1], Pauli[u+1], Pauli[v+1]])
            for p = 1:2^nb, q = p:2^nb
                rtpcoe,itpcoe = reim(tp[p,q])
                if rtpcoe != 0
                    rfbcons[p,q] += rtpcoe*mvar[Locb]
                elseif itpcoe != 0
                    ifbcons[p,q] += itpcoe*mvar[Locb]
                end
            end
        end
    end
    pos = [AffExpr(0) for i=1:2*2^nb, j=1:2*2^nb]
    pos[1:2^nb,1:2^nb] = rfbcons
    pos[2^nb+1:2*2^nb,2^nb+1:2*2^nb] = rfbcons
    pos[1:2^nb,2^nb+1:2*2^nb] = ifbcons' - ifbcons
    @constraint(model, Symmetric(pos) in PSDCone())
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
