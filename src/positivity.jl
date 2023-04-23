function posepsd6!(model, cons, tsupp, L; sites=[1;2;3;4;5;6])
    ltsupp = length(tsupp)
    Pauli = [[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    pos = @variable(model, [1:2^6, 1:2^6], PSD)
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3
        mon = mono([i,j,k,l,s,t], sites=sites)
        if !isz(mon)
            Locb = bfind(tsupp, ltsupp, reduce4(mon, L))
            tp = tensor([Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1], Pauli[s+1], Pauli[t+1]])
            @inbounds add_to_expression!(cons[Locb], sum(tp.*pos))
        end
    end
end

function posepsd8!(model, cons, tsupp, L)
    ltsupp = length(tsupp)
    Pauli = [[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    pos = @variable(model, [1:2^8, 1:2^8], PSD)
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3
        mon = mono([i,j,k,l,s,t,u,v])
        if !isz(mon)
            Locb = bfind(tsupp, ltsupp, reduce4(mon, L))
            tp = tensor([Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1], Pauli[s+1], Pauli[t+1], Pauli[u+1], Pauli[v+1]])
            @inbounds add_to_expression!(cons[Locb], sum(tp.*pos))
        end
    end
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
