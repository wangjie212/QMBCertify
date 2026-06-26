using Arblib

function certify_qmb(data, Nsites, J, numopt;
                           tol_gram::Real = 1e-12,
                           tol_dft::Real  = 1e-20,
                           snn::Bool=false, J2::Real=1.0,
                           check::Bool=false, eig_prec::Int=256)

    GC.gc()

    G1_num_blocks = data.GramMat[1]
    G2_num_blocks = data.GramMat[2]

    res = round_project_qmb(data, Nsites, J, numopt;
                                  tol_gram=tol_gram,
                                  tol_dft=tol_dft,
                                  snn=snn, J2=J2, check=check)

    G1_blocks_proj = res.G1_blocks_proj
    G2_blocks_proj = res.G2_blocks_proj
    LHS_vec        = res.LHS_rat

    max_rel = 0.0

    @inbounds for (i, block) in enumerate(G1_num_blocks)
        Gnum = ComplexF64.(block)  
        F, rel = frob(Gnum, G1_blocks_proj[i])
        max_rel = max(max_rel, rel)
    end

    @inbounds for (i, block) in enumerate(G2_num_blocks)
        Gnum = ComplexF64.(block)
        F, rel = frob(Gnum, G2_blocks_proj[i])
        max_rel = max(max_rel, rel)
    end

    eigminsg1 = Vector{Float64}()
    eigminsg2 = Vector{Float64}()
    eigminsg1_rat = Vector{Rat}()
    eigminsg2_rat = Vector{Rat}()

    @inbounds for i in 1:length(G1_num_blocks)
        Gproj = ComplexF64.(G1_blocks_proj[i])
        lambda_proj_rat, lambda_proj = rigorous_min_eig_bound(Hermitian(Gproj); prec=eig_prec)
        push!(eigminsg1_rat, lambda_proj_rat)
        push!(eigminsg1, Float64(lambda_proj))
    end

    lambda_min_g1 = minimum(eigminsg1)
    lambda_min_g1_rat = minimum(eigminsg1_rat)

    @inbounds for i in 1:length(G2_num_blocks)
        Gproj = ComplexF64.(G2_blocks_proj[i])
        lambda_proj_rat, lambda_proj = rigorous_min_eig_bound(Hermitian(Gproj); prec=eig_prec)
        push!(eigminsg2_rat, lambda_proj_rat)
        push!(eigminsg2, Float64(lambda_proj))
    end

    lambda_min_g2 = minimum(eigminsg2)
    lambda_min_g2_rat = minimum(eigminsg2_rat)

    shift = lambda_min_g1*length(data.basis[1]) + lambda_min_g2*length(data.basis[2])
    shift_rat = lambda_min_g1_rat*length(data.basis[1]) + lambda_min_g2_rat*length(data.basis[2])
    const_idx = bfind(data.tsupp, UInt16[])
    const_idx === nothing && error("Constant word [] not found in tsupp")
    oldbound_rat = -LHS_vec[const_idx]
    bound_rat = oldbound_rat + shift_rat
    bound = Float64(numopt) + shift

    println("Shift: ", shift)

    return (oldbound = numopt,
            newbound = bound,
            oldbound_rat = oldbound_rat,
            newbound_rat = bound_rat,
            shift    = shift,
            shift_rat = shift_rat,
            mineigs  = vcat(eigminsg1, eigminsg2),
            mineigs_rat = vcat(eigminsg1_rat, eigminsg2_rat),
            dims     = [length(data.basis[1]), length(data.basis[2])],)
end
