using Arblib
using ITensors, ITensorMPS

function build_LHS_corr_rat(data_corr, N;
                            J2::Real,
                            e_lb_rat::Rational{BigInt},
                            e_ub_rat::Rational{BigInt},
                            obj_supp::Vector{Vector{Int}},
                            obj_coe_rat::Vector{Rational{BigInt}},
                            corr_bound_rat::Rational{BigInt},
                            upper::Bool = true,
                            tol_mult::Real = 1e-12)

    tsupp = data_corr.tsupp
    T     = length(tsupp)
    LHS   = zeros(Rational{BigInt}, T)

    @inline function reduce_word(w::Vector{UInt16})
        QMBCertify.reduce!(w; L = N, realify = true)
    end

    const_idx = bfind(tsupp, UInt16[])
    const_idx === nothing && error("Constant word [] not found in tsupp")

    @assert length(obj_supp) == length(obj_coe_rat)
    for (w_int, c_rat) in zip(obj_supp, obj_coe_rat)
        w_u16 = UInt16.(w_int)
        w_nf, c_nf = reduce_word(w_u16)
        abs(imag(c_nf)) > 1e-10 &&
            error("Unexpected imaginary coefficient in objective reduction")
        loc = bfind(tsupp, w_nf)
        loc === nothing && continue
        r_nf = Rational{BigInt}(rationalize(real(c_nf); tol = tol_mult))
        LHS[loc] += c_rat * r_nf
    end

    λ_rat = upper ? (+corr_bound_rat) : (-corr_bound_rat)
    LHS[const_idx] += λ_rat

    H_supp = [[1; 4], [1; 7]]
    J2_rat = Rational{BigInt}(rationalize(J2; tol = tol_mult))
    H_coe_rat = Rational{BigInt}[3//4, (3//4) * J2_rat]

    E_vec = zeros(Rational{BigInt}, T)
    for (w_int, cH_rat) in zip(H_supp, H_coe_rat)
        w_u16 = UInt16.(w_int)
        w_nf, c_nf = reduce_word(w_u16)
        abs(imag(c_nf)) > 1e-10 &&
            error("Unexpected imaginary coefficient in H reduction")
        loc = bfind(tsupp, w_nf)
        loc === nothing && continue
        r_nf = Rational{BigInt}(rationalize(real(c_nf); tol = tol_mult))
        E_vec[loc] += cH_rat * r_nf
    end

    mu1_f, mu2_f = data_corr.multiplier
    mu1_rat = Rational{BigInt}(rationalize(Float64(mu1_f); tol = tol_mult))
    mu2_rat = Rational{BigInt}(rationalize(Float64(mu2_f); tol = tol_mult))

    @inbounds for t in 1:T
        LHS[t] -= mu1_rat * E_vec[t]
    end
    LHS[const_idx] += mu1_rat * e_lb_rat

    @inbounds for t in 1:T
        LHS[t] += mu2_rat * E_vec[t]
    end
    LHS[const_idx] -= mu2_rat * e_ub_rat

    return LHS
end

function certify_qmb_corr(
    N::Int,
    d_E::Int,
    d_corr::Int;
    J2::Real = 0.0,
    dist = 1,
    extra_E::Int = 0,
    extra_corr::Int = 0,
    tol_gram::Real   = 1e-15,
    tol_dft::Real    = 1e-12,
    tol_E::Real      = 1e-12,
    digits_dmrg::Int = 12,
    QUIET::Bool      = false,
    check::Bool      = true,
    eig_prec::Int = 256,
)
    J2f = Float64(J2)

    H_supp = [[1; 4], [1; 7]]
    coe_E  = Float64[0.75, 0.75 * J2f]

    opt_E, data_energy = GSB(
        H_supp, coe_E, N, d_E;
        QUIET = QUIET,
        extra = extra_E,
        lso   = 0,
        pso   = 0,
        rdm   = 0,
        lol   = N,
        Gram  = true,
    )

    energy_cert = certify_qmb(
        data_energy, N, coe_E[1], opt_E;
        tol_gram = tol_gram,
        tol_dft  = tol_dft,
        snn      = true,
        J2       = J2f,
        check    = check,
        eig_prec = eig_prec,
    )

    e_lb_rat = energy_cert.newbound_rat

    e_ub_rat_int = dmrg_heisenberg_rat(N, 1.0; J2 = J2, digits = digits_dmrg)
    e_ub_rat = Rat(e_ub_rat_int)

    e_lb = Float64(e_lb_rat)
    e_ub = Float64(e_ub_rat)

    println("Rigorous energy bounds per site (floated): [$e_lb, $e_ub]")

    energy_bounds = Float64[e_lb, e_ub]

    dists = dist isa Integer ? Int[dist] : collect(dist)

    function corr_branch_poly(d::Int, dir::Symbol)
        println("Running correlation branch dir=$(dir), dist=$(d)")

        coe_corr = (dir == :upper) ? Float64[-0.25] : Float64[+0.25]
        supp_corr = Vector{Vector{Int}}([[1; 3*(1 + d) - 2]])

        opt_corr, data_corr = GSB(
            supp_corr, coe_corr, N, d_corr;
            H_supp = H_supp,
            H_coe  = coe_E,
            energy = energy_bounds,
            J2     = J2f,
            QUIET  = QUIET,
            extra  = extra_corr,
            lso    = 0,
            pso    = 0,
            rdm    = 0,
            lol    = N,
            Gram   = true,
        )

        C_sdp = (dir == :upper) ? -Float64(opt_corr) : Float64(opt_corr)
        println("  SDP $(dir) bound (numeric) for dist=$(d): ", C_sdp)

        corr_bound_rat = Rat(rationalize(C_sdp; tol = tol_gram))

        obj_supp = supp_corr
        obj_coe_rat = if dir == :upper
            Rat[-1//4]
        else
            Rat[+1//4]
        end

        LHS_vec = build_LHS_corr_rat(
            data_corr, N;
            J2             = J2f,
            e_lb_rat       = e_lb_rat,
            e_ub_rat       = e_ub_rat,
            obj_supp       = obj_supp,
            obj_coe_rat    = obj_coe_rat,
            corr_bound_rat = corr_bound_rat,
            upper          = (dir == :upper),
            tol_mult       = tol_gram,
        )

        t_r = @elapsed begin
            G1_blocks, G2_blocks =
                round_grams(data_corr, N; tol_gram = tol_gram)
        end
        println("  Time to round small Gram blocks (corr, $(dir), dist=$(d)): ", t_r, " s")

        G1_blocks_before = deepcopy(G1_blocks)
        G2_blocks_before = deepcopy(G2_blocks)

        t_rhs0 = @elapsed begin
            RHS_before, N_before, gram_maps =
                build_rhs_and_N_rat(data_corr, N, G1_blocks, G2_blocks; tol_dft = tol_dft)

            resid_before = RHS_before .- LHS_vec

        end
        println("  Time to build RHS & N (pre-projection, $(dir), dist=$(d)): ", t_rhs0, " s")

        t_proj = @elapsed begin
            G1_blocks, G2_blocks = project_blocks_frob!(
                data_corr, N,
                G1_blocks, G2_blocks,
                LHS_vec, RHS_before, N_before, gram_maps;
                tol_dft = tol_dft,
            )
        end
        println("  Time to project small blocks (corr, $(dir), dist=$(d)): ", t_proj, " s")

        if check
            t_rhs = @elapsed begin
                RHS_chk = apply_A_to_grams_rat(
                    data_corr, N,
                    G1_blocks, G2_blocks,
                    gram_maps;
                    tol_dft = tol_dft,
                )
                resid = RHS_chk .- LHS_vec

            end
            println("  Time to rebuild RHS for check (corr, $(dir), dist=$(d)): ", t_rhs, " s")
        end

        eigminsg1 = Vector{Float64}()
        eigminsg2 = Vector{Float64}()
        eigminsg1_rat = Vector{Rat}()
        eigminsg2_rat = Vector{Rat}()

        @inbounds for block in G1_blocks
            Gproj = ComplexF64.(block)
            λ_proj_rat, λ_proj = rigorous_min_eig_bound(Hermitian(Gproj); prec = eig_prec)
            push!(eigminsg1_rat, λ_proj_rat)
            push!(eigminsg1, Float64(λ_proj))
        end

        @inbounds for block in G2_blocks
            Gproj = ComplexF64.(block)
            λ_proj_rat, λ_proj = rigorous_min_eig_bound(Hermitian(Gproj); prec = eig_prec)
            push!(eigminsg2_rat, λ_proj_rat)
            push!(eigminsg2, Float64(λ_proj))
        end

        lambda_min_g1 = minimum(eigminsg1)
        lambda_min_g2 = minimum(eigminsg2)
        lambda_min_g1_rat = minimum(eigminsg1_rat)
        lambda_min_g2_rat = minimum(eigminsg2_rat)

        println("  CORR $(dir), dist=$(d): Minimum eigenvalue over G1 small blocks = ",
                lambda_min_g1)

        println("  CORR $(dir), dist=$(d): Minimum eigenvalue over G2 small blocks = ",
                lambda_min_g2)

        # crude but rigorous shift: use global min eigenvalue per family
        shift_rat = lambda_min_g1_rat * length(data_corr.basis[1]) +
                    lambda_min_g2_rat * length(data_corr.basis[2])

        C_bound_corr_rat = if dir == :upper
            corr_bound_rat - shift_rat
        else
            corr_bound_rat + shift_rat
        end

        println("  CORR $(dir), dist=$(d): shift = ", Float64(shift_rat))
        println("  CORR $(dir), dist=$(d): old SDP bound = ", C_sdp)
        println("  CORR $(dir), dist=$(d): new rigorous bound = ",
                Float64(C_bound_corr_rat))

        return (
            C_sdp        = C_sdp,
            C_bound_rat  = C_bound_corr_rat,
            shift_corr   = shift_rat,
            eigmins_g1   = eigminsg1,
            eigmins_g2   = eigminsg2,
            eigmins_g1_rat = eigminsg1_rat,
            eigmins_g2_rat = eigminsg2_rat,
        )
    end

    ups  = [corr_branch_poly(d, :upper) for d in dists]
    lows = [corr_branch_poly(d, :lower) for d in dists]

    return (
        e_lb       = e_lb_rat,
        e_ub       = e_ub_rat,
        dist           = dists,
        C_num_upper    = [u.C_sdp       for u in ups],
        C_rat_upper  = [u.C_bound_rat for u in ups],
        shift_upper    = [u.shift_corr  for u in ups],
        C_num_lower    = [l.C_sdp       for l in lows],
        C_rat_lower  = [l.C_bound_rat for l in lows],
        shift_lower    = [l.shift_corr  for l in lows],
    )
end
