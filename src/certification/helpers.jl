const Rat = Rational{BigInt}

struct GramNMaps

    nmap_re_1::Vector{Vector{Vector{Dict{Int,Rat}}}}
    nmap_im_1::Vector{Vector{Vector{Dict{Int,Rat}}}}

    nmap_re_2::Vector{Vector{Vector{Dict{Int,Rat}}}}
    nmap_im_2::Vector{Vector{Vector{Dict{Int,Rat}}}}
end

function var_from_array(var, arr)
    result = DynamicPolynomials.Monomial(var, zeros(Int, length(var)))
    if length(arr) == 0
        return result
    end
    for i in arr
        result *= var[i]
    end
    return result
end

function rationalize_poly(p::DynamicPolynomials.Polynomial; tol::Real = 1e-20)
    # Already rational?
    if eltype(p.a) <: Rational{BigInt} || eltype(p.a) <: Complex{Rational{BigInt}}
        return p
    end

    # Real coefficients
    if eltype(p.a) <: Real
        coeffs = Rational{BigInt}.(rationalize.(p.a; tol=tol))
        return DynamicPolynomials.Polynomial(coeffs, p.x)
    end

    # Complex coefficients (common in QMB code)
    if eltype(p.a) <: Complex
        coeffs = Vector{Complex{Rational{BigInt}}}(undef, length(p.a))
        @inbounds for k in eachindex(p.a)
            c = p.a[k]
            re = Rational{BigInt}(rationalize(real(c); tol=tol))
            im = Rational{BigInt}(rationalize(imag(c); tol=tol))
            coeffs[k] = Complex{Rational{BigInt}}(re, im)
        end
        return DynamicPolynomials.Polynomial(coeffs, p.x)
    end

    return p
end


function rigorous_min_eig_bound(m::AbstractMatrix; prec::Int = 128)
    mc = AcbMatrix(m; prec=64)

    ev_approx, R_approx = Arblib.approx_eig_qr(mc; prec=prec)

    eps = Arblib.eig_global_enclosure(mc, ev_approx, R_approx; prec=prec)

    eig_ball = minimum(Arb.(real.(ev_approx); prec=prec)) - Arblib.Arb(eps; prec=prec)
    eig_min_rat = Rational{BigInt}(Arblib.lbound(eig_ball))
    return eig_min_rat, BigFloat(eig_min_rat)
end

function rigorous_min_eig(m::AbstractMatrix; prec::Int = 128)
    return rigorous_min_eig_bound(m; prec=prec)[2]
end

function make_vars_chain(N::Int)
    @ncpolyvar z[1:N] x[1:N] y[1:N]
    vars = Vector{typeof(z[1])}(undef, 3N)
    @inbounds for i in 1:N
        vars[3i-2] = z[i]; vars[3i-1] = x[i]; vars[3i] = y[i]
    end
    return vars
end

function heisenberg_chain_poly(J, vars, N::Int; snn::Bool=false, J2::Real=1.0)
    H = zero(vars[1]*vars[1])
    @inbounds for i in 1:N
        j  = (i == N) ? 1 : i + 1    
        Zi, Xi, Yi = vars[3i-2], vars[3i-1], vars[3i]
        Zj, Xj, Yj = vars[3j-2], vars[3j-1], vars[3j]

        H += (J/(3N)) * (Xi*Xj + Yi*Yj + Zi*Zj)

        if snn
            k  = (i >= N-1) ? (i - (N-2)) : i + 2 
            Zk, Xk, Yk = vars[3k-2], vars[3k-1], vars[3k]
            H += (J*J2/(3N)) * (Xi*Xk + Yi*Yk + Zi*Zk)
        end
    end
    return H
end

function array_from_var(m::DynamicPolynomials.Monomial, vars)
    idxs = Int[]
    @inbounds for (v, e) in zip(m.vars, m.z)
        pos = findfirst(==(v), vars) :: Int
        for _ in 1:e
            push!(idxs, pos)
        end
    end
    return idxs
end

monomial_to_u16(m::DynamicPolynomials.Monomial, vars) = UInt16.(array_from_var(m, vars))
u16_to_monomial(a::AbstractVector{<:Integer}, vars) = var_from_array(vars, Int.(a))

function poly_nf_QMBC(p::DynamicPolynomials.Polynomial, vars; L::Int=0, realify=false)
    acc = Dict{DynamicPolynomials.Monomial, ComplexF64}()
    @inbounds for (c, m) in zip(p.a, p.x)
        a = monomial_to_u16(m, vars)
        a_nf, c_nf = QMBCertify.reduce!(a; L=L, realify=realify)
        m_nf = u16_to_monomial(a_nf, vars)
        acc[m_nf] = get(acc, m_nf, 0.0+0im) + complex(c) * complex(c_nf)
    end
    q = nothing
    @inbounds for (m, c) in acc
        t = DynamicPolynomials.Polynomial([c], [m]); q = q === nothing ? t : (q + t)
    end
    return q === nothing ? zero(p) : q
end

@inline function sqrtN_rat(N::Int; tol_dft::Real = 1e-20)
    s_f = sqrt(float(N))
    return Rational{BigInt}(rationalize(s_f; tol = tol_dft))
end

function dft_products_rat(
    P::AbstractMatrix{ComplexF64},
    Imax::Int,
    Nsites::Int,
    tol_dft::Real,
)
    x_rat = Array{Rat}(undef, Imax, Nsites, Nsites)
    y_rat = Array{Rat}(undef, Imax, Nsites, Nsites)

    @inbounds for i in 1:Imax, u in 1:Nsites, v in 1:Nsites
        Pprod = P[i,u] * conj(P[i,v])
        x_f = real(Pprod)
        y_f = imag(Pprod)
        x_rat[i,u,v] = Rat(rationalize(x_f; tol = tol_dft))
        y_rat[i,u,v] = Rat(rationalize(y_f; tol = tol_dft))
    end

    return x_rat, y_rat
end

function build_rhs_and_N_rat(
    data,
    Nsites::Int,
    G1_blocks::Vector{<:AbstractMatrix{Complex{Rat}}},
    G2_blocks::Vector{<:AbstractMatrix{Complex{Rat}}};
    tol_dft::Real = 1e-20,
)
    tsupp = data.tsupp
    T     = length(tsupp)

    RHS  = zeros(Rat, T)
    Nmat = zeros(Rat, T, T)

    s    = Int.(length.(data.basis) ./ Nsites)
    s1, s2 = s
    Imax = Int(Nsites ÷ 2) + 1
    two  = 2//1
    sqrtN = sqrtN_rat(Nsites; tol_dft = tol_dft)

    P = dft_rows(Imax, Nsites) 

    x_rat, y_rat = dft_products_rat(P, Imax, Nsites, tol_dft)

    loc_cache = Dict{Tuple{Vararg{UInt16}}, Int}()

    @inline function loc_cached(w_nf::Vector{UInt16})
        key = tk(w_nf)
        get!(loc_cache, key) do
            loc = bfind(tsupp, w_nf)
            (loc === nothing ? 0 : loc)
        end
    end

    pair_cache_1 = Dict{Tuple{Int,Int}, Tuple{Int, Rat}}()
    pair_cache_2 = Dict{Tuple{Int,Int}, Tuple{Int, Rat}}()

    @inline function pair_loc_coeff(l::Int, idx_i::Int, idx_k::Int)
        cache   = (l == 1 ? pair_cache_1 : pair_cache_2)
        basis_l = data.basis[l]
        key     = (idx_i, idx_k)

        return get!(cache, key) do
            wi = basis_l[idx_i]
            wk = basis_l[idx_k]

            w_nf, c_nf = QMBCertify.reduce!([wi; wk]; L = Nsites, realify = true)
            c_nf == 0 && return (0, 0//1)

            loc = bfind(tsupp, w_nf)
            (loc === nothing || loc == 0) && return (0, 0//1)

            c_rat = Rat(rationalize(real(c_nf); tol = tol_dft))
            return (loc, c_rat)
        end
    end

    idxbuf = Int[]

    @inline function update_N_from_nmaps!(
        N::Matrix{Rat},
        nmap_re::Dict{Int,Rat},
        nmap_im::Dict{Int,Rat},
    )
        isempty(nmap_re) && isempty(nmap_im) && return

        empty!(idxbuf)

        @inbounds for t in keys(nmap_re)
            push!(idxbuf, t)
        end

        @inbounds for t in keys(nmap_im)
            found = false
            for u in idxbuf
                if u == t
                    found = true
                    break
                end
            end
            found || push!(idxbuf, t)
        end

        Lloc = length(idxbuf)
        @inbounds for a in 1:Lloc
            ta = idxbuf[a]
            va_re = get(nmap_re, ta, 0//1)
            va_im = get(nmap_im, ta, 0//1)
            N[ta, ta] += va_re*va_re + va_im*va_im
            for b in (a+1):Lloc
                tb = idxbuf[b]
                vb_re = get(nmap_re, tb, 0//1)
                vb_im = get(nmap_im, tb, 0//1)
                contrib = va_re*vb_re + va_im*vb_im
                N[ta, tb] += contrib
                N[tb, ta] += contrib
            end
        end
    end

    nmap_re_1 = [ [ [ Dict{Int,Rat}() for _ in 1:s1 ] for _ in 1:s1 ]
                  for _ in 1:Imax ]
    nmap_im_1 = [ [ [ Dict{Int,Rat}() for _ in 1:s1 ] for _ in 1:s1 ]
                  for _ in 1:Imax ]
    nmap_re_2 = [ [ [ Dict{Int,Rat}() for _ in 1:s2 ] for _ in 1:s2 ]
                  for _ in 1:Imax ]
    nmap_im_2 = [ [ [ Dict{Int,Rat}() for _ in 1:s2 ] for _ in 1:s2 ]
                  for _ in 1:Imax ]

    G11 = G1_blocks[1]
    g00 = G11[1,1]
    g00_re = getfield(g00, :re)

    empty_idx = bfind(tsupp, UInt16[])
    @assert empty_idx !== nothing && empty_idx != 0 "empty word not in tsupp"
    t_empty = empty_idx

    RHS[t_empty] += g00_re

    nmap_const_re = Dict{Int,Rat}(t_empty => 1//1)
    nmap_const_im = Dict{Int,Rat}()
    update_N_from_nmaps!(Nmat, nmap_const_re, nmap_const_im)

    @inbounds for k in 1:s1
        wi = data.basis[1][Nsites*(k-1) + 1]
        w_nf, c_nf = QMBCertify.reduce!(wi; L = Nsites)
        c_nf == 0 && continue

        loc = loc_cached(w_nf)
        (loc == 0) && continue

        g = G11[1, k+1]
        g_re = getfield(g, :re)

        coeff = two * sqrtN
        RHS[loc] += coeff * g_re

        nmap_re = Dict{Int,Rat}(loc => coeff)
        nmap_im = Dict{Int,Rat}()
        update_N_from_nmaps!(Nmat, nmap_re, nmap_im)
    end

    @inbounds for l in 1:2
        basis_l  = data.basis[l]
        G_blocks = (l == 1 ? G1_blocks : G2_blocks)
        s_l      = s[l]

        nmap_re_l = (l == 1 ? nmap_re_1 : nmap_re_2)
        nmap_im_l = (l == 1 ? nmap_im_1 : nmap_im_2)

        for i in 1:Imax
            Gi = G_blocks[i]
            maps_re_i = nmap_re_l[i]
            maps_im_i = nmap_im_l[i]

            for j in 1:s_l, k in j:s_l

                g = if l == 1 && i == 1
                    Gi[j+1, k+1]
                else
                    Gi[j,   k]
                end
                g_re = getfield(g, :re)
                g_im = getfield(g, :im)
                iszero(g_re) && iszero(g_im) && continue

                nmap_re = maps_re_i[j][k]
                nmap_im = maps_im_i[j][k]
                empty!(nmap_re)
                empty!(nmap_im)

                α = (j == k) ? (1//1) : two

                for u in 1:Nsites, v in 1:Nsites

                    idx_i = Nsites*(j-1) + u
                    idx_k = Nsites*(k-1) + v

                    loc, c_rat = pair_loc_coeff(l, idx_i, idx_k)
                    c_rat == 0 && continue

                    x = x_rat[i, u, v]
                    y = y_rat[i, u, v]

                    n_re_inc = α * c_rat * x
                    n_im_inc = (l == 1 && i == 1) ? (0//1) : (-α * c_rat * y)

                    nmap_re[loc] = get(nmap_re, loc, 0//1) + n_re_inc
                    nmap_im[loc] = get(nmap_im, loc, 0//1) + n_im_inc

                    RHS[loc] += n_re_inc * g_re + n_im_inc * g_im
                end

                update_N_from_nmaps!(Nmat, nmap_re, nmap_im)
            end
        end
    end

    gram_maps = GramNMaps(nmap_re_1, nmap_im_1, nmap_re_2, nmap_im_2)
    return RHS, Nmat, gram_maps
end

tk(w::Vector{UInt16}) = Tuple(w)

function coeff_vector_LHS_tsupp_rat(LHS_poly_rat::DynamicPolynomials.Polynomial, vars, tsupp; N::Int)
    v = zeros(Rational{BigInt}, length(tsupp))
    loc_cache = Dict{Tuple{Vararg{UInt16}}, Int}()
    red_cache = Dict{Tuple{Vararg{UInt16}}, Tuple{Vector{UInt16},ComplexF64}}()
    @inline function reduce_cached(a)
        key = tk(a)
        (wnf,c) = get!(red_cache, key) do
            QMBCertify.reduce!(a; L=N, realify=true)
        end
        return wnf, c
    end
    @inline function loc_cached(wnf)
        key = tk(wnf)
        get!(loc_cache, key) do
            loc = bfind(tsupp, wnf)
            (loc === nothing ? 0 : loc)
        end
    end
    @inbounds for (c, m) in zip(LHS_poly_rat.a, LHS_poly_rat.x)
        w_u16 = monomial_to_u16(m, vars)
        w_nf, c_nf = reduce_cached(w_u16)
        loc = loc_cached(w_nf)
        (loc == 0) && continue
        v[loc] += c * Rational{BigInt}(rationalize(real(c_nf)))
    end
    return v
end

@inline function dft_rows(m::Int, N::Int)
    P = Array{ComplexF64}(undef, m, N)
    @inbounds for i in 1:m, u in 1:N
        P[i,u] = (1/sqrt(N)) * cis(-2π * (i-1)*(u-1)/N)
    end
    return P
end


function project_blocks_frob!(
    data,
    Nsites::Int,
    G1_blocks::Vector{Matrix{Complex{Rat}}},
    G2_blocks::Vector{Matrix{Complex{Rat}}},
    LHS_vec::Vector{Rat},
    RHS::Vector{Rat},
    Nmat::Matrix{Rat},
    gram_maps::GramNMaps;
    tol_dft::Real = 1e-20,
)
    tsupp = data.tsupp
    s    = Int.(length.(data.basis) ./ Nsites)
    s1, s2 = s
    Imax = Int(Nsites ÷ 2) + 1

    r = LHS_vec .- RHS

    t_Ninv = @elapsed begin
        Delta = Nmat \ r   
    end
    println("Time to solve N Delta = r in rationals: $(t_Ninv) seconds")

    empty_idx = bfind(tsupp, UInt16[])
    @assert empty_idx !== nothing && empty_idx != 0
    t_empty = empty_idx

    Delta0 = Delta[t_empty]
    if Delta0 != 0//1
        G11 = G1_blocks[1]
        old = G11[1,1]
        re_old = getfield(old, :re)
        im_old = getfield(old, :im)
        G11[1,1] = Complex{Rat}(re_old + Delta0, im_old)
    end

    two  = 2//1
    sqrtN = sqrtN_rat(Nsites; tol_dft = tol_dft)
    G11 = G1_blocks[1]

    @inbounds for k in 1:s1
        wi = data.basis[1][Nsites*(k-1) + 1]
        w_nf, c_nf = QMBCertify.reduce!(wi; L = Nsites)
        c_nf == 0 && continue

        loc = bfind(tsupp, w_nf)
        (loc === nothing || loc == 0) && continue

        Deltat = Delta[loc]
        Deltat == 0//1 && continue

        coeff = two * sqrtN
        Deltag_re = Deltat * coeff

        old = G11[1, k+1]
        re_old = getfield(old, :re)
        im_old = getfield(old, :im)
        new_val = Complex{Rat}(re_old + Deltag_re, im_old)
        G11[1, k+1] = new_val
        G11[k+1, 1] = conj(new_val)
    end

    nmap_re_1 = gram_maps.nmap_re_1
    nmap_im_1 = gram_maps.nmap_im_1
    nmap_re_2 = gram_maps.nmap_re_2
    nmap_im_2 = gram_maps.nmap_im_2

    @inbounds for l in 1:2
        G_blocks = (l == 1 ? G1_blocks : G2_blocks)
        s_l      = (l == 1 ? s1 : s2)
        nmap_re_l = (l == 1 ? nmap_re_1 : nmap_re_2)
        nmap_im_l = (l == 1 ? nmap_im_1 : nmap_im_2)

        for i in 1:Imax
            Gi = G_blocks[i]
            maps_re_i = nmap_re_l[i]
            maps_im_i = nmap_im_l[i]

            for j in 1:s_l, k in j:s_l
                nmap_re = maps_re_i[j][k]
                nmap_im = maps_im_i[j][k]
                isempty(nmap_re) && isempty(nmap_im) && continue

                Deltag_re = 0//1
                Deltag_im = 0//1

                # Deltag_re from n_re
                for (t, n_re) in nmap_re
                    Deltag_re += Delta[t] * n_re
                end
                # Deltag_im from n_im
                for (t, n_im) in nmap_im
                    Deltag_im += Delta[t] * n_im
                end

                (Deltag_re == 0//1 && Deltag_im == 0//1) && continue

                if l == 1 && i == 1
                    old = Gi[j+1, k+1]
                    re_old = getfield(old, :re)
                    im_old = getfield(old, :im)
                    new_val = Complex{Rat}(re_old + Deltag_re, im_old + Deltag_im)
                    Gi[j+1, k+1] = new_val
                    (j != k) && (Gi[k+1, j+1] = conj(new_val))
                else
                    old = Gi[j, k]
                    re_old = getfield(old, :re)
                    im_old = getfield(old, :im)
                    new_val = Complex{Rat}(re_old + Deltag_re, im_old + Deltag_im)
                    Gi[j, k] = new_val
                    (j != k) && (Gi[k, j] = conj(new_val))
                end
            end
        end
    end

    return G1_blocks, G2_blocks
end

function apply_A_to_grams_rat(
    data,
    Nsites::Int,
    G1_blocks::Vector{<:AbstractMatrix{Complex{Rat}}},
    G2_blocks::Vector{<:AbstractMatrix{Complex{Rat}}},
    gram_maps::GramNMaps;
    tol_dft::Real = 1e-20,
)
    tsupp = data.tsupp
    T     = length(tsupp)
    RHS   = zeros(Rat, T)

    s    = Int.(length.(data.basis) ./ Nsites)
    s1, s2 = s
    Imax = Int(Nsites ÷ 2) + 1
    two  = 2//1
    sqrtN = sqrtN_rat(Nsites; tol_dft = tol_dft)

    # empty word index
    empty_idx = bfind(tsupp, UInt16[])
    @assert empty_idx !== nothing && empty_idx != 0
    t_empty = empty_idx

    # const/const
    g00 = G1_blocks[1][1,1]
    g00_re = getfield(g00, :re)
    RHS[t_empty] += g00_re

    @inbounds for k in 1:s1
        wi = data.basis[1][Nsites*(k-1) + 1]
        w_nf, c_nf = QMBCertify.reduce!(wi; L = Nsites)
        c_nf == 0 && continue
        loc = bfind(tsupp, w_nf)
        (loc === nothing || loc == 0) && continue

        g = G1_blocks[1][1, k+1]
        g_re = getfield(g, :re)
        coeff = two * sqrtN
        RHS[loc] += coeff * g_re
    end

    nmap_re_1 = gram_maps.nmap_re_1
    nmap_im_1 = gram_maps.nmap_im_1
    nmap_re_2 = gram_maps.nmap_re_2
    nmap_im_2 = gram_maps.nmap_im_2

    @inbounds for l in 1:2
        G_blocks = (l == 1 ? G1_blocks : G2_blocks)
        s_l      = (l == 1 ? s1 : s2)
        nmap_re_l = (l == 1 ? nmap_re_1 : nmap_re_2)
        nmap_im_l = (l == 1 ? nmap_im_1 : nmap_im_2)

        for i in 1:Imax
            Gi = G_blocks[i]
            maps_re_i = nmap_re_l[i]
            maps_im_i = nmap_im_l[i]

            for j in 1:s_l, k in j:s_l
                g = if l == 1 && i == 1
                    Gi[j+1, k+1]
                else
                    Gi[j,   k]
                end
                g_re = getfield(g, :re)
                g_im = getfield(g, :im)
                iszero(g_re) && iszero(g_im) && continue

                nmap_re = maps_re_i[j][k]
                nmap_im = maps_im_i[j][k]
                isempty(nmap_re) && isempty(nmap_im) && continue

                for (t, n_re) in nmap_re
                    RHS[t] += n_re * g_re
                end
                for (t, n_im) in nmap_im
                    RHS[t] += n_im * g_im
                end
            end
        end
    end

    return RHS
end

function frob(Gnum::AbstractMatrix{ComplexF64},
                  Grat::AbstractMatrix{Complex{Rational{BigInt}}};
                  precision::Integer=256)

    @assert size(Gnum) == size(Grat) "Size mismatch: $(size(Gnum)) vs $(size(Grat))"
    setprecision(precision) do
        # promote to Complex{BigFloat}
        GnumBF = Complex{BigFloat}.(big.(real(Gnum)), big.(imag(Gnum)))
        GratBF = Complex{BigFloat}.(BigFloat.(real(Grat)), BigFloat.(imag(Grat)))

        Delta = GnumBF .- GratBF

        abs2_sum(A) = foldl((s,x)->s + abs2(x), A; init=zero(BigFloat))
        nf  = sqrt(abs2_sum(Delta))
        nref = sqrt(abs2_sum(GnumBF))

        return (Float64(nf), Float64(nf / (nref == 0 ? one(BigFloat) : nref)))
    end
end

function _round_hermitian_to_rational(G::AbstractMatrix{ComplexF64}, tol::Real)
    n = size(G, 1)
    
    R = Matrix{Complex{Rational{BigInt}}}(undef, n, n)
    zero_val = Complex{Rational{BigInt}}(0//1, 0//1)

    @inbounds for j in 1:n
        for i in 1:j
            z = G[i, j]

            if z == (0.0 + 0.0im)
                val = zero_val
                R[i, j] = val
                R[j, i] = val
            else
                re_rat = Rational{BigInt}(rationalize(real(z); tol = tol))
                im_rat = Rational{BigInt}(rationalize(imag(z); tol = tol))
                val    = Complex{Rational{BigInt}}(re_rat, im_rat)
                R[i, j] = val
                R[j, i] = Complex{Rational{BigInt}}(re_rat, -im_rat)
            end
        end
    end
    return R
end

function round_project_qmb(data, Nsites, J, numopt;
                                 tol_gram::Real = 1e-12,
                                 tol_dft::Real  = 1e-20,
                                 snn::Bool=false, J2::Real=1.0,
                                 check::Bool=false)

    t_p = @elapsed begin
        vars      = make_vars_chain(Nsites)
        Hpoly     = heisenberg_chain_poly(J, vars, Nsites; snn=snn, J2=J2)
        LHS_poly  = poly_nf_QMBC(Hpoly - numopt, vars; L=Nsites, realify=true)

        LHS_rat_p = rationalize_poly(LHS_poly; tol=0)

        LHS_vec   = coeff_vector_LHS_tsupp_rat(LHS_rat_p, vars, data.tsupp; N=Nsites)
    end
    println("Time to build and rationalize LHS vector: ", t_p, " seconds")

    t_r = @elapsed begin
        G1_blocks, G2_blocks = round_grams(data, Nsites; tol_gram=tol_gram)
    end
    println("Time to round small Gram blocks: ", t_r, " seconds")

    G1_blocks_before = deepcopy(G1_blocks)
    G2_blocks_before = deepcopy(G2_blocks)

    t_rhs0 = @elapsed begin
        RHS_before, N_before, gram_maps =
            build_rhs_and_N_rat(data, Nsites, G1_blocks, G2_blocks; tol_dft = tol_dft)

        resid_before = RHS_before .- LHS_vec
    end
    println("Time to build RHS & N (pre-projection): ", t_rhs0, " seconds")


    t_proj = @elapsed begin
        G1_blocks, G2_blocks = project_blocks_frob!(
            data, Nsites,
            G1_blocks, G2_blocks,
            LHS_vec, RHS_before, N_before, gram_maps;
            tol_dft = tol_dft,
        )
    end
    println("Time to project blocks: ", t_proj, " seconds")

    if check
        t_rhs = @elapsed begin
            RHS_chk = apply_A_to_grams_rat(
                data, Nsites,
                G1_blocks, G2_blocks,
                gram_maps;
                tol_dft = tol_dft,
            )
            resid = RHS_chk .- LHS_vec

            println("CHECK  ||RHS-LHS||inf = ", Float64(maximum(abs.(resid))))
            println("CHECK  ||RHS-LHS||1 = ", Float64(sum(abs.(resid))))
        end
        println("Time to rebuild RHS for check: ", t_rhs, " seconds")
    end

    return (G1_blocks_proj    = G1_blocks,
            G2_blocks_proj    = G2_blocks,
            G1_blocks_before  = G1_blocks_before,
            G2_blocks_before  = G2_blocks_before,
            LHS_rat           = LHS_vec)
end


function round_grams(data, Nsites::Int; tol_gram::Real = 1e-12)

    G1_num_blocks = data.GramMat[1] 
    G2_num_blocks = data.GramMat[2]

    n1 = length(G1_num_blocks)
    n2 = length(G2_num_blocks)

    G1_blocks = Vector{Matrix{Complex{Rational{BigInt}}}}(undef, n1)
    G2_blocks = Vector{Matrix{Complex{Rational{BigInt}}}}(undef, n2)

    for i in 1:n1
        Gi = G1_num_blocks[i]

        Gc = ComplexF64.(Gi)
        G1_blocks[i] = _round_hermitian_to_rational(Gc, tol_gram)
    end
    for i in 1:n2
        Gi = G2_num_blocks[i]
        Gc = ComplexF64.(Gi)
        G2_blocks[i] = _round_hermitian_to_rational(Gc, tol_gram)
    end

    return G1_blocks, G2_blocks
end

tol_from_digits(d::Integer) = 10.0^(-(d-1))

function dmrg_heisenberg_rat(N, J; J2=0.0, digits = 12)
    J2 = Float64(J2)
    sites = siteinds("S=1/2", N)

    os = OpSum()
    for j = 1:N-1
        os += "Sx",j,"Sx",j+1
        os += "Sy",j,"Sy",j+1
        os += "Sz",j,"Sz",j+1
    end
    for j = 1:N-2
        os += J2,"Sx",j,"Sx",j+2
        os += J2,"Sy",j,"Sy",j+2
        os += J2,"Sz",j,"Sz",j+2
    end
    os += "Sx",1,"Sx",N
    os += "Sy",1,"Sy",N
    os += "Sz",1,"Sz",N
    os += J2,"Sx",1,"Sx",N-1
    os += J2,"Sy",1,"Sy",N-1
    os += J2,"Sz",1,"Sz",N-1
    os += J2,"Sx",2,"Sx",N
    os += J2,"Sy",2,"Sy",N
    os += J2,"Sz",2,"Sz",N

    H = MPO(os, sites)
    nsweeps = 4
    maxdim = [10, 20, 100, 180]
    cutoff = [1E-7]
    psi0 = randomMPS(sites, 2)

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
    energy_per_site = energy / N

    return rationalize(ceil(energy_per_site, digits = digits),
                       tol = tol_from_digits(digits))
end
