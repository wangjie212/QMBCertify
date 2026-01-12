using QMBCertify
using ITensorMPS
using ITensors

J2s = 0.2:0.2:2.0
e_lbs = Float64[]
e_ubs = Float64[]
C_bounds_upper = Float64[]
C_bounds_lower = Float64[]
shifts_upper = Float64[]
shifts_lower = Float64[]

for J2 in J2s
    println("=== J2 = $J2 ===")
    res = certify_qmb_corr(
        12,
        3,
        3;
        J2 = J2,
        dist = 1,
        extra_E = 3,
        extra_corr = 3,
        QUIET = true,
        tol_gram = 1e-13,
        tol_dft  = 1e-12
    )
    push!(C_bounds_upper, Float64(res.C_rat_upper[1]))
    push!(C_bounds_lower, Float64(res.C_rat_lower[1]))
    push!(shifts_upper, Float64(res.shift_upper[1]))
    push!(shifts_lower, Float64(res.shift_lower[1]))
    push!(e_lbs, Float64(res.e_lb))
    push!(e_ubs, Float64(res.e_ub))
end
using QMBCertify
using ITensorMPS
using ITensors

J2s = 0.2:0.2:2.0
e_lbs = Float64[]
e_ubs = Float64[]
C_bounds_upper = Float64[]
C_bounds_lower = Float64[]
shifts_upper = Float64[]
shifts_lower = Float64[]

for J2 in J2s
    println("=== J2 = $J2 ===")
    res = certify_qmb_corr(
        12,
        3,
        3;
        J2 = J2,
        dist = 1,
        extra_E = 3,
        extra_corr = 3,
        QUIET = true,
        tol_gram = 1e-13,
        tol_dft  = 1e-12
    )
    push!(C_bounds_upper, Float64(res.C_rat_upper[1]))
    push!(C_bounds_lower, Float64(res.C_rat_lower[1]))
    push!(shifts_upper, Float64(res.shift_upper[1]))
    push!(shifts_lower, Float64(res.shift_lower[1]))
    push!(e_lbs, Float64(res.e_lb))
    push!(e_ubs, Float64(res.e_ub))
end
