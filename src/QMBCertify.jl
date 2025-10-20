module QMBCertify

using JuMP
using MosekTools
using LinearAlgebra
using Dualization

export GSB, slabel, reduce!, mosek_para

mutable struct qmb_data
    correlation1
    correlation2
    correlation3
    basis
    tsupp
    GramMat
    moment
end

include("basicfunction.jl")
include("rdm_positivity.jl")
include("sdp.jl")

end
