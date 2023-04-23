module ManyBodySOS

using JuMP
using MosekTools
using LinearAlgebra

export GSE, slabel

include("basicfunction.jl")
include("sos.jl")
include("positivity.jl")

end
