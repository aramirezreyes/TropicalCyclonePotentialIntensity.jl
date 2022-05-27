using TropicalCyclonePotentialIntensity
using Test

using Unitful: @u_str, unit, ustrip, Quantity
using DelimitedFiles

@testset "TropicalCyclonePotentialIntensity.jl" begin
    include("physicsfunctions.jl")
    include("potentialintensity.jl")   
end
