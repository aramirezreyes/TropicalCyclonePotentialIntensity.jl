using TropicalCyclonePotentialIntensity
using Test

using Unitful: @u_str, unit, ustrip, Quantity
using NCDatasets: Dataset, variable

@testset "TropicalCyclonePotentialIntensity.jl" begin
    include("physicsfunctions.jl")
    include("potentialintensity.jl")    # Write your tests here.
    
end
