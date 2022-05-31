module TropicalCyclonePotentialIntensity

# Write your package code here.
using Unitful: @u_str, unit, ustrip, Quantity

export mixing_ratio_to_specific_humidity,
    specific_humidity_to_mixing_ratio,
    get_saturation_vapor_pressure,
    get_partial_vapor_pressure,
    get_mixing_ratio,
    get_specific_entropy,
    get_lifted_condensation_level,
    get_potential_temperature,
    get_virtual_temperature,
###Potential intensity
    get_buoyancy_of_lifted_parcel,
    get_cape_and_outflow_temp_from_sounding,
    get_potential_intensity_of_tropical_cyclone

include("physicalconstants.jl")
include("physicsfunctions.jl")
include("rootfinding.jl")
include("potentialintensity.jl")

const ϵ = Dryair.R / Watervapor.R


end
