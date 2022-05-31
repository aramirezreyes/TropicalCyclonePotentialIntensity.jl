pres = 1.0f0u"hPa" .* Float32.(readdlm("pressure_in_hpa")[:])
tabs = 1.0f0u"K" .* Float32.(readdlm("tabs_in_Kelvin")[:])
qv = 1.0f0u"kg/kg" .* Float32.(readdlm("specific_humidity_in_kgperkg")[:])
r = specific_humidity_to_mixing_ratio.(qv)

pparcel = pres[1]
tparcel = tabs[1]
rparcel = r[1]

#With these value, the fortran implementation yields:
include(joinpath(@__DIR__,"emanuel_potential_intensity_wrapper.jl"))
cape_fortran_implementation, temp_outflow_fortran_implementation = get_cape(ustrip(tparcel),ustrip(rparcel), ustrip(pparcel), ustrip.(tabs), ustrip.(r), ustrip.(pres))
@show cape_fortran_implementation, temp_outflow_fortran_implementation 

min_pres_fortran_implementation, max_speed_fortran_implementation = get_pcmin( ustrip(tparcel) .- 273.15f0,ustrip(pparcel),ustrip.(pres),ustrip.(tabs) .- 273.15f0, 1f3.*ustrip.(r)) 
@show min_pres_fortran_implementation, max_speed_fortran_implementation

#cape_fortran_implementation, temp_outflow_fortran_implementation = (1106.4801f0, 240.90143f0)
#min_pres_fortran_implementation, max_speed_fortran_implementation = (991.1094f0, 26.201181f0)
cape, temp_outflow, index_outflow = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,tabs,r,pres)
min_pres, max_speed = get_minimum_pressure_of_tropical_cyclone(tparcel, pparcel, pres, tabs, r)

@testset "Potential Intensity" begin
    @test unit(get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,tabs,r,pres)[1]) == u"K"
    @test isapprox(cape, cape_fortran_implementation*u"J/kg", rtol = 0.01) 
    @test isapprox(temp_outflow, temp_outflow_fortran_implementation*u"K", rtol = 0.01) 
    @test isapprox(min_pres, min_pres_fortran_implementation*u"hPa", rtol = 0.01)
    @test isapprox(max_speed, max_speed_fortran_implementation*u"m/s", rtol = 0.01) 
end