pres = Dataset(joinpath(@__DIR__,"testfiles/thermoprofile.nc")) do ds 
    1u"hPa" .* variable(ds, "PRES")[:,:]
end
tabs = Dataset(joinpath(@__DIR__,"testfiles/thermoprofile.nc")) do ds 
    1u"K" .* variable(ds, "TABS")[:,:]
end
qv = 1f-3u"kg/g" .* 1u"g/kg" .* Dataset(joinpath(@__DIR__,"testfiles/thermoprofile.nc")) do ds 
    variable(ds, "QV")[:,:] #was originally in g/kg
end
@info size(pres) size(qv) size(tabs)
r = specific_humidity_to_mixing_ratio.(qv)
timeindex = 1200
pparcel = pres[1,timeindex]
tparcel = tabs[1,timeindex]
rparcel = r[1,timeindex]

#I will create a similar profile but with a perturbation to see what happens

@test unit(get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,tabs[:,timeindex],r[:,timeindex],pres[:,timeindex])[1]) == u"K"

cape_this_implementation, temp_outflow_this_implementation, index_outflow_this_implementation = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,tabs[:,timeindex],r[:,timeindex],pres[:,timeindex])


#include(joinpath(@__DIR__,"emanuel_potential_intensity_wrapper.jl"))
#cape_fortan_implementation, temp_outflow_fortran_implementation = get_cape(ustrip(tparcel),ustrip(rparcel), ustrip(pparcel), ustrip.(tabs[:,timeindex]), ustrip.(r[:,timeindex]), ustrip.(pres[:,timeindex]))

#@test_broken isapprox(cape_this_implementation, cape_fortran_implementation*u"J/kg", rtol = 0.01)
#@test_broken isapprox(temp_outflow_this_implementation, temp_outflow_fortran_implementation*u"K", rtol = 0.01)

#min_pres_this_implementation, max_speed_this_implementation = get_minimum_pressure_of_tropical_cyclone(tparcel, pparcel, pres[:,timeindex], tabs[:,timeindex], r[:,timeindex])
#min_pres_fortran_implementation, max_speed_fortran_implementation = get_pcmin( ustrip(tparcel) .- 273.15f0,ustrip(pparcel),ustrip.(pres[:,timeindex]),ustrip.(tabs[:,timeindex]) .- 273.15f0, 1f3.*ustrip.(r[:,timeindex]) ) 
#@test isapprox(min_pres_this_implementation, min_pres_fortran_implementation*u"hPa", rtol = 0.01)
#@test_broken isapprox(max_speed_this_implementation, max_speed_fortran_implementation*u"m/s", rtol = 0.01)