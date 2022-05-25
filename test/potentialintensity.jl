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
    @info AvailablePotentialEnergyFramework.get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,tabs[:,timeindex],r[:,timeindex],pres[:,timeindex])
    include(joinpath(@__DIR__,"emanuel_potential_intensity_wrapper.jl"))
    @info get_cape(ustrip(tparcel),ustrip(rparcel), ustrip(pparcel), ustrip.(tabs[:,timeindex]), ustrip.(r[:,timeindex]), ustrip.(pres[:,timeindex]))
