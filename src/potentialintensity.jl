"""
    get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop=50)
"""
function get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop=50*unit(pparcel))
    n_valid_levels = findfirst(<(ptop),p)
    p = p[begin:n_valid_levels]
    t = t[begin:n_valid_levels]
    r = r[begin:n_valid_levels]
    tvirtual_diff_parcel_env = similar(t)
    parcel_sat_vapor_pressure = get_saturation_vapor_pressure(tparcel)
    parcel_vapor_pressure = get_partial_vapor_pressure(rparcel,pparcel)
    parcel_rh = min(parcel_vapor_pressure/parcel_sat_vapor_pressure  , 1.0)
    parcel_specific_entropy = get_specific_entropy(tparcel,rparcel,pparcel; adjust_for_ice_phase = true)
    #@info parcel_specific_entropy
    #@info parcel_vapor_pressure
    #@info parcel_sat_vapor_pressure
    parcel_lcl = get_lifted_condensation_level(tparcel,parcel_rh,pparcel)
    #@info parcel_lcl
#    @show parcel_lcl
    levels_below_lcl = findall(>=(parcel_lcl),p)
    levels_above_lcl = findall(<(parcel_lcl),p)
    
    for level in levels_below_lcl
        tlifted = tparcel*(p[level]/pparcel)^(Dryair.R/Dryair.cp)
        tvirtual_lifted = get_virtual_temperature(tlifted,rparcel,rparcel)
        tvirtual_env = get_virtual_temperature(t[level],r[level],r[level])
        tvirtual_diff_parcel_env[level] = tvirtual_lifted - tvirtual_env
        #@info level, tvirtual_diff_parcel_env[level]
    end

    for level in levels_above_lcl
        initial_guess = t[level]
        target_value = parcel_specific_entropy
        tlifted = find_root_newton_raphson(temp -> get_specific_entropy_emanuel(temp,rparcel,p[level]), temp -> ∂specific_entropy_∂temp_emanuel(temp,rparcel, p[level]); target_value, initial_guess)
        #@info level, tlifted
        saturation_vapor_pressure_lifted = get_saturation_vapor_pressure(tlifted)
        mixing_ratio_lifted = get_mixing_ratio(saturation_vapor_pressure_lifted,p[level])
        tvirtual_lifted = get_virtual_temperature(tlifted, rparcel, mixing_ratio_lifted)
        tvirtual_env = get_virtual_temperature(t[level],r[level], r[level])
        tvirtual_diff_parcel_env[level] = tvirtual_lifted - tvirtual_env
        #@info level, tvirtual_diff_parcel_env[level]
    end
    return tvirtual_diff_parcel_env
end



"""
    get_minimum_pressure_of_tropical_cyclone(sea_surface_temp,sea_surface_pressure, pressure, temperature, mixing_ratio; ckovercd = 0.9, reversible_ascent=1, dissipative_heating = true)
temperatures in kelvin

"""
function get_minimum_pressure_of_tropical_cyclone(sea_surface_temperature,sea_surface_pressure, pressure, temperature, mixing_ratio; ck_over_cd = 0.9, reversible_ascent=true, dissipative_heating = true, vreduc = 0.8)

    initial_level_for_lifting = 1
    exponent_central_pressure = 2.0
    saturation_vapor_pressure_surface = get_saturation_vapor_pressure(sea_surface_temperature)
#
#   ***   Find environmental CAPE *** 
#
    tparcel=temperature[initial_level_for_lifting];
    rparcel=mixing_ratio[initial_level_for_lifting]
    pparcel=pressure[initial_level_for_lifting]
    cape_env, outflow_temp_env , index_level_of_neutral_buoyancy = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,temperature,mixing_ratio,pressure)

    pressure_at_rmax = 950.0u"hPa"
    pressure_at_rmax_old = 0.0u"hPa"
    saturation_cape_at_rmax = 0.0u"J/kg"
    cape_at_rmax = 0.0u"J/kg"
    average_virtual_temp = 0u"K"
    temp_ratio = 0.0
    niter = 1
    #@info cape_env
    while (abs(pressure_at_rmax_old-pressure_at_rmax)) > 0.2u"hPa" && niter < 200
       #   @info pressure_at_rmax
       #@info cape_at_rmax
       #@info saturation_cape_at_rmax
#   ***  Find CAPE at radius of maximum winds   ***
#       #this one depends on the pressure
        
        pparcel_approx=min(pressure_at_rmax,1000.0u"hPa") #these two are the ones we are iterating over
        rparcel_approx = epsilon*mixing_ratio[initial_level_for_lifting]*sea_surface_pressure / (pparcel*(epsilon+mixing_ratio[initial_level_for_lifting]) - mixing_ratio[initial_level_for_lifting]*sea_surface_pressure) #what in the name of god is this? it is not documented

        cape_at_rmax, outflow_temp_at_rmax, index_level_of_neutral_buoyancy = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel_approx,pparcel_approx,temperature,mixing_ratio,pressure)
        #
        #  ***  Find saturation CAPE at radius of maximum winds   ***
        #
        #this one also depends on the iterated pressure
        rparcel_sat=get_mixing_ratio(saturation_vapor_pressure_surface,pparcel_approx)
        saturation_cape_at_rmax, temp_outflow, index_lnb = get_cape_and_outflow_temp_from_sounding(sea_surface_temperature,rparcel_sat,pparcel_approx,temperature,mixing_ratio,pressure);
        temp_ratio = dissipative_heating ? sea_surface_temperature/temp_outflow : 1.0

        #
        #  ***  Estimate of pressure at radius of maximum winds  ***
        #
        specific_humidity_lowest_level = mixing_ratio_to_specific_humidity(mixing_ratio[initial_level_for_lifting])
        specific_humidity_sat = mixing_ratio_to_specific_humidity(rparcel_sat)
        virtual_temp_lowest_level=get_virtual_temperature(temperature[initial_level_for_lifting],specific_humidity_lowest_level , specific_humidity_lowest_level )
        virtual_temp_parcel_sst=get_virtual_temperature(sea_surface_temperature,specific_humidity_sat, specific_humidity_sat)
        average_virtual_temp = 0.5 * (virtual_temp_lowest_level + virtual_temp_parcel_sst)
        CAT=cape_at_rmax-cape_env + 0.5 * ck_over_cd * temp_ratio *(saturation_cape_at_rmax - cape_at_rmax)
        CAT=max(CAT,0.0u"J/kg")
        
        pressure_at_rmax_old = pressure_at_rmax
        pressure_at_rmax = sea_surface_pressure*exp(-CAT / (Dryair.R * average_virtual_temp) )

        niter = niter + 1
    end

    reduction_factor=0.5(1.0 + 1.0/exponent_central_pressure)
    CAT=(cape_at_rmax-cape_env)+reduction_factor*ck_over_cd*temp_ratio*(saturation_cape_at_rmax-cape_at_rmax)
    CAT=max(CAT,0.0u"J/kg")
    
    # Calculate the minimum pressure at the eye of the storm
    # BE02 EQN. 4
    min_pressure_at_center = sea_surface_pressure*exp(-CAT/(Dryair.R*average_virtual_temp))
                 
    # Calculate the potential intensity at the radius of maximum winds
    # BE02 EQN. 3, reduced by some fraction (default 20%) to account for the reduction 
    # of 10-m winds from gradient wind speeds (Emanuel 2000, Powell 1980)
    vmax=vreduc*sqrt(ck_over_cd*temp_ratio*(saturation_cape_at_rmax-cape_at_rmax))

    return min_pressure_at_center, vmax
end


function get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,t,r,p,ptop=50*unit(pparcel))
    buoyancy_profile = get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop)
#    @show buoyancy_profile
    negative_area=0.0*unit(buoyancy_profile[1]*Dryair.R)
    positive_area=0.0*unit(buoyancy_profile[1]*Dryair.R)
    level_neutral_buoyancy = searchsortedlast(buoyancy_profile, 0.0unit(buoyancy_profile[1]); rev=true)
    #@info level_neutral_buoyancy
    for level in 2:(level_neutral_buoyancy + 1)
        area=Dryair.R*(buoyancy_profile[level]+buoyancy_profile[level-1])*(p[level-1]-p[level])/(p[level]+p[level-1])
        positive_area += max(area,0.0*unit(area))
        negative_area -= min(area,0.0*unit(area))
    end
    outflow_temp = t[level_neutral_buoyancy]
    return positive_area - negative_area, outflow_temp, level_neutral_buoyancy
end
