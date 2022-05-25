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
    parcel_specific_entropy = get_specific_entropy(tparcel,rparcel,pparcel)
    @info parcel_specific_entropy
    @info parcel_vapor_pressure
    @info parcel_sat_vapor_pressure
    parcel_lcl = get_lifted_condensation_level(tparcel,parcel_rh,pparcel)
    @info parcel_lcl
#    @show parcel_lcl
    levels_below_lcl = findall(>=(parcel_lcl),p)
    levels_above_lcl = findall(<(parcel_lcl),p)
    
    #These two must populate buoyancy of lifted parcel_get_vapor_pressure
    #this would be adiabatic lifting, easy enough
    for level in levels_below_lcl
        tlifted = tparcel*(p[level]/pparcel)^(Dryair.R/Dryair.cp)
        rlifted = rparcel
        tvirtual_lifted = get_virtual_temperature(tlifted,rparcel,rparcel)
        tvirtual_env = get_virtual_temperature(t[level],r[level],r[level])
        tvirtual_diff_parcel_env[level] = tvirtual_lifted - tvirtual_env
        #@info level, tvirtual_diff_parcel_env[level]
    end

    #We start with environmental values of temperature, mixing ratio, entropy etc
    #Our goal: to find the temperature Tx such that s_approx(Tx) ≈ s_parcel s being total entropy
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
        @info level, tvirtual_diff_parcel_env[level]
    end
    return tvirtual_diff_parcel_env
end



"""
    get_minimum_pressure_of_tropical_cyclone(sea_surface_temp,sea_surface_pressure, pressure, temperature, mixing_ratio; ckovercd = 0.9, reversible_ascent=1, dissipative_heating = true)
temperatures in kelvin

"""
function get_minimum_pressure_of_tropical_cyclone(sea_surface_temp,sea_surface_pressure, pressure, temperature, mixing_ratio; ckovercd = 0.9, reversible_ascent=1, dissipative_heating = true)

    initial_level_for_lifting = 1
    exponent_central_pressure = 2.0
    vreduc = 0.8

    t0 = 230.0
    vmax = 0.0
    pmin = 0.0
    ifl = 0 #what is this
    saturation_vapor_pressure0=6.112.*exp(17.67*sea_surface_temperature/(243.5+sea_surface_temperature));
    # initial values (for what?)
    ifl=1;
    NP=0;
    min_pressure=970.0;
    min_pressure_old=min_pressure;
    min_pressure_neew=0.0;
#
#   ***   Find environmental CAPE *** 
#
    tparcel=temperature[initial_level_for_lifting];
    rparcel=mixing_ratio[initial_level_for_lifting]
    pparcel=pressure[initial_level_for_lifting]
    cape_env, outflow_temp_env, index_level_of_neutral_buoyancy = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,temperature,mixing_ratio,pressure)

    while (abs(min_pressure_new-min_pressure_old)) > 0.2
#
#   ***  Find CAPE at radius of maximum winds   ***
#
      tparcel=temperature[initial_level_for_lifting]
      pparcel_approx=min(min_pressure,1000.0) #these two are the ones we are iterating over
      rparcel_approx = epsilon*mixing_ratio[initial_level_for_lifting]*sea_surface_pressure / (pparcel*(epsilon+mixing_ratio[initial_level_for_lifting]) - mixing_ratio[initial_level_for_lifting]*sea_surface_pressure) #what in the name of god is this? it is not documented
      cape_at_rmax, outflow_temp_at_rmax, index_level_of_neutral_buoyancy = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel_approx,pparcel_approx,temperature,mixing_ratio,pressure);
#
#  ***  Find saturation CAPE at radius of maximum winds   ***
#
      tparcel=sea_surface_temperature;
      pparcel_sat=min(min_pressure,1000.0)
      rparcel_sat=get_mixing_ratio(saturation_vapor_presure0,pparcel_approx)
      saturation_cape_at_rmax, temp_outflow, intdex_lnb = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel_sat,pparcel_sat,temperature,mixing_ratio,pressure);
      temp_ratio=sea_surface_temperature/temp_outflow;
       
#
#  ***  Initial estimate of minimum pressure   ***
#
        virtual_temp_parcel_approx=get_virtual_temperature(temperature[initial_level_for_lifting],mixing_ratio_parcel_approx)
        virtual_temp_parcel_sst=get_virtual_temperature(sea_surface_temperature,specific_humidity_parcel_approx)
        average_virtual_temp = 0.5.*(virtual_temperature_parcel_approx + virtual_temperature_parcel_sst);
	CAT=cape_at_rmax-cape_env+0.5.*ck_over_cd.*temp_ratio.*(saturation_cape_at_rmax-cape_at_rmax);
	CAT=max(CAT,0.0);
	min_pressure_new=sea_surface_temperature*exp(-CAT./(287.04.*average_virtual_temp));

	min_pressure_old = min_pressure;
	min_pressure = min_pressure_new;
	niter = niter + 1;
    end
end


function get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,t,r,p,ptop=50*unit(pparcel))
    buoyancy_profile = get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop)
#    @show buoyancy_profile
    negative_area=0.0*unit(buoyancy_profile[1]*Dryair.R)
    positive_area=0.0*unit(buoyancy_profile[1]*Dryair.R)
    level_neutral_buoyancy = searchsortedlast(buoyancy_profile, 0.0unit(buoyancy_profile[1]); rev=true)
    @info level_neutral_buoyancy
    for level in 2:(level_neutral_buoyancy + 1)
        area=Dryair.R*(buoyancy_profile[level]+buoyancy_profile[level-1])*(p[level-1]-p[level])/(p[level]+p[level-1])
        positive_area += max(area,0.0*unit(area))
        negative_area -= min(area,0.0*unit(area))
    end
    outflow_temp = t[level_neutral_buoyancy]
    return positive_area - negative_area, outflow_temp, level_neutral_buoyancy
end
