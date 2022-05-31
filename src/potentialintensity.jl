"""
    get_buoyancy_of_lifted_parcel(tparcel <: Real,rparcel <: Real,pparcel <: Real,t :: Array{<: Real},r :: Array{<: Real},p  :: Array{<: Real},ptop=50)
Get the buoyancy profile defined for each level `z` as the difference between the temperature of a parcel lifted from the level `p = pparcel` to the level `p(z)`.
    `Buoyancy(z) = Tv_lifted(z) - Tv_env(z)` 

This computes the lifting condensation level (LCL) to decide if the parcel should be lifted adiabatically or not. If not adiabatically, it uses the Newton-Raphson method to find the temperature at a level while conserving specific entropy.
    `Buoyancy(z) = Tv_lifted(z) - Tv_env(z)` 
"""
function get_buoyancy_of_lifted_parcel(tparcel, rparcel, pparcel, t, r, p, ptop=59*u"hPa")
    n_valid_levels = findfirst(<(ptop),p)
    p = p[begin:n_valid_levels]
    t = t[begin:n_valid_levels]
    r = r[begin:n_valid_levels]
    tvirtual_diff_parcel_env = similar(t)
    parcel_sat_vapor_pressure = get_saturation_vapor_pressure(tparcel)
    parcel_vapor_pressure = get_partial_vapor_pressure(rparcel,pparcel)
    parcel_rh = min(parcel_vapor_pressure/parcel_sat_vapor_pressure  , 1.0)
    parcel_specific_entropy = get_specific_entropy(tparcel,rparcel,pparcel; adjust_for_ice_phase = true)
    parcel_lcl = get_lifted_condensation_level(tparcel,parcel_rh,pparcel)
    levels_below_lcl = findall(>=(parcel_lcl),p)
    levels_above_lcl = findall(pres -> ptop < pres < parcel_lcl,p)
    
    for level in levels_below_lcl
        tlifted = tparcel*(p[level]/pparcel)^(Dryair.R/Dryair.cp)
        tvirtual_lifted = get_virtual_temperature(tlifted,rparcel,rparcel)
        tvirtual_env = get_virtual_temperature(t[level],r[level],r[level])
        tvirtual_diff_parcel_env[level] = tvirtual_lifted - tvirtual_env
    end

    for level in levels_above_lcl
        initial_guess = t[level]
        target_value = parcel_specific_entropy
        tlifted = find_root_newton_raphson(temp -> get_specific_entropy_emanuel(temp,rparcel,p[level]), temp -> ∂specific_entropy_∂temp_emanuel(temp,rparcel, p[level]); target_value, initial_guess)
        saturation_vapor_pressure_lifted = get_saturation_vapor_pressure(tlifted)
        mixing_ratio_lifted = get_mixing_ratio(saturation_vapor_pressure_lifted,p[level])
        tvirtual_lifted = get_virtual_temperature(tlifted, rparcel, mixing_ratio_lifted)
        tvirtual_env = get_virtual_temperature(t[level],r[level], r[level])
        tvirtual_diff_parcel_env[level] = tvirtual_lifted - tvirtual_env
    end
    return tvirtual_diff_parcel_env
end

function get_buoyancy_of_lifted_parcel(tparcel :: Real, rparcel :: Real ,pparcel :: Real, t :: Array{ <: Real} ,r :: Array{ <: Real},p :: Array{ <: Real}, ptop=59)
    ustrip.(get_buoyancy_of_lifted_parcel(1u"K" * tparcel , 1u"kg/kg" * rparcel ,1u"hPa" * pparcel, 1u"K" .* t , 1u"kg/kg" .* r, 1u"hPa" .* p, 1u"hPa" * ptop ))
end


"""
    get_potential_intensity_of_tropical_cyclone(sea_surface_temp <: Real,
    sea_surface_pressure <: Real, 
    pressure :: Array{<: Real}, 
    temperature :: Array{<: Real},
     mixing_ratio :: Array{<: Real}; 
    ckovercd = 0.9, reversible_ascent=1, dissipative_heating = true)

Compute the minimum pressure at the center and the maximum wind speed of a tropical cyclone using Emanuel's potential intensity theory.
"""
function get_potential_intensity_of_tropical_cyclone(sea_surface_temperature,sea_surface_pressure, pressure, temperature, mixing_ratio; ck_over_cd = 0.9, reversible_ascent=true, dissipative_heating = true, vreduc = 0.8)

    initial_level_for_lifting = 1
    exponent_central_pressure = 2.0
    saturation_vapor_pressure_surface = get_saturation_vapor_pressure(sea_surface_temperature)
#
#   ***   Find environmental CAPE *** 
#
    tparcel=temperature[initial_level_for_lifting];
    rparcel=mixing_ratio[initial_level_for_lifting]
    pparcel=pressure[initial_level_for_lifting]
    cape_env, ~ , index_level_of_neutral_buoyancy = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,temperature,mixing_ratio,pressure)

    pressure_at_rmax = 950.0u"hPa"
    pressure_at_rmax_old = 0.0u"hPa"
    saturation_cape_at_rmax = 0.0u"J/kg"
    cape_at_rmax = 0.0u"J/kg"
    average_virtual_temp = 0u"K"
    temp_ratio = 0.0
    niter = 1
    mixing_ratio_lowest_level = mixing_ratio[initial_level_for_lifting]
    temperature_lowest_level = temperature[initial_level_for_lifting]
    virtual_temp_lowest_level = get_virtual_temperature(temperature_lowest_level, mixing_ratio_lowest_level, mixing_ratio_lowest_level)
    #@show cape_env
    while (abs(pressure_at_rmax_old-pressure_at_rmax)) > 0.2u"hPa" && niter < 200
     
        #These three is where the iteration happens
        pparcel_approx=min(pressure_at_rmax,1000.0u"hPa") #these two are the ones we are iterating over
        rparcel_approx = ϵ*mixing_ratio_lowest_level*sea_surface_pressure / (pparcel_approx*(ϵ+mixing_ratio_lowest_level) - mixing_ratio_lowest_level*sea_surface_pressure) #what in the name of god is this? it is not documented
        rparcel_sat=get_mixing_ratio(saturation_vapor_pressure_surface,pparcel_approx)
        #   ***  Find CAPE at radius of maximum winds   ***             
        cape_at_rmax, ~ , index_level_of_neutral_buoyancy = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel_approx,pparcel_approx,temperature,mixing_ratio,pressure)
        #  ***  Find saturation CAPE at radius of maximum winds   ***
        saturation_cape_at_rmax, temp_outflow, index_lnb = 
        get_cape_and_outflow_temp_from_sounding(sea_surface_temperature,
            get_mixing_ratio(saturation_vapor_pressure_surface,pparcel_approx),
            pparcel_approx,
            temperature,
            mixing_ratio,
            pressure)
        #  ***  Estimate of pressure at radius of maximum winds  ***
        virtual_temp_parcel_sst=get_virtual_temperature(sea_surface_temperature, rparcel_sat, rparcel_sat)
        average_virtual_temp = 0.5 * (virtual_temp_lowest_level + virtual_temp_parcel_sst)


        temp_ratio = dissipative_heating ? sea_surface_temperature/temp_outflow : 1.0
        CAT=cape_at_rmax-cape_env + 0.5 * ck_over_cd * temp_ratio *(saturation_cape_at_rmax - cape_at_rmax)
        CAT=max(CAT,0.0u"J/kg")
        
        pressure_at_rmax_old = pressure_at_rmax
        pressure_at_rmax = sea_surface_pressure*exp(-CAT / (Dryair.R * average_virtual_temp) )

        niter = niter + 1
    end

    reduction_factor=0.5(1.0 + 1.0/exponent_central_pressure)
    CAT=(cape_at_rmax-cape_env)+reduction_factor*ck_over_cd*temp_ratio*(saturation_cape_at_rmax-cape_at_rmax)
    #@info CAT
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

function get_potential_intensity_of_tropical_cyclone(sea_surface_temperature :: Real,sea_surface_pressure :: Real, pressure :: Array{<: Real}, temperature :: Array{<: Real}, mixing_ratio :: Array{<: Real}; ck_over_cd = 0.9, reversible_ascent=true, dissipative_heating = true, vreduc = 0.8)
    return ustrip.(get_potential_intensity_of_tropical_cyclone(1u"K" * sea_surface_temperature, 1u"hPa" * sea_surface_pressure, 1u"hPa" .* pressure, 1u"K" .* temperature, 1u"kg/kg" .* mixing_ratio; ck_over_cd, reversible_ascent, dissipative_heating, vreduc))
end

"""
    get_cape_and_outflow_temp_from_sounding(tparcel :: Real, 
    rparcel :: Real, 
    pparcel :: Real,
    t Array{<: Real},
    r Array{<: Real},
    p Array{<: Real}, ptop=50u"hPa")
Compute cape, outflow temperature and index of neutral buoyancy from thermodynamic profiles.
"""
function get_cape_and_outflow_temp_from_sounding(tparcel, rparcel, pparcel, t, r, p, ptop=50u"hPa")
    buoyancy_profile = get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop)
    negative_area=0.0*unit(buoyancy_profile[1]*Dryair.R)
    positive_area=0.0*unit(buoyancy_profile[1]*Dryair.R)
    level_neutral_buoyancy = searchsortedlast(buoyancy_profile, 0.0u"K"; rev=true)
    for level in 2:(level_neutral_buoyancy + 1)
        area=Dryair.R*(buoyancy_profile[level]+buoyancy_profile[level-1])*(p[level-1]-p[level])/(p[level]+p[level-1])
        positive_area += max(area,0.0u"J/kg")
        negative_area -= min(area,0.0u"J/kg")
    end
    outflow_temp = t[level_neutral_buoyancy]
    #Add buoyancy of parcel with respect to first level? (lcl may not be a level in the profile)
    parcel_buoyancy_area = Dryair.R *(pparcel - p[1])/(pparcel + p[1])
    positive_area += parcel_buoyancy_area*max(buoyancy_profile[1],0.0*unit(buoyancy_profile[1]))
    negative_area -= parcel_buoyancy_area*min(buoyancy_profile[1],0.0*unit(buoyancy_profile[1]))
    #Add residual above inb and t0
    ## This is unsafe, need to check for bounds
    pres_neutral_buoyancy = (p[level_neutral_buoyancy + 1]*buoyancy_profile[level_neutral_buoyancy] - p[level_neutral_buoyancy]*buoyancy_profile[level_neutral_buoyancy + 1]) / (buoyancy_profile[level_neutral_buoyancy] - buoyancy_profile[level_neutral_buoyancy + 1]) 
    residual_area = Dryair.R * buoyancy_profile[level_neutral_buoyancy]*(p[level_neutral_buoyancy] - pres_neutral_buoyancy)/(p[level_neutral_buoyancy] + pres_neutral_buoyancy)
    
    outflow_temp = (outflow_temp * (pres_neutral_buoyancy - p[level_neutral_buoyancy + 1]) + t[level_neutral_buoyancy + 1] * (p[level_neutral_buoyancy] - pres_neutral_buoyancy)) / (p[level_neutral_buoyancy] - p[level_neutral_buoyancy + 1])

    return positive_area - negative_area + residual_area, outflow_temp, level_neutral_buoyancy
end


function get_cape_and_outflow_temp_from_sounding(tparcel :: Real,rparcel :: Real,pparcel :: Real,t :: Array{<: Real},r :: Array{<: Real},p :: Array{<: Real},ptop=50)
    ustrip.(get_cape_and_outflow_temp_from_sounding(1u"K" * tparcel,1u"kg/kg" * rparcel,1u"hPa" * pparcel, 1u"K" .* t, 1u"kg/kg" .*r,1u"hPa" .* p,1u"hPa" * ptop))
end