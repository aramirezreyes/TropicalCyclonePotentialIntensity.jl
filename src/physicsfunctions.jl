

"""
    get_potential_temperature(temperature, pressure, reference_pressure)
Compute potential temperature from temperature and pressure.
"""
function get_potential_temperature(temperature, pressure, reference_pressure)
    exponent = ustrip(Dryair.R / Dryair.cp)
    return temperature * (reference_pressure/pressure)^exponent
end

function get_potential_temperature(temperature :: Quantity, pressure :: Quantity, reference_pressure :: Quantity)
    exponent = Dryair.R / Dryair.cp
    return temperature * (reference_pressure/pressure)^exponent
end

"""
    get_virtual_temperature(temperature, specific_humidity)
Compute virtual temperature from temperature and specific humidity.
"""
function get_virtual_temperature(temperature, specific_humidity)
    return (one(temperature) + one(temperature)/1000*epsilon*specific_humidity)*temperature
end

function get_virtual_temperature(temperature :: Quantity, specific_humidity :: Quantity)
    return (one(temperature) + one(temperature)/1000*epsilon*specific_humidity)*temperature
end

"""
    get_virtual_temperature(temperature,mixing_ratio_total_water,mixing_ratio_water_vapor)
Receive temperature (K) and mixing ratios of total water and water vapor (unitless g/g) and compute the virtual temperature
"""
function get_virtual_temperature(temperature,mixing_ratio_total_water,mixing_ratio_water_vapor)
    return temperature*(1 + mixing_ratio_water_vapor/epsilon)/(1 + mixing_ratio_total_water)
end

"""
    specific_humidity_to_mixing_ratio(specific_humidity)
Take a specific humidity (unitless g/g) and return a mixing ratio
"""
function specific_humidity_to_mixing_ratio(specific_humidity)
return mixing_ratio = specific_humidity / (1 - specific_humidity)

end

"""
    mixing_ratio_to_specific_humidity(mixing_ratio)
Take a mixing ratio (unitless g/g) and return a specific humidity
"""
function mixing_ratio_to_specific_humidity(mixing_ratio)
    return q = mixing_ratio / (1 + mixing_ratio)
end



"""
    get_saturation_vapor_pressure(T)
Receive temperature T in Kelvin and compute the saturation vapor pressure in hPa from the August-Roche-Magnus formula that approximates the solution to the Clausius-Clapeyron relationship (Wikipedia contributors. (2020, December 19). Clausius–Clapeyron relation. In Wikipedia, The Free Encyclopedia. Retrieved 06:57, December 20, 2020, from https://en.wikipedia.org/w/index.php?title=Clausius%E2%80%93Clapeyron_relation&oldid=995159175)
"""
function get_saturation_vapor_pressure(T)
    return 6.112*exp(17.67 * (T-273.15) / (243.5 + (T - 273.15)))
end

function get_saturation_vapor_pressure(T :: Quantity)
    return 6.112u"hPa"*exp(17.67 * (T-273.15u"K") / (243.5u"K" + (T - 273.15u"K")))
end


"""
    get_partial_vapor_pressure(mixing_ratio,pressure)
Receive a water vapor mixing ratio (unitless g/g) and environmental pressure and compute the partial pressure of water vapor in the same units as the input pressure.
"""
function get_partial_vapor_pressure(mixing_ratio,pressure)
    return mixing_ratio*pressure/(epsilon + mixing_ratio)
end

"""
    get_mixing_ratio(water_vapor_partial_pressure,env_pressure)
Receive a water vapor mixing ratio (unitless g/g) and environmental pressure and compute the partial pressure of water vapor in the same units as the incoming pressure.
"""
function get_mixing_ratio(water_vapor_partial_pressure,env_pressure)
    return epsilon*water_vapor_partial_pressure/(env_pressure - water_vapor_partial_pressure)
end

"""
    get_specific_entropy(temperature,mixing_ratio,pressure)
Receive temperature in Kelvin, water vapor mixing ratio (unitless g/g) and pressure (hPa) and compute the specific entropy of a parcel using equation in Emmanuel's (E94, EQN. 4.5.9)
"""
function get_specific_entropy(temperature,mixing_ratio,pressure ; adjust_for_ice_phase = false)
    # Adjust for ice phase is a modified value of liquid water specific heat capacity to compensate for the lack of explicit ice phase when lifting a parcel(Personal communication with Kerry Emanuel on April 22 2022 - Argel Ramirez Reyes
    alv = Liquidwater.Lv + (Watervapor.cp - Dryair.cp)*(temperature - 273.15f0u"K")
    adjusted_cl = adjust_for_ice_phase ? Liquidwater.cp - 1690u"J/kg/K" : Liquidwater.cp
    vapor_pressure = get_partial_vapor_pressure(mixing_ratio,pressure)
    saturation_vapor_pressure = get_saturation_vapor_pressure(temperature)
    RH = min(vapor_pressure/saturation_vapor_pressure,1.0)
    specific_entropy =  (Dryair.cp + mixing_ratio * adjusted_cl) *
        log(temperature/unit(temperature)) - Dryair.R * log((pressure - vapor_pressure)/unit(pressure)) +
        alv * mixing_ratio / temperature - mixing_ratio * Watervapor.R * log(RH)
end



"""
    get_lifted_condensation_level(temperature,relative_humidity,pressure)   
Receive temperature in Kelvin, relative humidity (unitless) and pressure (hPa) and compute the lifted condensation level based on Emanuel's E94 "calcsound.f" code at http://texmex.mit.edu/pub/emanuel/BOOK/
"""
function get_lifted_condensation_level(temperature,relative_humidity,pressure) 
    return pressure * (relative_humidity^(temperature/(1669.0-122.0*relative_humidity-temperature)))
end

function get_lifted_condensation_level(temperature :: Quantity ,relative_humidity ,pressure :: Quantity) 
    return pressure * (relative_humidity^(temperature/(1669.0u"K"-122.0u"K"*relative_humidity-temperature)))
end




##### Specific to potential intensity

"""

"""
function ∂specific_entropy_∂temp(temperature, mixing_ratio)
    ∂specific_entropy_∂temp = (Dryair.cp + mixing_ratio * Liquidwater.cp)/temperature - Liquidwater.Lv * mixing_ratio / temperature^2 
end

function ∂specific_entropy_∂temp_emanuel(temperature, mixing_ratio, pressure)
    CL = Liquidwater.cp - 1690.0f0u"J/kg/K" # This is a modified value of liquid water specific heat capacity to compensate for the lack of explicit ice phase when lifting a parcel(Personal communication with Kerry Emanuel on April 22 2022 - Argel Ramirez Reyes
    alv = Liquidwater.Lv + (Watervapor.cp - CL)*(temperature - 273.15u"K")
    saturation_vapor_pressure = get_saturation_vapor_pressure(temperature)
    saturation_mixing_ratio = get_mixing_ratio(saturation_vapor_pressure, pressure)
    ∂specific_entropy_∂temp = (Dryair.cp + mixing_ratio * Liquidwater.cp + alv^2 * saturation_mixing_ratio /(Watervapor.R*temperature^2))/temperature

end


"""
    get_specific_entropy_emanuel(temperature,mixing_ratio,pressure)
Receive temperature in Kelvin, water vapor mixing ratio (unitless g/g) and pressure (hPa) and compute the specific entropy of a parcel using equation in Emmanuel's (E94, EQN. 4.5.9)
"""
function get_specific_entropy_emanuel(temperature,mixing_ratio,pressure)
    CL = Liquidwater.cp - 1690.0f0u"J/kg/K" # This is a modified value of liquid water specific heat capacity to compensate for the lack of explicit ice phase (Personal communication with Kerry Emanuel
    alv = Liquidwater.Lv + (Watervapor.cp - CL)*(temperature - 273.15u"K")
    saturation_vapor_pressure = get_saturation_vapor_pressure(temperature)
    #vapor_pressure = get_partial_vapor_pressure(saturation_vapor_pressure,pressure)
    saturation_mixing_ratio = get_mixing_ratio(saturation_vapor_pressure, pressure)
    specific_entropy =  (Dryair.cp + mixing_ratio * CL) *
        log(temperature/unit(temperature)) - Dryair.R * log((pressure - saturation_vapor_pressure)/unit(pressure)) +
        alv * saturation_mixing_ratio / temperature 
end