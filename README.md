# TropicalCyclonePotentialIntensity

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://www.argelramirezreyes.com/TropicalCyclonePotentialIntensity.jl/dev/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://www.argelramirezreyes.com/TropicalCyclonePotentialIntensity.jl/dev/)
[![Coverage](https://codecov.io/gh/aramirezreyes/TropicalCyclonePotentialIntensity.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aramirezreyes/TropicalCyclonePotentialIntensity.jl)

# TropicalCyclonePotentialIntensity.jl
This package implements the functions of [Daniel Gliford](https://github.com/dgilford)'s [TCPyPI](https://github.com/dgilford/tcpyPI/blob/master/tcpyPI/pi.py) for the [JuliaLanguage](https://github.com/JuliaLang/julia). This are widely used routines in the Tropical Cyclone community. I hope this helps widen the atmospheric science ecosystem for Julia.

# Installation

This package is very lightweight and it is registered in the Julia General Registry so to install it in the Julia REPL just  hit the ```]``` key to enter Pkg mode, then type
```julia
Pkg> add TropicalCyclonePotentialIntensity
```

# Usage

TropicalCyclonePotentialIntensity.jl has three main functions:

In here, tparcel, pparcel and rparcel are the temperature (in Kelvin), pressure (in hPa) and mixing ratio (kg/kg) of a parcel and tabs, r and pres are environmental profiles (should be equally sized arrays).

```julia
get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,tabs,r,pres)

get_potential_intensity_of_tropical_cyclone(tparcel, pparcel, pres, tabs, r)

get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,tabs,r,pres)
```

# Units
This package uses [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) to manage units. We recommend to use TropicalCyclonePotentialIntensity.jl with unitful quantities to help maintain your results consistent, although it should run without real numbers as well.

## Package under development

This package is under development by Argel Ramírez Reyes. It is definitely not perfect. Bugs are expected as this has not been widely tested. Bug reports or collaborations are highly appreciated.


### Main reference

Gilford, D. M., 2021: pyPI (v1.3): Tropical Cyclone Potential Intensity Calculations in Python. Geoscientific Model Development, 14, 2351–2369, https://doi.org/10.5194/gmd-14-2351-2021.