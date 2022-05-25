
const epsilon       = 18.016/28.966
const g             = 10u"m/s/s" #acceleration of gravity

struct Substance{T} 
    cp :: Union{Nothing,T}
    cv :: Union{Nothing,T}
    R  :: Union{Nothing,T}
    Lv :: Union{Nothing,T}
    Lf :: Union{Nothing,T}
end

Substance{T}(;cp = nothing, cv = nothing, R = nothing, Lv = nothing, Lf = nothing) where T = Substance{T}(cp,cv,R,Lv,Lf) 

const Dryair = Substance{Quantity}(
    cp = 1005.7u"J/kg/K", #J/kg/k at 1013 hPa
    cv = 718.0u"J/kg/K",
    R  = 287.05u"J/kg/K" # J/kg/k
)

const Liquidwater = Substance{Quantity}(
     Lv = 2.501e6u"J/kg", #J/kg
     Lf = 3.33e5u"J/kg",
     cp = 4190.0u"J/kg/K" #j/kg/k
)

const Watervapor = Substance{Quantity}(
    R = 461.52u"J/kg/K", #j/kg/K
    cp = 1870.0u"J/kg/K"
)
