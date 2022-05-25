@testset "Thermodynamics" begin
    @test mixing_ratio_to_specific_humidity(specific_humidity_to_mixing_ratio(0.5)) ≈ 0.5
    @test get_saturation_vapor_pressure(273.15u"K") == 6.112u"hPa"
    @test get_saturation_vapor_pressure(273.15) == 6.112
    @test get_partial_vapor_pressure(0,1000u"hPa") == 0u"hPa"
    @test get_partial_vapor_pressure(1,1000u"hPa") == 1000u"hPa"/(18.016/28.966 + 1.0)
    @test get_mixing_ratio(0u"hPa",1000u"hPa") == 0
    @test get_mixing_ratio(get_partial_vapor_pressure(0.5,1000.0), 1000.0) == 0.5
    @test get_partial_vapor_pressure(0,1000) == 0
    @test get_partial_vapor_pressure(1,1000) == 1000/(18.016/28.966 + 1.0)
    @test get_mixing_ratio(0,1000) == 0
    @test get_mixing_ratio(get_partial_vapor_pressure(0.5,1000.0), 1000.0) == 0.5
    @test unit(get_specific_entropy(300u"K",0.2,1000u"hPa"))== u"J/K/kg"
    @test get_potential_temperature(300,1000,1000) == 300
    @test get_potential_temperature(300,1010,1000) < 300
    @test get_potential_temperature(300,900,1000) > 300
    @test get_virtual_temperature(300,0) == 300
    @test get_virtual_temperature(300,0) == 300
    @test get_virtual_temperature(300,10) > 300
    @test get_potential_temperature(300u"K",1000u"hPa",1000u"hPa") == 300u"K"
    @test get_potential_temperature(300u"K",1010u"hPa",1000u"hPa") < 300u"K"
    @test get_potential_temperature(300u"K",900u"hPa",1000u"hPa") > 300u"K"
    @test get_virtual_temperature(300u"K",0u"g/kg") == 300u"K"
    @test get_virtual_temperature(300u"K",0u"g/kg") == 300u"K"
    @test get_virtual_temperature(300u"K",10u"g/kg") > 300u"K"

end
