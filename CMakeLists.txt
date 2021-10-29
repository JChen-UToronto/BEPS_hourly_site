cmake_minimum_required(VERSION 3.17)
project(BEPS_Hourly_site C)

set(CMAKE_C_STANDARD 99)

add_executable(BEPS_Hourly_DS bepsmain_pnt.c
        beps.h DB.h debug.h soil.h
        aerodynamic_conductance.c calc_temp_leaf.c evaporation_canopy.c evaporation_soil.c
        init_soil.c inter_prg.c lai_calc.c meteo_pack.c netRadiation.c photosyn_gs.c
        plant_respir.c rainfall.c readcoef.c readparam.c s_coszs.c sensible_heat.c snowpack.c
        soil_thermal_regime.c soil_water_stress.c soilresp.c surface_temp.c transpiration.c
        updatesoilmoisture.c vcmax_jmax_mod.c
        )