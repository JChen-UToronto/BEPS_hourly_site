/// @file meteo_pack.c
/// @brief This function will calculate all the meteorological variables based on input
/// @author Edited by XZ Luo
/// @date May 19, 2015


# include "beps.h"


/// @brief Function to calculate meteorological variables based on input
/// @details default input is temperature (C) and relative humidity (0-100)
///          output is an array, named as meteo_pack_output []
/// @details [input] meteo_pack_output [1]= air_density kg/m3 \n
///                  meteo_pack_output [2]= specific heat of air J/kg/C \n
///                  meteo_pack_output [3]= VPD kPa \n
///                  meteo_pack_output [4]= slope of vapor pressure to temperature kPa/C \n
///                  meteo_pack_output [5]= psychrometer constant kPa/C \n
///                  meteo_pack_output [6]= saturate water vapor potential kPa \n
///                  meteo_pack_output [7]= actual water vapor potential kPa \n
///                  meteo_pack_output [8]= specific humidity g/g
/// @param temp               temperature
/// @param rh                 relative humidity
/// @param meteo_pack_output  meteorological variables array
/// @return void
void meteo_pack(double temp, double rh, double* meteo_pack_output)
{
    double density_air, cp_air, vpd, slope_vapor, psy, e_saturate, e_actual, sp_humidity;
    double latent_water;  //latent heat of water J/kg
    density_air=1.292;
    e_saturate=0.61078*exp(17.3*temp/(237.3+temp));
    e_actual=e_saturate*rh/100;
    vpd=e_saturate-e_actual;
    sp_humidity=0.622*e_actual/(101.35-0.378*e_actual);
    cp_air=1004.65*(1+0.84*sp_humidity);
    slope_vapor=2503.0/pow((temp+237.3),2)*exp(17.27*temp/(temp+237.3));
    psy=0.066;
    latent_water=(2.501-0.00237*temp)*1000000;
    //psy=cp_air*101.13/(0.622*latent_water);

    meteo_pack_output[1]= density_air;
    meteo_pack_output[2]= cp_air;
    meteo_pack_output[3]= vpd;
    meteo_pack_output[4]= slope_vapor;
    meteo_pack_output[5]= psy;
    meteo_pack_output[6]= e_saturate;
    meteo_pack_output[7]= e_actual;
    meteo_pack_output[8]= sp_humidity;
}
