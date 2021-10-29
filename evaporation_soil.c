/// @file evaporation_soil.c
/// @brief This module will calculate evaporation from ground surface/top soil,
///        and the evaporation of snow and pond water on surface.
/// @author Edited by XZ Luo
/// @date May 25, 2015


#include "beps.h"


/// @brief Function to calculate evaporation from ground surface/top soil,
///        and the evaporation of snow and pond water on surface.
/// @details [input] air temperature; ground surface temperature;
///                  relative humidity of ground (BEPS takes it as the air RH);
///                  percentage of snow cover on ground; depth of water; depth of snow
///                  soil water content on first soil layer; porosity of first soil layer
/// @details [output] evaporation from soil surface;
///                   depth of water and snow on ground after evaporation and sublimation
/// @param temp_air         air temperature
/// @param temp_g           ground temperature
/// @param rh_air           relative humidity of air
/// @param netRad_g         net radiation on ground
/// @param Gheat_g          aerodynamic conductance of heat on ground surface
/// @param percent_snow_g   percentage of snow on ground
/// @param depth_water      depth of water on ground, after rainfall and snowfall stage 1, before evaporation.
///                         output after subtracting evaporation
/// @param depth_snow       depth of snow on ground, ...
/// @param mass_water_g     mass of water on ground, output after subtracting evaporation
/// @param mass_snow_g      mass of snow on ground, ...
/// @param density_snow     density of snow, from snowpack stage1
/// @param swc_g            soil water content (from last step)
/// @param porosity_g       porosity on ground
/// @param evapo_soil       evaporation from soil
/// @param evapo_water_g    evaporation from pond water
/// @param evapo_snow_g     evaporation from snow on surface
/// @return void
void evaporation_soil(double temp_air, double temp_g, double rh_air, double netRad_g, double Gheat_g,
                      double* percent_snow_g,double* depth_water, double* depth_snow, double* mass_water_g, double* mass_snow_g,
                      double density_snow, double swc_g, double porosity_g,
                      double* evapo_soil, double* evapo_water_g, double* evapo_snow_g)
{
    double meteo_pack_output[10];
    double density_air_g;  // density of air
    double cp_air_g;  // specific heat of air
    double vpd_g;  // vpd of air
    double slope_vapor_g;  // slope of vapor
    double psy_air_g;  // psychrometer on ground
    double Gwater_g; // conductance of water on soil surface
    double latent_water, latent_snow;
    double density_water;
    double length_step;
    //double change_depth_water,change_depth_snow; // change of water depth in this step

    /********************************/
    meteo_pack (temp_g, rh_air, meteo_pack_output);
    density_air_g = meteo_pack_output[1];
    cp_air_g = meteo_pack_output[2];
    vpd_g = meteo_pack_output[3];
    slope_vapor_g = meteo_pack_output[4];
    psy_air_g = meteo_pack_output[5];

    latent_water=(2.501-0.00237*temp_air)*1000000;
    latent_snow=2.83*1000000;

    density_water=1025.0;
    /********************************/

    Gwater_g=1/(4.0*exp(8.2-4.2*swc_g/porosity_g));
    length_step=kstep;

    // get the percentage of snow
    if (*depth_snow>0.02) *percent_snow_g=1;
    else *percent_snow_g=(*mass_snow_g)/(0.025*density_snow);
    *percent_snow_g=max(*percent_snow_g,0);
    *percent_snow_g=min(*percent_snow_g,1);

    // when there are pond water on ground, there is evaporation from the water
    if (*depth_water>0 && *depth_snow == 0)
    {
        *evapo_water_g=1/(latent_water )*(slope_vapor_g *(netRad_g*0.8-0)+density_air_g *cp_air_g *vpd_g*Gheat_g)/(slope_vapor_g+psy_air_g*(1+(Gheat_g)/0.01));
    }
    else
    {
        *evapo_water_g = 0;
    }

    *evapo_water_g = max(-0.002/length_step,*evapo_water_g);
    //*evapo_water_g = max(0,*evapo_water_g);
    if(*evapo_water_g >0) *evapo_water_g =min(*evapo_water_g,(*depth_water)*density_water/length_step);

    *depth_water=*depth_water-(*evapo_water_g)/density_water*length_step;
    *depth_water=max(0,*depth_water);
    *mass_water_g=*mass_water_g-(*evapo_water_g)*length_step;

    // when there are snow on ground, there is only evaporation from the snow
    if (*depth_snow >0)
    {
        *evapo_snow_g=1/(latent_snow )*(slope_vapor_g *(netRad_g*0.8-0)+density_air_g *cp_air_g *vpd_g*Gheat_g)/(slope_vapor_g+psy_air_g*(1+(Gheat_g)/0.01))*(*percent_snow_g);
    }
    else
    {
        *evapo_snow_g = 0;
    }

    *evapo_snow_g = max(-0.002/length_step,*evapo_snow_g);
    //*evapo_snow_g = max(0,*evapo_snow_g);
    if(*evapo_snow_g >0) *evapo_snow_g =min(*evapo_snow_g,*mass_snow_g/length_step);

    *mass_snow_g=*mass_snow_g-(*evapo_snow_g)*length_step;
    *mass_snow_g=max(*mass_snow_g,0);

    if(*mass_snow_g>0)
        *depth_snow=*depth_snow-(*evapo_snow_g)/density_snow*length_step;
    else
        *depth_snow=0;


    // if evaporation from pond and snow have not consumed all energy, then no evaporation from soil surface
	// change_depth_water=(*evapo_water_g/density_water);
    // change_depth_snow=(*evapo_snow_g/density_snow);

//
    if ( *depth_water >0 || *depth_snow >0 )
        *evapo_soil=0;
    else
    {
        *evapo_soil=(1-(*percent_snow_g))*1/(latent_water )*(slope_vapor_g*(netRad_g-0)+density_air_g*cp_air_g*vpd_g*Gheat_g)/(slope_vapor_g +psy_air_g*(1 + Gheat_g/Gwater_g));
        *evapo_soil = max(0,*evapo_soil);
    }

}
