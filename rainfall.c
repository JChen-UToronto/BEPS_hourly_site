/// @file rainfall.c
/// @brief This module will calculate the water remained on canopy surface after evaporation in this step
///        (used for next step)
/// @details [rainfall_stage1] happens before evaporation of intercepted water from canopy (supply)
/// @details [rainfall_stage2] happens after evaporation of intercepted water from canopy (demand)
/// @note rainfall on ground is considered in stage1,
///       and then considered in surface water module (or soil moisture module)
/// @author XZ Luo
/// @date May 25, 2015


#include "beps.h"


/// @brief Function of rainfall stage1.
/// @details [rainfall_stage1] happens before evaporation of intercepted water from canopy (supply)
/// @details [input] air temperature, precipitation (m/s),
///                  remain of water on leaves from last step (kg/m2) per leaf area
///                  leaf area index of overstorey and understorey, excluding stem.
///                  length of this step (s), if time step is 10min, then it is set as 600,
///                  air temperature and humidity
/// @details [output] percentage of canopy covered by rainfall,
///                   overstorey and understorey (provided to evaporation_canopy),
///                   mass of water available for evaporation on canopy in this step
///                   precipitation on ground
/// @details [optical output] intercepted mass of rainfall in this step
/// @param  temp_air           air temperature (Celsius)
/// @param  precipitation      precipitation rate (m/s)
/// @param  mass_water_o_last  remains of water from last step, overstory
/// @param  mass_water_u_last  remains of water from last step, understory
/// @param  lai_o              leaf area index, overstory
/// @param  lai_u              leaf area index, understory
/// @param  clumping           clumping index
/// @param  mass_water_o       mass of water on leaves (kg/m2) per ground area, overstory
/// @param  mass_water_u       mass of water on leaves (kg/m2) per ground area, understory
/// @param  percent_water_o    the fraction of canopy covered by liquid water and snow, overstory
/// @param  percent_water_u    the fraction of canopy covered by liquid water and snow, understory
/// @param  precipitation_g    precipitation on ground
/// @return void
void rainfall_stage1(double temp_air, double precipitation, double mass_water_o_last, double mass_water_u_last,
                     double  lai_o, double lai_u, double clumping,
                     double* mass_water_o, double* mass_water_u,
                     double* percent_water_o, double* percent_water_u, double* precipitation_g)
{
    double precipitation_o, precipitation_u;   // rate of precipitation on overstorey, understorey and ground

    double massMax_water_o, massMax_water_u;   // Maximum mass of water could be intercepted per leaf area, in kg/m2
    double massStep_water_o, massStep_water_u; // mass of water intercepted in this step per leaf area, in kg/m2
    double density_water;

    double length_step;

    length_step=kstep;

    if (temp_air < 0|| temp_air == 0) // temperature > 0, otherwise it is snow fall
        precipitation=0;

    density_water=1025.0;
    *precipitation_g=0;

    // overstorey
    precipitation_o=precipitation;
    *mass_water_o=mass_water_o_last+precipitation_o*length_step*density_water*(1-exp(-lai_o*clumping));
    massMax_water_o=0.1*lai_o;

    *mass_water_o=max(0,*mass_water_o);
    *mass_water_o=min(massMax_water_o, *mass_water_o);

    massStep_water_o=(*mass_water_o)-mass_water_o_last;
    massStep_water_o=max(0,massStep_water_o);

    *percent_water_o = (*mass_water_o)/massMax_water_o;
    *percent_water_o = min(1,*percent_water_o);


    // understorey
    precipitation_u=precipitation_o-(massStep_water_o)/density_water/length_step;
    *mass_water_u=mass_water_u_last+precipitation_u*length_step*density_water*(1-exp(-lai_u*clumping));
    massMax_water_u=0.1*lai_u;

    *mass_water_u=max(0,*mass_water_u);
    *mass_water_u=min(massMax_water_u, *mass_water_u);

    massStep_water_u=(*mass_water_u)-mass_water_u_last;
    massStep_water_u=max(0,massStep_water_u);

    *percent_water_u = (*mass_water_u)/massMax_water_u;
    *percent_water_u = min(1,*percent_water_u);


    // ground
    *precipitation_g=precipitation_u-(massStep_water_u)/density_water/length_step;
}


/// @brief Function of rainfall stage2.
/// @details [rainfall_stage2] happens after evaporation of intercepted water from canopy (demand)
/// @details [input] mass of water on leaves after precipitation in this step,
///                  evaporation from leaves in this step
/// @details [output] mass of water on leaves after the evaporation on leaves in this step
///                   (this value is transferred to next step)
/// @param  evapo_water_o evaporation of intercepted rain in this step, overstorey, kg/m2/s = mm/s
/// @param  evapo_water_u evaporation of intercepted rain in this step, understorey, kg/m2/s = mm/s
/// @param  mass_water_o  supply of rain on leaves, overstory, already added precipitation in this step
/// @param  mass_water_u  supply of rain on leaves, understory, already added precipitation in this step
/// @return void
void rainfall_stage2(double evapo_water_o, double evapo_water_u,
                     double* mass_water_o, double* mass_water_u)
{
    double length_step; // length of step
    length_step=kstep;  // 6min or 360s per step

    *mass_water_o=*mass_water_o-evapo_water_o*length_step;
    *mass_water_o=max(0,*mass_water_o);

    *mass_water_u=*mass_water_u-evapo_water_u*length_step;
    *mass_water_u=max(0,*mass_water_u);
}

