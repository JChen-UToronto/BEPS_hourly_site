/// @file surface_temp.c
/// @brief This module will simulate surface temperature in each step, as well as heat flux form surface to soil layers.
/// @details As it is an interface between ground, air and soil,
///          the core idea is to separate the interface as different layers by depth of snow,
///          then calculate the temperature gradient and at last calculate the heat flux from ground surface to soil.
/// @details Original beps would use Xg_snow[kkk] at some places after snow melt & frozen,
///          now we uniformly use the value before snow melt & frozen.
/// @author Edited by XZ Luo
/// @date June 1, 2015


#include "beps.h"


/// @brief Function to simulate surface temperature, and heat flux from surface to soil layers
/// @param temp_air             air temperature (Celsius degree)
/// @param rh_air               relative humidity (0-100)
/// @param depth_snow           depth of snow (m)
/// @param depth_water          depth of water on ground (m)
/// @param capacity_heat_soil1  heat capacity of layer1 soil (J/m2/K)
/// @param capacity_heat_soil0  heat capacity of layer2 soil (J/m2/K)
/// @param Gheat_g              aerodynamic conductance of heat on ground (m/s)
/// @param depth_soil1          depth of soil in layer1 (m)
/// @param density_snow         density of snow (kg/m3)
/// @param tempL_u              leaf temperature, understory (Celsius degree)
/// @param netRad_g             net radiation on ground (W/m2)
/// @param evapo_soil           evaporation from soil surface (mm/s)
/// @param evapo_water_g        evaporation from pond water on ground (mm/s)
/// @param evapo_snow_g         evaporation from snow pack on ground (mm/s)
/// @param lambda_soil1         thermal conductivity of layer1 soil (W/m/K)
/// @param percent_snow_g       percentage of snow coverage on ground (0-1)
/// @param heat_flux_soil1      heat flux from layer1 soil to the next soil layer (W/m2)
/// @param temp_ground_last     temperature of ground, from last step
/// @param temp_soil1_last      temperature of layer1 soil, from last step
/// @param temp_any0_last       temperature of any layer right above the soil, from last step
/// @param temp_snow_last       temperature of snow, from last step
/// @param temp_soil0_last      temperature of soil0, from last step
/// @param temp_snow1_last      temperature of snow layer 2, from last step
/// @param temp_snow2_last      temperature of snow layer 3, from last step
/// @param temp_ground          ground temperature at this step
/// @param temp_any0            temperature of any layer right above the soil
///                             could be a mixture of snow temperature and soil surface temperature
/// @param temp_snow            snow temperature at this step
/// @param temp_soil0           temperature of soil surface right above the soil, the part not covered by snow
/// @param temp_snow1           temperature of snow layer 2, used in depth_snow>0.05 m
/// @param temp_snow2           temperature of snow layer 3, used in depth_snow>0.05 m
/// @param heat_flux            heat flux from ground to soil
/// @return void
void surface_temperature(double temp_air,double rh_air, double depth_snow, double depth_water,
                         double capacity_heat_soil1, double capacity_heat_soil0, double Gheat_g,
                         double depth_soil1, double density_snow,double tempL_u, double netRad_g,
                         double evapo_soil, double evapo_water_g, double evapo_snow_g, double lambda_soil1,
                         double percent_snow_g, double heat_flux_soil1, double temp_ground_last,
                         double temp_soil1_last, double temp_any0_last, double temp_snow_last,
                         double temp_soil0_last, double temp_snow1_last, double temp_snow2_last,
                         double* temp_ground, double* temp_any0, double* temp_snow,
                         double* temp_soil0, double* temp_snow1, double* temp_snow2, double* heat_flux)
{
    double length_step;
    double meteo_pack_output[10];
    double density_air, cp_air;   // density of air, specific heat of air
    double cp_ice;                // specific heat of ice
    double latent_water, latent_snow;
    double Gg;                    // radiation available for heating the ground;
    double lambda_snow;           // thermal conductivity of snow in this step
    double heat_flux_soil, heat_flux_snow;  // heat flux through the soil and snow fraction on ground, separately.
    double heat_flux_snow1, heat_flux_snow2;

    double ra_g; // aerodynamic resistance of heat

    double ttt;  // temporary variables

    length_step=kstep;

    cp_ice=2228.261;
    latent_water=(2.501-0.00237*temp_air)*1000000;
    latent_snow=2.83*1000000;
    meteo_pack (temp_air, rh_air, meteo_pack_output);
    density_air = meteo_pack_output[1];
    cp_air = meteo_pack_output[2];
    ra_g=1/Gheat_g;

    // thermal conductivity of snow
    lambda_snow=0.021+4.2*density_snow/10000+2.2*pow(density_snow,3)*pow(10,-9);

    // available energy on ground for
    Gg=netRad_g - evapo_snow_g*latent_snow-(evapo_water_g+evapo_soil)*latent_water;

    // case 1, snow depth<2cm, snow temperature, ground surface temperature, soil surface temperature are the same
    if(depth_snow<=0.02)
    {
        ttt=capacity_heat_soil1*0.02/length_step;
        *temp_ground=(temp_ground_last*(ttt)*ra_g*depth_soil1+Gg*ra_g*depth_soil1+density_air*cp_air*temp_air*depth_soil1+ ra_g*lambda_soil1*temp_soil1_last);
        *temp_ground=*temp_ground/(density_air*cp_air*depth_soil1+ ra_g*lambda_soil1+ttt*ra_g*depth_soil1);

        *temp_ground=max(temp_ground_last-25,*temp_ground);
        *temp_ground=min(temp_ground_last+25,*temp_ground);

        *temp_any0 = *temp_ground;
        *temp_snow  =*temp_any0;
        *temp_soil0 =*temp_any0;
        *temp_snow1 =*temp_any0;
        *temp_snow2 =*temp_any0;

        *heat_flux=2*lambda_soil1*((*temp_any0)-temp_soil1_last)/depth_soil1;
        *heat_flux=min(100,*heat_flux);
        *heat_flux=max(-100,*heat_flux);
    }

    // case 2: depth of snow larger than 2cm, smaller than 5 cm.
    // snow fraction on ground decide the snow temperature based on energy balance
    // soil fraction on ground decide the soil surface temperature based on energy balance
    // snow and soil fraction works in parallel to determine the ground surface temperature
    else if(depth_snow>0.02  && depth_snow<=0.05)
    {
        ttt=capacity_heat_soil1*0.02/length_step; // for soil fraction part

        *temp_soil0=(temp_soil0_last*ttt*ra_g*depth_soil1+Gg*ra_g*depth_soil1+density_air*cp_air*temp_air*depth_soil1+ 2*ra_g*lambda_soil1*temp_soil1_last)/
                    (density_air*cp_air*depth_soil1+ 2*ra_g*lambda_soil1+ttt*ra_g*depth_soil1);

        *temp_soil0=max(temp_air-25,(*temp_soil0));
        *temp_soil0=min(temp_air+25,(*temp_soil0));

        ttt=cp_ice*density_snow*depth_snow/length_step;  /*for snow part*/
        *temp_snow=(temp_snow_last*ttt*ra_g*depth_snow+Gg*ra_g*depth_snow+density_air*cp_air*tempL_u*depth_snow+ ra_g*lambda_snow*temp_any0_last)/
                   (density_air*cp_air*depth_snow+ ra_g*lambda_snow+ttt*ra_g*depth_snow);

        *temp_snow=max(temp_air-25,*temp_snow);
        *temp_snow=min(temp_air+25,*temp_snow);

        ttt=(lambda_soil1*temp_soil1_last/depth_soil1+(*temp_snow)*lambda_snow+0.02*capacity_heat_soil1/length_step*temp_any0_last)/
            (lambda_soil1/depth_soil1+lambda_snow/depth_snow+0.02*capacity_heat_soil1/length_step);
        *temp_any0=(*temp_soil0)*(1-percent_snow_g)+ttt*percent_snow_g;

        heat_flux_snow = lambda_snow/(depth_snow+0.5*depth_soil1)*((*temp_snow)-temp_soil1_last);
        heat_flux_soil = heat_flux_snow*((*temp_any0)-temp_soil1_last)/depth_soil1;

        *heat_flux=heat_flux_snow*percent_snow_g+heat_flux_soil*(1-percent_snow_g);

        *heat_flux=min(100,*heat_flux);
        *heat_flux=max(-100,*heat_flux);

        // starting to melt
        if(*temp_snow>zero && temp_snow_last<=zero && depth_snow>zero)
        {
            *temp_snow=0;
        }

        // starting to frozen
        if(*temp_snow<zero && temp_snow_last>=zero && depth_water>zero)
        {
            *temp_snow=0;
        }

        // percent_snow_g =min(1.0,Wg_snow[kkk] / (0.05 * rho_snow[kkk])); // use the fraction before
        *temp_ground =(*temp_snow)*percent_snow_g+(*temp_soil0)*(1-percent_snow_g);
        *temp_ground=max(temp_air-25,*temp_ground);
        *temp_ground=min(temp_air+25,*temp_ground);

        *temp_snow1=*temp_snow;
        *temp_snow2=*temp_snow;
    }

        // case 3, snow depth is larger than 5 cm
        // snow coverage on ground is 100%
        // the first layer of snow is set as 2 cm
        // the second layer of snow is set as 2 cm, too
        // the depth of third snow layer is depth_snow-0.04;
    else if(depth_snow>0.05)
    {
        ttt=cp_ice*density_snow*0.02/length_step;

        *temp_snow=(temp_snow_last*ttt*ra_g*0.04+Gg*ra_g*0.02+density_air*cp_air*temp_air*0.04+ra_g*lambda_snow*temp_snow1_last)/
                   (density_air*cp_air*0.04+ ra_g*lambda_snow+ttt*ra_g*0.04);
        *temp_snow=max(temp_air-25,*temp_snow);
        *temp_snow=min(temp_air+25,*temp_snow);

        heat_flux_snow=lambda_snow*(*temp_snow-temp_snow1_last)/0.04; // why 0.04 here?

        *heat_flux=heat_flux_snow;
        *heat_flux=min(100,*heat_flux);
        *heat_flux=max(-100,*heat_flux);

        heat_flux_snow1 = lambda_snow*(temp_snow1_last-temp_snow2_last)/(depth_snow-0.02);
        *temp_snow1=temp_snow1_last+((*heat_flux)- heat_flux_snow1)/(cp_ice*density_snow*0.02 )*length_step ;
        heat_flux_snow2 = (temp_snow2_last-temp_any0_last)/(0.5*(depth_snow-0.04)/lambda_snow+0.02/lambda_soil1); //??
        *temp_snow2=temp_snow2_last+(heat_flux_snow1- heat_flux_snow2)/(cp_ice*density_snow*(depth_snow-0.04))*length_step ;

        *temp_any0 = temp_any0_last+(heat_flux_snow2- heat_flux_soil1) / (capacity_heat_soil0 * 0.02)*length_step;
        *temp_soil0=*temp_any0;
        /* starting to melt*/
        if(*temp_snow>zero && temp_snow_last<=zero && depth_snow>zero)
        {
            *temp_snow=0;
        }

        /* starting to frozen */
        if(*temp_snow<zero && temp_snow_last>=zero && depth_water>zero)
        {
            *temp_snow=0;
        }

        *temp_ground=*temp_snow;
    }
}

/*==========================*
theoretical top layer: temp_any0
heat_flux_snow2
snow layer 3: temp_snow2
heat_flux_snow1
snow layer 2: temp_snow1
heat_flux_snow
snow layer 1: temp_snow
*/
