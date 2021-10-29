/// @file evaporation_canopy.c
/// @brief This module calculates evaporation and sublimation from canopy, from overstorey understorey sunlit and shaded
/// @author Edited by XZ Luo
/// @date May 25, 2015


#include "beps.h"


/// @brief Function to calculate evaporation and sublimation from canopy
/// @details [input] temperature of sunlit and shaded leaves from other storey (leaf temperature module);
///                  temperature of air; relative humidity;
///                  aerodynamic conductance of water (snow) for sunlit shaded leaves from overstorey and understorey;
///                  percentage of overstorey or understorey covered by water or snow;
///                  leaf area index, sunlit and shaded, overstorey and understorey (from leaf area index module);
/// @details [output] evaporation of water and snow from overstorey and understorey
/// @param tempL_o_sunlit    temperature of leaves, overstory, sunlit (leaf temperature module)
/// @param tempL_o_shaded    temperature of leaves, overstory, shaded
/// @param tempL_u_sunlit    temperature of leaves, understory, sunlit
/// @param tempL_u_shaded    temperature of leaves, understory, shaded
/// @param temp_air          air temperature
/// @param rh_air            relative humidity
/// @param Gwater_o_sunlit   aerodynamic conductance of water (snow) for overstory, sunlit leaves
/// @param Gwater_o_shaded   aerodynamic conductance of water (snow) for overstory, shaded leaves
/// @param Gwater_u_sunlit   aerodynamic conductance of water (snow) for understory, sunlit leaves
/// @param Gwater_u_shaded   aerodynamic conductance of water (snow) for understory, shaded leaves
/// @param lai_o_sunlit      leaf area index, overstory, sunlit (from leaf area index module)
/// @param lai_o_shaded      leaf area index, overstory, shaded
/// @param lai_u_sunlit      leaf area index, understory, sunlit
/// @param lai_u_shaded      leaf area index, understory, shaded
/// @param percent_water_o   percentage of overstorey covered by water
/// @param percent_water_u   percentage of understorey covered by water
/// @param percent_snow_o    percentage of overstorey covered by snow
/// @param percent_snow_u    percentage of understorey covered by snow
/// @param evapo_water_o     evaporation of water from overstorey
/// @param evapo_water_u     evaporation of water from understorey
/// @param evapo_snow_o      evaporation of snow from overstorey
/// @param evapo_snow_u      evaporation of snow from understorey
/// @return void
void evaporation_canopy(double tempL_o_sunlit, double tempL_o_shaded, double tempL_u_sunlit, double tempL_u_shaded,
                        double temp_air, double rh_air,
                        double Gwater_o_sunlit, double Gwater_o_shaded, double Gwater_u_sunlit, double Gwater_u_shaded,
                        double lai_o_sunlit, double lai_o_shaded, double lai_u_sunlit, double lai_u_shaded,
                        double percent_water_o, double percent_water_u, double percent_snow_o, double percent_snow_u,
                        double* evapo_water_o, double* evapo_water_u, double* evapo_snow_o, double* evapo_snow_u)
{
    double LHw_o_sunlit, LHw_o_shaded, LHw_u_sunlit, LHw_u_shaded; // latent heat from leaves W/m2, caused by evaporation of intercepted rain
    double LHs_o_sunlit, LHs_o_shaded, LHs_u_sunlit, LHs_u_shaded; // latent heat from leaves W/m2, caused by evaporation of intercepted snow
    double meteo_pack_output[10];
    double density_air, cp_air, vpd_air, slope_vapor_air, psy_air;
    double latent_water, latent_snow;

    meteo_pack (temp_air, rh_air, meteo_pack_output);
    density_air = meteo_pack_output [1];
    cp_air = meteo_pack_output [2];
    vpd_air = meteo_pack_output [3];
    slope_vapor_air = meteo_pack_output [4];
    psy_air = meteo_pack_output [5];

    latent_water=(2.501-0.00237*temp_air)*1000000;
    latent_snow=2.83*1000000;

    // leaf level latent heat caused by evaporation or sublimation
    LHw_o_sunlit =percent_water_o*(vpd_air+slope_vapor_air *(tempL_o_sunlit -temp_air ))*density_air*cp_air*Gwater_o_sunlit /psy_air;
    LHw_o_shaded =percent_water_o*(vpd_air+slope_vapor_air *(tempL_o_shaded -temp_air ))*density_air*cp_air*Gwater_o_shaded /psy_air;

    LHw_u_sunlit =percent_water_u*(vpd_air+slope_vapor_air *(tempL_u_sunlit -temp_air ))*density_air*cp_air*Gwater_u_sunlit /psy_air;
    LHw_u_shaded =percent_water_u*(vpd_air+slope_vapor_air *(tempL_u_shaded -temp_air ))*density_air*cp_air*Gwater_u_shaded /psy_air;

    LHs_o_sunlit =percent_snow_o*(vpd_air+slope_vapor_air *(tempL_o_sunlit -temp_air ))*density_air*cp_air*Gwater_o_sunlit /psy_air;
    LHs_o_shaded =percent_snow_o*(vpd_air+slope_vapor_air *(tempL_o_shaded -temp_air ))*density_air*cp_air*Gwater_o_shaded /psy_air;

    LHs_u_sunlit =percent_snow_u*(vpd_air+slope_vapor_air *(tempL_u_sunlit -temp_air ))*density_air*cp_air*Gwater_u_sunlit /psy_air;
    LHs_u_shaded =percent_snow_u*(vpd_air+slope_vapor_air *(tempL_u_shaded -temp_air ))*density_air*cp_air*Gwater_u_shaded /psy_air;


    /*******************************************/
    *evapo_water_o =1/(latent_water )*(LHw_o_sunlit *lai_o_sunlit +LHw_o_shaded *lai_o_shaded );
    *evapo_water_u =1/(latent_water )*(LHw_u_sunlit *lai_u_sunlit +LHw_u_shaded *lai_u_shaded );

    *evapo_snow_o =1/(latent_snow )*(LHs_o_sunlit *lai_o_sunlit +LHs_o_shaded *lai_o_shaded );
    *evapo_snow_u =1/(latent_snow )*(LHs_u_sunlit *lai_u_sunlit +LHs_u_shaded *lai_u_shaded );

    // to eliminate negative ET value
//    *evapo_water_o=max(0,*evapo_water_o);
//    *evapo_water_u=max(0,*evapo_water_u);
//    *evapo_snow_o=max(0,*evapo_snow_o);
//    *evapo_snow_u=max(0,*evapo_snow_u);
}
