/// @file transpiration.c
/// @brief This module calculates transpiration, for overstorey and understorey, sunlit and shaded
/// @author Edited by XZ Luo
/// @date May 20, 2015


#include "beps.h"


/// @brief Function to calculate transpiration
/// @details A transformation of Penman-Monteith equation is used here.
///          It could be regarded as a mass transfer process.
///          Water vapor inside cells are required by VPD from air and VPD on leaf surface.
/// @details [input] temperature of sunlit and shaded leaves from other storey (leaf temperature module);
///                  temperature of air; relative humidity;
///                  conductance of water for sunlit shaded leaves from overstorey and understorey;
///                  leaf area index, sunlit and shaded, overstorey and understorey (from leaf area index module);
/// @details [output] transpiration from overstorey and understorey
/// @param tempL_o_sunlit  temperature of leaf, overstory, sunlit
/// @param tempL_o_shaded  temperature of leaf, overstory, shaded
/// @param tempL_u_sunlit  temperature of leaf, understory, sunlit
/// @param tempL_u_shaded  temperature of leaf, understory, shaded
/// @param temp_air        air temperature
/// @param rh_air          relative humidity of air
/// @param Gtrans_o_sunlit total conductance of water
///                        tandem of stomatal conductance and aerodynamic conductance, overstory, sunlit
/// @param Gtrans_o_shaded ..., overstory, shaded
/// @param Gtrans_u_sunlit ..., understory, sunlit
/// @param Gtrans_u_shaded ..., understory, shaded
/// @param lai_o_sunlit    leaf area index, overstory, sunlit
/// @param lai_o_shaded    leaf area index, overstory, shaded
/// @param lai_u_sunlit    leaf area index, understory, sunlit
/// @param lai_u_shaded    leaf area index, understory, shaded
/// @param trans_o         transpiration from overstory
/// @param trans_u         transpiration from understory
/// @return void
void transpiration (double tempL_o_sunlit, double tempL_o_shaded, double tempL_u_sunlit, double tempL_u_shaded,
                    double temp_air, double rh_air,
                    double Gtrans_o_sunlit, double Gtrans_o_shaded, double Gtrans_u_sunlit, double Gtrans_u_shaded,
                    double lai_o_sunlit, double lai_o_shaded, double lai_u_sunlit, double lai_u_shaded,
                    double* trans_o, double* trans_u)
{
    double LHt_o_sunlit, LHt_o_shaded, LHt_u_sunlit, LHt_u_shaded;  // latent heat from leaves W/m2
    double meteo_pack_output[10];
    double density_air, cp_air, vpd_air, slope_vapor_air, psy_air;
    double latent_water;

    meteo_pack (temp_air, rh_air, meteo_pack_output);
    density_air = meteo_pack_output [1];
    cp_air = meteo_pack_output [2]; // specific heat of moist air above canopy
    vpd_air = meteo_pack_output [3];
    slope_vapor_air = meteo_pack_output [4]; // slope of saturated vapor potential to temperature
    psy_air = meteo_pack_output [5];  // psychrometer constant

    latent_water=(2.501-0.00237*temp_air)*1000000;


    LHt_o_sunlit =(vpd_air+slope_vapor_air *(tempL_o_sunlit -temp_air ))*density_air*cp_air*Gtrans_o_sunlit /psy_air;
    LHt_o_shaded =(vpd_air+slope_vapor_air *(tempL_o_shaded -temp_air ))*density_air*cp_air*Gtrans_o_shaded /psy_air;
    LHt_u_sunlit =(vpd_air+slope_vapor_air *(tempL_u_sunlit -temp_air ))*density_air*cp_air*Gtrans_u_sunlit /psy_air;
    LHt_u_shaded =(vpd_air+slope_vapor_air *(tempL_u_shaded -temp_air ))*density_air*cp_air*Gtrans_u_shaded /psy_air;


    *trans_o =1/(latent_water )*(LHt_o_sunlit *lai_o_sunlit +LHt_o_shaded *lai_o_shaded );
    *trans_u =1/(latent_water )*(LHt_u_sunlit *lai_u_sunlit +LHt_u_shaded *lai_u_shaded );
}
