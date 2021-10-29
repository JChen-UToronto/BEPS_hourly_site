/// @file sensible_heat.c
/// @brief This module will calculate sensible heat from overstorey, understorey and ground
/// @author Edited by XZ Luo
/// @date May 23, 2015


# include "beps.h"


/// @brief Function to calculate sensible heat
/// @details [input] temperature of sunlit and shaded leaves from other storey (leaf temperature module);
///                  temperature of air; relative humidity;
///                  temperature of ground (soil heat flux module);
///                  aerodynamic heat conductance of sunlit shaded leaves from overstorey and understorey;
///                  aerodynamic heat conductance of ground;
///                  leaf area index, sunlit and shaded, overstorey and understorey (from leaf area index module);
/// @details [output] sensible heat from overstorey, understorey and ground
/// @param tempL_o_sunlit    temperature of leaves, overstory, sunlit
/// @param tempL_o_shaded    temperature of leaves, overstory, shaded
/// @param tempL_u_sunlit    temperature of leaves, understory, sunlit
/// @param tempL_u_shaded    temperature of leaves, understory, shded
/// @param temp_g            temperature of ground
/// @param temp_air          air temperature
/// @param rh_air            relative humidity of air
/// @param Gheat_o_sunlit    aerodynamic resistance of heat, overstory, sunlit
/// @param Gheat_o_shaded    aerodynamic resistance of heat, overstory, shaded
/// @param Gheat_u_sunlit    aerodynamic resistance of heat, understory, sunlit
/// @param Gheat_u_shaded    aerodynamic resistance of heat, understory, shaded
/// @param Gheat_g           aerodynamic resistance of heat, ground
/// @param lai_o_sunlit      leaf area index, overstory, sunlit
/// @param lai_o_shaded      leaf area index, overstory, shaded
/// @param lai_u_sunlit      leaf area index, understory, sunlit
/// @param lai_u_shaded      leaf area index, understory, shaded
/// @param SH_o              sensible heat, overstory
/// @param SH_u              sensible heat, understory
/// @param SH_g              sensible heat, ground
/// @return void
void sensible_heat(double tempL_o_sunlit, double tempL_o_shaded, double tempL_u_sunlit, double tempL_u_shaded,
                   double temp_g, double temp_air, double rh_air,
                   double Gheat_o_sunlit, double Gheat_o_shaded, double Gheat_u_sunlit, double Gheat_u_shaded, double Gheat_g,
                   double lai_o_sunlit, double lai_o_shaded, double lai_u_sunlit, double lai_u_shaded,
                   double* SH_o, double* SH_u, double* SH_g)
{
    double SH_o_sunlit, SH_o_shaded, SH_u_sunlit, SH_u_shaded; // sensible heat from leaves
    double meteo_pack_output[10];
    double density_air0, cp_air, vpd;

    meteo_pack (temp_air, rh_air, meteo_pack_output);
    density_air0 =meteo_pack_output[1];
    cp_air = meteo_pack_output[2];
    vpd = meteo_pack_output[3];

    /********************************************/
    SH_o_sunlit=(tempL_o_sunlit-temp_air)*density_air0*cp_air*Gheat_o_sunlit;
    SH_o_shaded=(tempL_o_shaded-temp_air)*density_air0*cp_air*Gheat_o_shaded;

    SH_u_sunlit=(tempL_u_sunlit-temp_air)*density_air0*cp_air*Gheat_u_sunlit;
    SH_u_shaded=(tempL_u_shaded-temp_air)*density_air0*cp_air*Gheat_u_shaded;
    /********************************************/

    *SH_o = SH_o_sunlit*lai_o_sunlit+SH_o_shaded*lai_o_shaded;
    *SH_u = SH_u_sunlit*lai_u_sunlit+SH_u_shaded*lai_u_shaded;

    *SH_o=max(-200, *SH_o);
    *SH_u=max(-200, *SH_u);

    *SH_g=(temp_g-temp_air)*density_air0*cp_air*Gheat_g;
}
