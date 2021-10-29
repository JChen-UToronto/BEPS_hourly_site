/// @file calc_temp_leaf.c
/// @brief Subroutine to calculate the sunlit and shaded leaf temperatures for overstory and understory leave.
/// @authors Written and refactored by Liming He (liming.he@gmail.com)
/// @authors Original contributor: Weimin Ju
/// @date Last update: Sept. 15, 2015
/// @date Created on May 15, 2015


#include "beps.h"


/// @brief Function to calculate leaf temperature four components
///        (sunlit and shaded leaves, overstory and understory)
/// @details [output] Tc_o_sunlit,Tc_o_shaded,Tc_u_sunlit,Tc_u_shaded
/// @param Tair                 air temperature
/// @param slope                the slope of saturation vapor pressure-temperature curve
/// @param psychrometer         psychrometer constant, 0.066 kPa K
/// @param VPD_air              vapor pressure deficit
/// @param Cp_ca                specific heat of moist air in kJ/kg/K
/// @param Gw_o_sunlit          total conductance for water from the intercellular space of the leaves
///                             to the reference height above the canopy, overstory, sunlit
/// @param Gw_o_shaded          ..., overstory, shaded
/// @param Gw_u_sunlit          ..., understory, sunlit
/// @param Gw_u_shaded          ..., understory, shaded
/// @param Gww_o_sunlit         total conductance for water from the surface of the leaves
///                             to the reference height above the canopy, overstory, sunlit
/// @param Gww_o_shaded         ..., overstory, shaded
/// @param Gww_u_sunlit         ..., understory, sunlit
/// @param Gww_u_shaded         ..., understory, shaded
/// @param Gh_o_sunlit          total conductance for heat transfer from the leaf surface
///                             to the reference height above the canopy, overstory, sunlit
/// @param Gh_o_shaded          ..., overstory, shaded
/// @param Gh_u_sunlit          ..., understory, sunlit
/// @param Gh_u_shaded          ..., understory, shaded
/// @param Xcs_o                the fraction of canopy covered by snow, overstory
/// @param Xcl_o                the fraction of canopy covered by liquid water, overstory
/// @param Xcs_u                the fraction of canopy covered by snow, understory
/// @param Xcl_u                the fraction of canopy covered by liquid water, understory
/// @param radiation_o_sun      net radiation on leaves, overstory, sunlit
/// @param radiation_o_shaded   net radiation on leaves, overstory, shaded
/// @param radiation_u_sun      net radiation on leaves, understory, sunlit
/// @param radiation_u_shaded   net radiation on leaves, understory, shaded
/// @param Tc_o_sunlit          the effective canopy temperature in Kalvin, overstory, sunlit
/// @param Tc_o_shaded          the effective canopy temperature in Kalvin, overstory, shaded
/// @param Tc_u_sunlit          the effective canopy temperature in Kalvin, understory, sunlit
/// @param Tc_u_shaded          the effective canopy temperature in Kalvin, understory, shaded
/// @return void
void Leaf_Temperatures(double Tair, double slope, double psychrometer, double VPD_air, double Cp_ca,
                       double Gw_o_sunlit, double Gw_o_shaded, double Gw_u_sunlit, double Gw_u_shaded,
                       double Gww_o_sunlit, double Gww_o_shaded, double Gww_u_sunlit, double Gww_u_shaded,
                       double Gh_o_sunlit, double Gh_o_shaded, double Gh_u_sunlit, double Gh_u_shaded,
                       double Xcs_o, double Xcl_o, double Xcs_u, double Xcl_u,
                       double radiation_o_sun, double radiation_o_shaded, double radiation_u_sun, double radiation_u_shaded,
                       double *Tc_o_sunlit, double *Tc_o_shaded, double *Tc_u_sunlit, double *Tc_u_shaded)
{
    *Tc_o_sunlit = Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca,
                                    Gw_o_sunlit, Gww_o_sunlit, Gh_o_sunlit, Xcs_o, Xcl_o, radiation_o_sun);

    *Tc_o_shaded = Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca,
                                    Gw_o_shaded, Gww_o_shaded, Gh_o_shaded, Xcs_o, Xcl_o, radiation_o_shaded);

    *Tc_u_sunlit = Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca,
                                    Gw_u_sunlit, Gww_u_sunlit, Gh_u_sunlit, Xcs_u, Xcl_u, radiation_u_sun);

    *Tc_u_shaded = Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca,
                                    Gw_u_shaded, Gww_u_shaded, Gh_u_shaded, Xcs_u, Xcl_u, radiation_u_shaded);
}



/// @brief Subroutine to calculate leaf temperature
/// @param Tair           air temperature
/// @param slope          the slope of saturation vapor pressure-temperature curve
/// @param psychrometer   psychrometer constant, 0.066 kPa K
/// @param VPD_air        vapor pressure deficit
/// @param Cp_ca          specific heat of moist air in kJ/kg/K
/// @param Gw             total conductance for water from the intercellular space of the leaves
///                       to the reference height above the canopy
/// @param Gww            total conductance for water from the surface of the leaves
///                       to the reference height above the canopy
/// @param Gh             total conductance for heat transfer from the leaf surface
///                       to the reference height above the canopy
/// @param Xcs            the fraction of canopy covered by snow
/// @param Xcl            the fraction of canopy covered by liquid water
/// @param radiation      net radiation on leaves
/// @return [double Tc] the effective canopy temperature in Kalvin
double Leaf_Temperature(double Tair, double slope, double psychrometer, double VPD_air, double Cp_ca,
                        double Gw, double Gww, double Gh, double Xcs, double Xcl, double radiation)
{
    double p_star;
    double Tc;

    double R = 1 / Gw + 1 / (Gww* (Xcs + Xcl));
    p_star = (Gw + Gww* (Xcs + Xcl) ) / psychrometer;

    Tc = Tair + (radiation - VPD_air * rho_a * Cp_ca * p_star) / (rho_a * Cp_ca * (Gh + slope * p_star) ); //temperature of sunlit leaves, overstorey

    Tc = max(Tair - 3.0, Tc);
    Tc = min(Tair + 5.0, Tc);

    return Tc;
}
