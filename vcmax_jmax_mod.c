/// @file vcamx_jmax_mode.c
/// @brief Subroutine to calculate the Vcmax and Jmax for sunlit and shaded big-leaf
/// @details [Reference] (1) Chen, J. M., G. Mo, J. Pisek, F. Deng, M. Ishozawa, D. Chan, 2012.
///                          Effects of foliage clumping on global terrestrial gross primary
///                          productivity. Global Biogeochemical Cycles, VOL. 26, GB1019, 18,
///                          doi:10.1029/2010GB003996
/// @details [Reference] (2) Medlyn, B.E. et al., 1999. Effects of elevated [CO2] on photosynthesis
///                          in European forest species: a meta-analysis of model parameters.
///                          Plant, Cell & Environment, 22(12): 1475-1495.
/// @authors Written and refactored by Liming He (liming.he@gmail.com)
/// @authors Original contributor: Gang Mo
/// @date  Last update: May 5, 2015
/// @date  First created: Mar. 18, 2014


#include "beps.h"


/// @brief Function to calculate the Vcmax and Jmax for sunlit and shaded leaf
/// @note Vcmax0 is for the leaf vcmax at top of the canopy. LHE
/// @note Just to clarify, in this version, Vcmax0 is still the average leaf Vcmax25,
///       Vcmax0 * slope_Vcmax_N * leaf_N is the top leaves Vcmax25. XL. 20190403.
/// @param  lai_o          overstory lai
/// @param  clumping       clumping index
/// @param  Vcmax0         maximum capacity of Rubisco at 25C-Vcmax
/// @param  slope_Vcmax_N  slope of Vcmax-N curve
/// @param  leaf_N         leaf Nitrogen content	mean value + 1 SD g/m2
/// @param  CosZs          cosine solar zenith angle
/// @param  Vcmax_sunlit   Vcmax of sunlit leaf
/// @param  Vcmax_shaded   Vcmax of shaded leaf
/// @param  Jmax_sunlit    Jmax of sunlit leaf
/// @param  Jmax_shaded    Jmax of shaded leaf
/// @return void
void Vcmax_Jmax(double lai_o, double clumping, double Vcmax0,
                double slope_Vcmax_N, double leaf_N, double CosZs,
                double *Vcmax_sunlit, double *Vcmax_shaded, double *Jmax_sunlit, double *Jmax_shaded)
{
    double Kn = 0.3;  // Kn = 0.713/2.4;
    double G_theta = 0.5;
    double K,expr1,expr2,expr3;

    if(lai_o < 0.001) lai_o = 0.001; // to avoid error in the code, a minimum is setup here. LHE.

    if (CosZs>0)
    {
        K = G_theta*clumping/CosZs; // always > 0 LHE

        expr1 = 1 - exp(-K*lai_o);  // always > 0 LHE
        expr2 = 1 - exp(-(Kn + K)*lai_o);
        expr3 = 1 - exp(-Kn*lai_o);

        if(expr1>0)
            *Vcmax_sunlit = Vcmax0 * slope_Vcmax_N * leaf_N * K*expr2 / (Kn + K) / expr1; // LHE. May 15, 2014@PGB205b

        else *Vcmax_sunlit = Vcmax0;

        if (K > 0 && lai_o > expr1/K)
            //Vcmax_shaded = Vcmax0*parameter1[47]*parameter1[46]*(expr3/Kn-expr2*clumping/(Kn+K))/(lai_o-expr1*2*CosZs); /* modified by Ting, Mar,28,2013 */

            *Vcmax_shaded = Vcmax0 * slope_Vcmax_N * leaf_N * (expr3 / Kn - expr2 / (Kn + K)) / (lai_o - expr1/K); // original code.
            //*Vcmax_shaded = Vcmax0 * slope_Vcmax_N * leaf_N * (expr3 / Kn - expr2 * clumping / (Kn + K)) / (lai_o - 2 * CosZs * expr1); // LHE May 15, 2014@PGB205b.

        else
            *Vcmax_shaded = Vcmax0;

    }
    else
    {
        *Vcmax_sunlit = Vcmax0;
        *Vcmax_shaded = Vcmax0;
    }

    *Jmax_sunlit = *Vcmax_sunlit * 2.39 - 14.2;
    *Jmax_shaded = *Vcmax_shaded * 2.39 - 14.2;
}
