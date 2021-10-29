/// @file lai_cacl.c
/// @brief recalculate sunlit and shaded leaf area index
/// @author Bin Chen
/// @date May, 2015

#include "beps.h"

/// @brief Function to recalculate sunlit and shaded leaf area index
/// @param  stem_o        overstory woody area
/// @param  stem_u        understory woody area
/// @param  LC            land cover type
/// @param  CosZs         cosine solar zenith angle
/// @param  lai_o         overstory lai
/// @param  clumping      clumping index
/// @param  lai_u         understory lai
/// @param  lai_o_sunlit  overstory sunlit lai
/// @param  lai_o_shaded  overstory shaded lai
/// @param  lai_u_sunlit  understory sunlit lai
/// @param  lai_u_shaded  understory shaded lai
/// @param  PAI_o_sunlit  overstory sunlit lai
/// @param  PAI_o_shaded  overstory shaded lai
/// @param  PAI_u_sunlit  understory sunlit lai
/// @param  PAI_u_shaded  understory shaded lai
/// @return void
void lai2(double stem_o,double stem_u,int LC,double CosZs,double lai_o,double clumping,double lai_u,
          double* lai_o_sunlit,double* lai_o_shaded,double* lai_u_sunlit,double* lai_u_shaded,
          double* PAI_o_sunlit,double* PAI_o_shaded,double* PAI_u_sunlit,double* PAI_u_shaded)
{
	 if(CosZs>0)
		 *PAI_o_sunlit = 2*CosZs*(1-exp(-0.5*clumping*(lai_o+stem_o)/CosZs));
	 else
		 *PAI_o_sunlit = 0;

	 *PAI_o_shaded = (lai_o+stem_o)-*PAI_o_sunlit;
	 
	 if(CosZs>0)
		 *PAI_u_sunlit = 2*CosZs*(1-exp(-0.5*clumping*(lai_o+stem_o+lai_u+stem_u)/CosZs))-*PAI_o_sunlit;
	 else
		 *PAI_u_sunlit = 0;
	 
	 *PAI_u_shaded = (lai_u+stem_u)-*PAI_u_sunlit;
	 
	 if(CosZs>0)
		 *lai_o_sunlit = 2*CosZs*(1-exp(-0.5*clumping*lai_o/CosZs));
	 else
		 *lai_o_sunlit = 0;

     //*lai_o_shaded = max(0,lai_o-*PAI_o_sunlit);  // original
     *lai_o_shaded = max(0,lai_o-*lai_o_sunlit);  // edited by J. Leng
	 
	 if(CosZs>0)
		 *lai_u_sunlit = 2*CosZs*(1-exp(-0.5*clumping*(lai_o+lai_u)/CosZs))-*lai_o_sunlit;
	 else
		 *lai_u_sunlit = 0;

     //*lai_u_shaded = max(0,lai_u-*PAI_u_sunlit);  // original
     *lai_u_shaded = max(0,lai_u-*lai_u_sunlit);  // edited by J. Leng
}