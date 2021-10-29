/// @file netRadiation.c
/// @brief This module will calculate net radiation at both canopy level and leaf level
/// @author Edited by XZ Luo
/// @date May 23, 2015


# include "beps.h"


/// @brief Function to calculate net radiation at canopy level and leaf level
/// @details [input] global solar radiation,
///                  cosine value for solar zenith angle, albedo of leaves
///                  albedo of snow, percentage of snow cover,
///                  leaf area index overstorey and understorey,
///                  temperature of overstorey, understorey and ground (contain snow?)
///                  temperature of air (C), relative humidity (0-100)
/// @details [output] net radiation for canopy, overstorey, understorey and ground;
///                   net radiation on sunlit, shaded leaves of overstorey and understorey.
/// @param  shortRad_global           global short radiation
/// @param  CosZs                     cosine value of solar zenith angle
/// @param  temp_o                    temperature of overstorey
/// @param  temp_u                    temperature of understory
/// @param  temp_g                    temperature of ground
/// @param  lai_o                     leaf area index of overstory, without stem
/// @param  lai_u                     leaf area index of understory, without stem
/// @param  lai_os                    leaf area index of overstory, with stem
/// @param  lai_us                    leaf area index of understory, with stem
/// @param  lai_o_sunlit              sunlit leaves LAI with consideration of stem, overstory
/// @param  lai_o_shaded              shaded leaves LAI with consideration of stem, overstory
/// @param  lai_u_sunlit              sunlit leaves LAI with consideration of stem, understory
/// @param  lai_u_shaded              shaded leaves LAI with consideration of stem, understory
/// @param  clumping                  clumping index
/// @param  temp_air                  air temperature
/// @param  rh                        relative humidity
/// @param  albedo_snow_v             albedo of snow in this step, visible
/// @param  albedo_snow_n             albedo of snow in this step, near infrared
/// @param  percentArea_snow_o        percentage of snow on overstorey (by area)
/// @param  percentArea_snow_u        percentage of snow on understorey (by area)
/// @param  percent_snow_g            percentage of snow on ground (by mass)
/// @param  albedo_v_o                albedo of overstory, visible, not considering snow, decided by land cover
/// @param  albedo_n_o                albedo of overstory, near infrared
/// @param  albedo_v_u                albedo of understory, visible
/// @param  albedo_n_u                albedo of understory, near infrared
/// @param  albedo_v_g                albedo of ground, visible
/// @param  albedo_n_g                albedo of ground, near infrared
/// @param  netRad_o                  net radiation on overstorey
/// @param  netRad_u                  net radiation on understorey
/// @param  netRad_g                  net radiation on ground
/// @param  netRadLeaf_o_sunlit       net radiation at the leaf level, overstory sunlit, for ET calculation
/// @param  netRadLeaf_o_shaded       net radiation at the leaf level, overstory shaded
/// @param  netRadLeaf_u_sunlit       net radiation at the leaf level, understory sunlit
/// @param  netRadLeaf_u_shaded       net radiation at the leaf level, understory shaded
/// @param  netShortRadLeaf_o_sunlit  net shortwave radiation at leaf level, overstory sunlit, for GPP calculation
/// @param  netShortRadLeaf_o_shaded  net shortwave radiation at leaf level, overstory shaded
/// @param  netShortRadLeaf_u_sunlit  net shortwave radiation at leaf level, understory sunlit
/// @param  netShortRadLeaf_u_shaded  net shortwave radiation at leaf level, understory shaded
/// @return void
void netRadiation(double shortRad_global, double CosZs, double temp_o, double temp_u, double temp_g,
                  double lai_o, double lai_u, double lai_os, double lai_us, // LAI of overstorey and understorey, with and without stem
                  double lai_o_sunlit, double lai_o_shaded, double lai_u_sunlit, double lai_u_shaded,
                  double clumping, double temp_air, double rh,
                  double albedo_snow_v, double albedo_snow_n, double percentArea_snow_o, double percentArea_snow_u, double percent_snow_g,
                  double albedo_v_o, double albedo_n_o, double albedo_v_u, double albedo_n_u, double albedo_v_g, double albedo_n_g,
                  double* netRad_o, double* netRad_u, double* netRad_g,
                  double* netRadLeaf_o_sunlit, double* netRadLeaf_o_shaded, double* netRadLeaf_u_sunlit, double* netRadLeaf_u_shaded,
                  double* netShortRadLeaf_o_sunlit, double* netShortRadLeaf_o_shaded, double* netShortRadLeaf_u_sunlit, double* netShortRadLeaf_u_shaded)
{
    double netShortRad_o, netShortRad_u, netShortRad_g; // net short wave radiation on overstorey, understorey and ground
    double netShortRad_o_dir, netShortRad_o_df,netShortRad_u_dir, netShortRad_u_df, netShortRad_g_dir, netShortRad_g_df; // direct and diffuse part of net solar radiation
    double shortRad_dir, shortRad_df; //direct and diffuse radiation on top of the canopy
    double netLongRadLeaf_o_sunlit, netLongRadLeaf_o_shaded, netLongRadLeaf_u_sunlit, netLongRadLeaf_u_shaded;

    double netLongRad_o, netLongRad_u, netLongRad_g; // net long wave radiation
    double shortRadLeaf_o_dir, shortRadLeaf_u_dir, shortRadLeaf_o_df, shortRadLeaf_u_df;
    double albedo_o, albedo_u, albedo_g; //albedo of overstorey, understorey and ground (considering snow)
    double albedo_v_os, albedo_n_os, albedo_v_us, albedo_n_us, albedo_v_gs, albedo_n_gs; // albedo of three parts in visible and NIR band, (considering snow)

    double e_actual; // saturated and actual water vapor potential
    double meteo_pack_output[10]; // used to get actual vapor pressure

    double emissivity_air, emissivity_o, emissivity_u, emissivity_g; //emissivity of air, overstorey, understorey and ground
    double longRad_air, longRad_o, longRad_u, longRad_g; // long wave radiation emitted by different part
    double sb_constant=5.67 / 100000000; // stephen-boltzman constant
    double cosQ_o, cosQ_u; // indicators to describe leaf distribution angles in canopy. slightly related with LAI
    double gap_o_dir, gap_u_dir, gap_o_df, gap_u_df; //gap fraction of direct and diffuse radiation for overstorey and understorey (diffuse used for diffuse solar radiation and longwave radiation)
    double gap_os_dir, gap_us_dir, gap_os_df, gap_us_df; // like above, considering stem
    double ratio_cloud; // a simple ratio to differentiate diffuse and direct radiation


    // calculate albedo of canopy in this step
    albedo_v_os=albedo_v_o*(1-percentArea_snow_o)+albedo_snow_v*percentArea_snow_o;
    albedo_n_os=albedo_n_o*(1-percentArea_snow_o)+albedo_snow_n*percentArea_snow_o;
    albedo_v_us=albedo_v_u*(1-percentArea_snow_u)+albedo_snow_v*percentArea_snow_u;
    albedo_n_us=albedo_n_u*(1-percentArea_snow_u)+albedo_snow_n*percentArea_snow_u;

    albedo_o=0.5*(albedo_v_os+albedo_n_os);
    albedo_u=0.5*(albedo_v_us+albedo_n_us);

    // calculate albedo of ground in this step
    albedo_v_gs=albedo_v_g*(1-percent_snow_g)+albedo_snow_v*percent_snow_g;
    albedo_n_gs=albedo_n_g*(1-percent_snow_g)+albedo_snow_n*percent_snow_g;
    albedo_g=0.5*(albedo_v_gs+albedo_n_gs);

    // separate global solar radiation into direct and diffuse one
    if (CosZs < 0.001) // solar zenith angle small, all diffuse radiation
        ratio_cloud=0;
    else
        ratio_cloud=shortRad_global/(1367*CosZs);

    if (ratio_cloud > 0.8)
        shortRad_df=0.13*shortRad_global;
    else
        shortRad_df=(0.943 + 0.734 *ratio_cloud - 4.9*pow((ratio_cloud),2)+ 1.796 *pow((ratio_cloud), 3 ) + 2.058*pow ((ratio_cloud),4))*shortRad_global;

    shortRad_df=min(shortRad_df,shortRad_global);
    shortRad_df=max(shortRad_df,0);

    shortRad_dir=shortRad_global-shortRad_df;

    // fraction at each layer of canopy, direct and diffuse. use Leaf only lai here
    gap_o_dir=exp(-0.5*clumping*lai_o/CosZs);
    gap_u_dir=exp(-0.5*clumping*lai_u/CosZs);

    gap_os_dir=exp(-0.5*clumping*lai_os/CosZs); // considering stem
    gap_us_dir=exp(-0.5*clumping*lai_us/CosZs);

    cosQ_o=0.537+0.025*lai_o;
	cosQ_u=0.537+0.025*lai_u;
    gap_o_df=exp(-0.5*clumping*lai_o/cosQ_o);
    gap_u_df=exp(-0.5*clumping*lai_u/cosQ_u);

    gap_os_df=exp(-0.5*clumping*lai_os/cosQ_o); // considering stem
    gap_us_df=exp(-0.5*clumping*lai_us/cosQ_u);

    // emissivity of each part
    meteo_pack (temp_air, rh, meteo_pack_output);
    e_actual=meteo_pack_output [7];

    emissivity_air =1-exp(-(pow(e_actual*10.0,(temp_air+273.15)/1200.0)));
	emissivity_air =min(1,emissivity_air);
    emissivity_air =max(0.7,emissivity_air);

    emissivity_o=0.98;
    emissivity_u=0.98;
    emissivity_g=0.96;

    // net short direct radiation on canopy and ground
    if (shortRad_global>zero && CosZs >zero ) // only happens in day time, when sun is out
        {
            netShortRad_o_dir=shortRad_dir*((1-albedo_o)-(1-albedo_u)*gap_o_dir);
            netShortRad_u_dir=shortRad_dir*gap_o_dir*((1-albedo_u)-(1-albedo_g)*gap_u_dir);
            netShortRad_g_dir=shortRad_dir*gap_o_dir*gap_u_dir*(1-albedo_g);
        }
        else
        {
            netShortRad_o_dir=0;
            netShortRad_u_dir=0;
            netShortRad_g_dir=0;
        }

    // net short diffuse radiation on canopy and ground
    if (shortRad_global>zero && CosZs >zero ) // only happens in day time, when sun is out
        {
            netShortRad_o_df=shortRad_df*((1-albedo_o)-(1-albedo_u)*gap_o_df)+0.21*clumping*shortRad_dir*(1.1-0.1*lai_o)*exp(-CosZs);
            netShortRad_u_df=shortRad_df*gap_o_df*((1-albedo_u)-(1-albedo_g)*gap_u_df)+0.21*clumping*shortRad_dir*gap_o_dir*(1.1-0.1*lai_u)*exp(-CosZs);
            netShortRad_g_df=shortRad_df*gap_o_df*gap_u_df*(1-albedo_g);
        }
        else
        {
            netShortRad_o_df=0;
            netShortRad_u_df=0;
            netShortRad_g_df=0;
        }

    // total net shortwave radiation at canopy level
    netShortRad_o=netShortRad_o_dir+netShortRad_o_df;
    netShortRad_u=netShortRad_u_dir+netShortRad_u_df;
    netShortRad_g=netShortRad_g_dir+netShortRad_g_df;


    // net long wave radiation on canopy and ground
    longRad_air=emissivity_air*sb_constant*pow((temp_air+273.15),4);
    longRad_o=emissivity_o*sb_constant*pow((temp_o+273.15),4);
    longRad_u=emissivity_u*sb_constant*pow((temp_u+273.15),4);
    longRad_g=emissivity_g*sb_constant*pow((temp_g+273.15),4);

    netLongRad_o=(emissivity_o*(longRad_air+longRad_u*(1-gap_u_df)+longRad_g*gap_u_df)-2*longRad_o )
            *(1-gap_o_df )+emissivity_o*(1-emissivity_u )*(1-gap_u_df )*(longRad_air*gap_o_df
            +longRad_o*(1-gap_o_df));

    netLongRad_u=(emissivity_u*(longRad_air*gap_o_df+longRad_o*(1-gap_o_df)+longRad_g )-2*longRad_u )
            *(1-gap_u_df )+(1-emissivity_g )*((longRad_air*gap_o_df+longRad_o*(1-gap_o_df ) )*gap_u_df
            +longRad_u*(1-gap_u_df))
            +emissivity_u*(1-emissivity_o )*(longRad_u*(1-gap_u_df )+longRad_g*gap_u_df)*(1-gap_o_df);

    netLongRad_g=emissivity_g*((longRad_air*gap_o_df+longRad_o*(1-gap_o_df ) )*gap_u_df
               +longRad_u*(1-gap_u_df) )
               -longRad_g+(1-emissivity_u )*longRad_g*(1-gap_u_df);

    // total net radiation for overstorey, understorey and ground.
    *netRad_o=netShortRad_o+netLongRad_o;
    *netRad_u=netShortRad_u+netLongRad_u;
    *netRad_g=netShortRad_g+netLongRad_g;

    // leaf level net radiation updated way
    // reference Chen 2012 clumping index paper


        if (shortRad_global>zero && CosZs >zero ) // only happens in day time, when sun is out
        {
            shortRadLeaf_o_dir=0.5*shortRad_dir/CosZs;
            shortRadLeaf_o_dir=min(shortRadLeaf_o_dir,0.7*1362);
            shortRadLeaf_u_dir=shortRadLeaf_o_dir;

            shortRadLeaf_o_df=(shortRad_df-shortRad_df*gap_os_df)/lai_os+0.07*shortRad_dir*(1.1-0.1*lai_os)*exp(-CosZs);
            shortRadLeaf_u_df=(shortRad_df*gap_o_df-shortRad_df*gap_o_df*gap_us_df)/lai_us // pay attention to the type of gap fraction used here
									+0.05*shortRad_dir*gap_o_dir*(1.1-0.1*lai_us)*exp(-CosZs);
        }
        else
        {
            shortRadLeaf_o_dir=0;
            shortRadLeaf_u_dir=0;
            shortRadLeaf_o_df=0;
            shortRadLeaf_u_df=0;
        }

	//overstorey sunlit leaves

    if (lai_o_sunlit>0)
		{
			*netShortRadLeaf_o_sunlit=(shortRadLeaf_o_dir+shortRadLeaf_o_df)*(1-albedo_o); // diffuse
			netLongRadLeaf_o_sunlit=netLongRad_o/lai_os; // leaf level net long
			*netRadLeaf_o_sunlit = *netShortRadLeaf_o_sunlit + netLongRadLeaf_o_sunlit;
		}
	else
		{
			*netShortRadLeaf_o_sunlit=(shortRadLeaf_o_dir+shortRadLeaf_o_df)*(1-albedo_o);
			netLongRadLeaf_o_sunlit=netLongRad_o;
			*netRadLeaf_o_sunlit = *netShortRadLeaf_o_sunlit + netLongRadLeaf_o_sunlit;
		}

	// overstorey shaded leaf
    if(lai_o_shaded>0)
		{
			*netShortRadLeaf_o_shaded=shortRadLeaf_o_df*(1-albedo_o); // diffuse
			netLongRadLeaf_o_shaded=netLongRad_o/lai_os;

			*netRadLeaf_o_shaded = *netShortRadLeaf_o_shaded + netLongRadLeaf_o_shaded;
		}
    else
		{
			*netShortRadLeaf_o_shaded=shortRadLeaf_o_df*(1-albedo_o); // diffuse
			netLongRadLeaf_o_shaded=netLongRad_o;
			*netRadLeaf_o_shaded = *netShortRadLeaf_o_shaded + netLongRadLeaf_o_shaded;
		}


	// understorey sunlit leaf
    if(lai_u_sunlit>0)
		{
			*netShortRadLeaf_u_sunlit=(shortRadLeaf_u_dir+shortRadLeaf_u_df)*(1-albedo_u);
			netLongRadLeaf_u_sunlit=netLongRad_u/lai_us;

			*netRadLeaf_u_sunlit = *netShortRadLeaf_u_sunlit + netLongRadLeaf_u_sunlit;
		}
    else
		{
			*netShortRadLeaf_u_sunlit=(shortRadLeaf_u_dir+shortRadLeaf_u_df)*(1-albedo_u);
			netLongRadLeaf_u_sunlit=netLongRad_u;

			*netRadLeaf_u_sunlit = *netShortRadLeaf_u_sunlit + netLongRadLeaf_u_sunlit;
		}

	// understorey shaded leaf
    if(lai_u_shaded>0)
		{
			*netShortRadLeaf_u_shaded=shortRadLeaf_u_df*(1-albedo_u);
			netLongRadLeaf_u_shaded=netLongRad_u/lai_us;

			*netRadLeaf_u_shaded = *netShortRadLeaf_u_shaded + netLongRadLeaf_u_shaded;
		}
    else
		{
			*netShortRadLeaf_u_shaded=shortRadLeaf_u_df*(1-albedo_u);
			netLongRadLeaf_u_shaded=netLongRad_u;

			*netRadLeaf_u_shaded = *netShortRadLeaf_u_shaded + netLongRadLeaf_u_shaded;
		}


//  leaf level net radiation: original way, use canopy radiation divided by LAI.
//  here calculate net radiation for all leaf component again use the original method.
//  These code could be commented, because not right, I just put it here to validate model.
//    if(lai_o_sunlit>0)
//        *netRadLeaf_o_sunlit = netShortRad_o_dir/lai_o_sunlit + netLongRad_o/lai_os;
//    else
//        *netRadLeaf_o_sunlit = netShortRad_o_dir + netLongRad_o;
//
//    if(lai_o_shaded>0)
//        *netRadLeaf_o_shaded = netShortRad_o_df/lai_o_shaded + netLongRad_o/lai_os;
//    else
//        *netRadLeaf_o_shaded = netShortRad_o_df + netLongRad_o;
//
//
//    if(lai_u_sunlit>0)
//        *netRadLeaf_u_sunlit = netShortRad_u_dir/lai_u_sunlit + netLongRad_u/lai_us;
//    else
//        *netRadLeaf_u_sunlit = netShortRad_u_dir + netLongRad_u;
//
//
//    if(lai_u_shaded>0)
//        *netRadLeaf_u_shaded = netShortRad_u_df/lai_u_shaded + netLongRad_u/lai_us;
//    else
//        *netRadLeaf_u_shaded =netShortRad_u_df + netLongRad_u;
}
