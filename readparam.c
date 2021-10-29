/// @file readparam.c
/// @brief Set parameters according to land cover types
/// @author Gang Mo
/// @date Apr., 2011

#include "beps.h"

/// @brief Function to set parameters
/// @param  lc          land cover type
///                     1-ENF 2-DNF 6-DBF 9-EBF 13-Shrub 40-C4 Plants default-Others
/// @param  parameter1  parameter array
/// @return void
void readparam(short lc, double* parameter1)
{
    /*****  land cover types  *****/
    parameter1[4]=lc;

    switch(lc)
    {
        /*****  conifer evergreen forest-ENF  *****/
        case 1:
            parameter1[2]=0.62;		       /* clumping_index  */
            parameter1[6]=100;		       /* light_compensate_point                    30*/
            parameter1[7]=1000;		       /* light_saturation_point                    1000*/
            parameter1[8]=4.5;		       /* LAI_max_overstory                         3.3*/
            parameter1[9]=2.4;		       /* LAI_max_understory                        2.4*/
            parameter1[10]=0.01;	       /* LAI_min_overstory                         0.0*/
            parameter1[11]=0.01;	       /* LAI_min_understory                        0.0*/
            parameter1[12]=0.87;	       /* albedo_new_snow                           0.87*/
            parameter1[13]=0.94;	       /* albedo_new_snow_vis                       0.94*/
            parameter1[14]=0.8;		       /* albedo_new_snow_nir                       0.80*/
            parameter1[15]=100.0;	       /* density_now_snow                          100.0*/
            parameter1[16]=1.33;	       /* z00 roughness length for heat =0.1 canopy height=1.5 */
            parameter1[17]=2700;	       /* specific_heat_overstory                   2700.0*/
            parameter1[18]=35;		       /* mass_overstory(kg/m^2)                    25*/
            parameter1[19]=2700;	       /* specific_heat_understory                  2700.0*/
            parameter1[20]=10;		       /* mass_understory(kg/m^2)                   10*/
            parameter1[21]=0.6;		       /* root_depth(m)                             0.8*/
            parameter1[22]=0.035;	       /* albedo_canopy_vis                         0.09*/
            parameter1[23]=0.23;	       /* albedo_c_nir_observed                     0.19*/
            parameter1[24]=0.10;	       /* albedo_saturated_soil                     0.10*/
            parameter1[25]=0.35;	       /* albedo_dry_soil                           0.35*/
            parameter1[26]=0.5;		       /* drainage_class                            1.0*/
            parameter1[27]=0.95;	       /* decay_rate_of_root_distribution           0.97*/
            parameter1[28]=150;		       /* minimum_stomatal_resistance(s/m)          100*/
            parameter1[29]=20;		       /* canopy_height                             23*/
            parameter1[30]=3;		       /* height of understory                      3*/
            parameter1[31]=30;		       /* the_height_to_measure_wind_speed          37*/
            parameter1[32]=0.05;	       /* the_depth_litter                          0.05*/

            // the following parameters were added by MEM for the Ball and Berry calculations
            parameter1[33]=8;		       /* k Ball			*/
            parameter1[34]=0.0175;		   /* intercept_for_H2O_ball_berry		*/
            parameter1[35]=0.0175/1.6;	   /* intercept_for_C2O_ball_berry		*/
            parameter1[36]=62.5;		   /* maximum capacity of Rubisco at 25C-Vcmax	*/
            parameter1[37]=2.39*parameter1[36] - 14.2;		/* Jmax	at 25 C	*/
            parameter1[38]=28.0;	       /* mb / coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants*/

            parameter1[39]=0.0015/RTIMES;  /* leaf resp co.  kg C-1 d-1 kg-1    */
            parameter1[40]=0.0020/RTIMES;  /* stem resp co.   kg C-1 d-1 kg-1      */
            parameter1[41]=0.0020/RTIMES;  /* root resp co.   kg C-1 d-1 kg-1   */
            parameter1[44]=0.003/RTIMES;   /* fine root resp co.   kg C-1 d-1 kg-1	*/
            parameter1[45]=2.3;            /* Q10=2.3 constant for exp. resp.	*/

            // the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
            parameter1[46]=3.10+1.35;	   /* leaf Nitrogen content	mean value + 1 SD g/m2 */
            parameter1[47]=20.72/62.5;	   /* slope of Vcmax-N curve  		*/

            break;

        /*****  conifer deciduous forest-DNF  *****/
        case 2:
            parameter1[2]=0.68;		       /* clumping_index  */
            parameter1[6]=100;		       /* light_compensate_point                    30*/
            parameter1[7]=1000;		       /* light_saturation_point                    1000*/
            parameter1[8]=4.5;		       /* LAI_max_overstory                         3.3*/
            parameter1[9]=2.4;		       /* LAI_max_understory                        2.4*/
            parameter1[10]=0.01;	       /* LAI_min_overstory                         0.0*/
            parameter1[11]=0.01;	       /* LAI_min_understory                        0.0*/
            parameter1[12]=0.87;	       /* albedo_new_snow                           0.87*/
            parameter1[13]=0.94;	       /* albedo_new_snow_vis                       0.94*/
            parameter1[14]=0.8;		       /* albedo_new_snow_nir                       0.80*/
            parameter1[15]=100.0;	       /* density_now_snow                          100.0*/
            parameter1[16]=1.33;	       /* z00 roughness length for heat =0.1 canopy height=1.5 */
            parameter1[17]=2700;	       /* specific_heat_overstory                   2700.0*/
            parameter1[18]=35;		       /* mass_overstory(kg/m^2)                    25*/
            parameter1[19]=2700;	       /* specific_heat_understory                  2700.0*/
            parameter1[20]=10;		       /* mass_understory(kg/m^2)                   10*/
            parameter1[21]=0.6;		       /* root_depth(m)                             0.8*/
            parameter1[22]=0.035;	       /* albedo_canopy_vis                         0.09*/
            parameter1[23]=0.23;	       /* albedo_c_nir_observed                     0.19*/
            parameter1[24]=0.10;	       /* albedo_saturated_soil                     0.10*/
            parameter1[25]=0.35;	       /* albedo_dry_soil                           0.35*/
            parameter1[26]=0.5;		       /* drainage_class                            1.0*/
            parameter1[27]=0.95;	       /* decay_rate_of_root_distribution           0.97*/
            parameter1[28]=150;		       /* minimum_stomatal_resistance(s/m)          100*/
            parameter1[29]=20;		       /* canopy_height                             23*/
            parameter1[30]=3;		       /* height of understory                      3*/
            parameter1[31]=30;		       /* the_height_to_measure_wind_speed          37*/
            parameter1[32]=0.05;	       /* the_depth_litter                          0.05*/

            // the following parameters were added by MEM for the Ball and Berry calculations
            parameter1[33]=8;		       /* k Ball			*/
            parameter1[34]=0.0175;		   /* intercept_for_H2O_ball_berry		*/
            parameter1[35]=0.0175/1.6;	   /* intercept_for_C2O_ball_berry		*/
            parameter1[36]=39.1;	   	   /* maximum capacity of Rubisco at 25C-Vcmax	*/
            parameter1[37]=2.39*parameter1[36] - 14.2;		/* Jmax	at 25 C	*/
            parameter1[38]=28.0;	       /* mb / coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants*/

            parameter1[39]=0.0015/RTIMES;  /* leaf resp co.  kg C-1 d-1 kg-1    */
            parameter1[40]=0.0020/RTIMES;  /* stem resp co.   kg C-1 d-1 kg-1      */
            parameter1[41]=0.0020/RTIMES;  /* root resp co.   kg C-1 d-1 kg-1   */
            parameter1[44]=0.003/RTIMES;   /* fine root resp co.   kg C-1 d-1 kg-1	*/
            parameter1[45]=2.3;            /* Q10=2.3 constant for exp. resp.	*/

            // the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
            parameter1[46]=1.81+0.64;	   /* leaf Nitrogen content	mean value + 1 SD g/m2 */
            parameter1[47]=22.05/39.1;	   /* slope of Vcmax-N curve  		*/
            break;

        /*****  broadleaf deciduous forest-DBF  *****/
        case 6:
            parameter1[2]=0.7;             /* clumping_index  0.7 */
            parameter1[6]=100;             /* light_compensate_point                    30*/
            parameter1[7]=1000;            /* light_saturation_point                    1000*/
            parameter1[8]=4.5;             /* LAI_max_overstory                         3.3*/
            parameter1[9]=2.4;             /* LAI_max_understory                        2.4*/
            parameter1[10]=0.01;           /* LAI_min_overstory                         0.0*/
            parameter1[11]=0.01;           /* LAI_min_understory                        0.0*/
            parameter1[12]=0.87;           /* albedo_new_snow                           0.87*/
            parameter1[13]=0.94;           /* albedo_new_snow_vis                       0.94*/
            parameter1[14]=0.8;            /* albedo_new_snow_nir                       0.80*/
            parameter1[15]=100.0;          /* density_now_snow                          100.0*/
            parameter1[16]=1.53;           /* z00  roughness length for heat             */
            parameter1[17]=2700;           /* specific_heat_overstory                   2700.0*/
            parameter1[18]=40;             /* mass_overstory(kg/m^2)                    25*/
            parameter1[19]=2700;           /* specific_heat_understory                  2700.0*/
            parameter1[20]=10;             /* mass_understory(kg/m^2)                   10*/
            parameter1[21]=0.8;            /* root_depth(m)                             0.8*/
            parameter1[22]=0.04;           /* albedo_canopy_vis                         0.09*/
            parameter1[23]=0.25;           /* albedo_c_nir_observed                     0.19*/
            parameter1[24]=0.10;           /* albedo_saturated_soil                     0.10*/
            parameter1[25]=0.35;           /* albedo_dry_soil                           0.35*/
            parameter1[26]=0.5;            /* drainage_class                            1.0*/
            parameter1[27]=0.97;           /* decay_rate_of_root_distribution           0.97*/
            parameter1[28]=200;            /* minimum_stomatal_resistance(s/m)          100*/
            parameter1[29]=23;             /* canopy_height                             23*/
            parameter1[30]=3;              /* height of understory                      3*/
            parameter1[31]=30;             /* the_height_to_measure_wind_speed          37*/
            parameter1[32]=0.05;           /* the_depth_litter                          0.05*/

            // the following parameters were added by MEM for the Ball and Berry calculations
            parameter1[33]=8;		       /* k Ball			*/
            parameter1[34]=0.0175;		   /* intercept_for_H2O_ball_berry		*/
            parameter1[35]=0.0175/1.6;	   /* intercept_for_C2O_ball_berry		*/
            parameter1[36]=57.7;		   /* maximum capacity of Rubisco at 25C-Vcmax	*/
            parameter1[37]=2.39*parameter1[36] - 14.2;		/* Jmax	at 25 C	*/
            parameter1[38]=28.0;	       /* mb / coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants*/

            parameter1[39]=0.015/RTIMES;   /* leaf resp co.  kg C-1 d-1 kg-1    */
            parameter1[40]=0.0035/RTIMES;  /* stem resp co.   kg C-1 d-1 kg-1      */
            parameter1[41]=0.0025/RTIMES;  /* root resp co.   kg C-1 d-1 kg-1   */
            parameter1[44]=0.003/RTIMES;   /* fine root resp co.   kg C-1 d-1 kg-1	*/
            parameter1[45]=2.3;            /* Q10=2.3 constant for exp. resp.	*/

            // the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
            parameter1[46]=1.74+0.71;	   /* leaf Nitrogen content	mean value + 1 SD g/m2 */
            parameter1[47]=33.79/57.7;	   /* slope of Vcmax-N curve  		*/
            break;

        /*****  broadleaf evergreen forest-EBF  *****/
        case 9:
            parameter1[2]=0.63;            /* clumping_index  0.8 */
            parameter1[6]=100;             /* light_compensate_point                    30*/
            parameter1[7]=1000;            /* light_saturation_point                    1000*/
            parameter1[8]=4.5;             /* LAI_max_overstory                         3.3*/
            parameter1[9]=2.4;             /* LAI_max_understory                        2.4*/
            parameter1[10]=0.01;           /* LAI_min_overstory                         0.0*/
            parameter1[11]=0.01;           /* LAI_min_understory                        0.0*/
            parameter1[12]=0.87;           /* albedo_new_snow                           0.87*/
            parameter1[13]=0.94;           /* albedo_new_snow_vis                       0.94*/
            parameter1[14]=0.8;            /* albedo_new_snow_nir                       0.80*/
            parameter1[15]=100.0;          /* density_now_snow                          100.0*/
            parameter1[16]=1.53;           /* z00  roughness length for heat                */
            parameter1[17]=2700;           /* specific_heat_overstory                   2700.0*/
            parameter1[18]=40;             /* mass_overstory(kg/m^2)                    25*/
            parameter1[19]=2700;           /* specific_heat_understory                  2700.0*/
            parameter1[20]=10;             /* mass_understory(kg/m^2)                   10*/
            parameter1[21]=0.8;            /* root_depth(m)                             0.8*/
            parameter1[22]=0.04;           /* albedo_canopy_vis                         0.09*/
            parameter1[23]=0.25;           /* albedo_c_nir_observed                     0.19*/
            parameter1[24]=0.10;           /* albedo_saturated_soil                     0.10*/
            parameter1[25]=0.35;           /* albedo_dry_soil                           0.35*/
            parameter1[26]=0.5;            /* drainage_class                            1.0*/
            parameter1[27]=0.97;           /* decay_rate_of_root_distribution           0.97*/
            parameter1[28]=200;            /* minimum_stomatal_resistance(s/m)          100*/
            parameter1[29]=23;             /* canopy_height                             23*/
            parameter1[30]=3;              /* height of understory                      3*/
            parameter1[31]=30;             /* the_height_to_measure_wind_speed          37*/
            parameter1[32]=0.05;           /* the_depth_litter                          0.05*/

            // the following parameters were added by MEM for the Ball and Berry calculations
            parameter1[33]=8;		       /* k Ball		            	*/
            parameter1[34]=0.0175;		   /* intercept_for_H2O_ball_berry		*/
            parameter1[35]=0.0175/1.6;	   /* intercept_for_C2O_ball_berry		*/
            parameter1[36]=29;		       /* maximum capacity of Rubisco at 25C-Vcmax	*/
            parameter1[37]=2.39*parameter1[36] - 14.2;		/* Jmax	at 25 C	*/
            parameter1[38]=28.0;	       /* mb / coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants*/

            parameter1[39]=0.015/RTIMES;   /* leaf resp co.  kg C-1 d-1 kg-1    */
            parameter1[40]=0.0035/RTIMES;  /* stem resp co.   kg C-1 d-1 kg-1      */
            parameter1[41]=0.0025/RTIMES;  /* root resp co.   kg C-1 d-1 kg-1   */
            parameter1[44]=0.003/RTIMES;   /* fine root resp co.   kg C-1 d-1 kg-1	*/
            parameter1[45]=2.3;            /* Q10=2.3 constant for exp. resp.	*/

            // the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
            parameter1[46]=2.17+0.8;	   /* leaf Nitrogen content	mean value + 1 SD g/m2 */
            parameter1[47]=14.02/29.0;	   /* 	slope of Vcmax-N curve  		*/
            break;

        /*****  Shrub-SH  *****/
        case 13:
            parameter1[2]=0.7;             /*   clumping_index                          0.8 */
            parameter1[6]=100;             /* light_compensate_point                    30*/
            parameter1[7]=1000;            /* light_saturation_point                    1000*/
            parameter1[8]=3.3;             /* LAI_max_overstory                         3.3*/
            parameter1[9]=0.01;            /* LAI_max_understory                        2.4*/
            parameter1[10]=0.01;           /* LAI_min_overstory                         0.0*/
            parameter1[11]=0.01;           /* LAI_min_understory                        0.0*/
            parameter1[12]=0.87;           /* albedo_new_snow                           0.87*/
            parameter1[13]=0.94;           /* albedo_new_snow_vis                       0.94*/
            parameter1[14]=0.8;            /* albedo_new_snow_nir                       0.80*/
            parameter1[15]=100.0;          /* density_now_snow                          1000*/
            parameter1[16]=0.3;            /* z00 roughness length for heat              */
            parameter1[17]=2700;           /* specific_heat_overstory                   2700.0*/
            parameter1[18]=25;             /* mass_overstory(kg/m^2)                    25*/
            parameter1[19]=2700;           /* specific_heat_understory                  2700.0*/
            parameter1[20]=0;              /* mass_understory(kg/m^2)                   10*/
            parameter1[21]=0.5;            /* root_depth(m)                             0.8*/
            parameter1[22]=0.045;          /* albedo_canopy_vis                         0.09*/
            parameter1[23]=0.28;           /* albedo_c_nir_observed                     0.19*/
            parameter1[24]=0.10;           /* albedo_saturated_soil                     0.10*/
            parameter1[25]=0.35;           /* albedo_dry_soil                           0.35*/
            parameter1[26]=0.5;            /* drainage_class                            1.0*/
            parameter1[27]=0.95;           /* decay_rate_of_root_distribution           0.97*/
            parameter1[28]=100;            /* minimum_stomatal_resistance(s/m)          100*/
            parameter1[29]=4;              /* canopy_height                             23*/
            parameter1[30]=0;              /* understory_height                           */
            parameter1[31]=30;             /* the_height_to_measure_wind_speed          37*/
            parameter1[32]=0.05;           /* the_depth_litter                          0.05*/

            // the following parameters were added by MEM for the Ball and Berry calculations
            parameter1[33]=8;		       /* k Ball		    	*/
            parameter1[34]=0.0175;		   /* intercept_for_H2O_ball_berry		*/
            parameter1[35]=0.0175/1.6;	   /* intercept_for_C2O_ball_berry		*/
            parameter1[36]=61.7*0.5 +54*0.5;/* maximum capacity of Rubisco at 25C-Vcmax	*/
            parameter1[37]=2.39*parameter1[36] - 14.2;		/* Jmax	at 25 C	*/
            parameter1[38]=28.0;	       /* mb / coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants*/

            parameter1[39]=0.001/RTIMES;   /*  leaf resp co.  kg C-1 d-1 kg-1    */
            parameter1[40]=0.002/RTIMES;   /*  stem resp co.   kg C-1 d-1 kg-1      */
            parameter1[41]=0.0015/RTIMES;  /*  root resp co.   kg C-1 d-1 kg-1   */
            parameter1[44]=0.003/RTIMES;   /*  fine root resp co.   kg C-1 d-1 kg-1	*/
            parameter1[45]=2.3;            /* Q10=2.3 constant for exp. resp.	*/

            // the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
            parameter1[46]=(2.03+1.05+1.69+0.62)*0.5;		/* leaf Nitrogen content	mean value + 1 SD g/m2 */
            parameter1[47]=(32.09/61.7+33.14/54.0)*0.5;		/* 	slope of Vcmax-N curve  		*/
            break;

        /*****  C4 plants  *****/
        case 40:
            parameter1[2]=0.73;	    	  /* clumping_index                            0.8 */
            parameter1[6]=100;		      /* light_compensate_point                    30*/
            parameter1[7]=1000;		      /* light_saturation_point                    1000*/
            parameter1[8]=4.5;		      /* LAI_max_overstory                         3.3*/
            parameter1[9]=0.01;		      /* LAI_max_understory                        2.4*/
            parameter1[10]=0.01;	      /*  LAI_min_overstory                        0.0*/
            parameter1[11]=0.01;	      /* LAI_min_understory                        0.0*/
            parameter1[12]=0.87;  	      /* albedo_new_snow                           0.87*/
            parameter1[13]=0.94;	      /* albedo_new_snow_vis                       0.94*/
            parameter1[14]=0.8;		      /* albedo_new_snow_nir                       0.80*/
            parameter1[15]=100.0;   	  /* density_now_snow                          100.0*/
            parameter1[16]=0.04;	      /*  z00 roughness length for heat    ../10   */
            parameter1[17]=2700;    	  /* specific_heat_overstory                   2700.0*/
            parameter1[18]=2.5;	    	  /* mass_overstory(kg/m^2)                    25*/
            parameter1[19]=2700;	      /* specific_heat_understory                  2700.0*/
            parameter1[20]=0.1;		      /* mass_understory(kg/m^2)                   10*/
            parameter1[21]=0.3;		      /* root_depth(m)                             0.8*/
            parameter1[22]=0.055;	      /* albedo_canopy_vis                         0.09*/
            parameter1[23]=0.3;		      /* albedo_c_nir_observed                     0.19*/
            parameter1[24]=0.10;	      /* albedo_saturated_soil                     0.10*/
            parameter1[25]=0.35;	      /* albedo_dry_soil                           0.35*/
            parameter1[26]=0.5;		      /* drainage_class                            1.0*/
            parameter1[27]=0.95;	      /* decay_rate_of_root_distribution           0.97*/
            parameter1[28]=200;		      /* minimum_stomatal_resistance(s/m)          100*/
            parameter1[29]=4;		      /* canopy_height                             23*/
            parameter1[30]=0.1;		      /* height of understory                      3*/
            parameter1[31]=30;		      /* the_height_to_measure_wind_speed          37*/
            parameter1[32]=0.05;	      /* the_depth_litter                          0.05*/

            // the following parameters were added by MEM for the Ball and Berry calculations
            parameter1[33]=4;		      /* k Ball			*/
            parameter1[34]=0.0175;		  /* intercept_for_H2O_ball_berry		*/
            parameter1[35]=0.0175/1.6;	  /* intercept_for_C2O_ball_berry		*/
            parameter1[36]=30;		      /* maximum capacity of Rubisco at 25C-Vcmax	*/
            parameter1[37]=2.39*parameter1[36] - 14.2;		/* Jmax	at 25 C	*/
            parameter1[38]=28.0;	      /* mb / coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants*/

            parameter1[39]=0.001/RTIMES;  /* leaf resp co.  kg C-1 d-1 kg-1    */
            parameter1[40]=0.002/RTIMES;  /*  stem resp co.   kg C-1 d-1 kg-1      */
            parameter1[41]=0.0015/RTIMES; /* root resp co.   kg C-1 d-1 kg-1   */
            parameter1[44]=0.003/RTIMES;  /*  fine root resp co.   kg C-1 d-1 kg-1	*/
            parameter1[45]=2.3;           /* Q10=2.3 constant for exp. resp.	*/

            // the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
            parameter1[46]=2.375;		  /* use crop's leaf Nitrogen content	mean value + 1 SD g/m2 */
            parameter1[47]=0.31;		  /* use crop's	slope of Vcmax-N curve  		*/
            break;

        /*****  other land cover types  *****/
        default:
            parameter1[2]=0.8;		      /* clumping_index                            0.8 */
            parameter1[6]=100;		      /* light_compensate_point                    30*/
            parameter1[7]=1000;		      /* light_saturation_point                    1000*/
            parameter1[8]=4.5;		      /* LAI_max_overstory                         3.3*/
            parameter1[9]=0.01;		      /* LAI_max_understory                        2.4*/
            parameter1[10]=0.01;	      /* LAI_min_overstory                         0.0*/
            parameter1[11]=0.01;	      /* LAI_min_understory                        0.0*/
            parameter1[12]=0.87;	      /* albedo_new_snow                           0.87*/
            parameter1[13]=0.94;	      /* albedo_new_snow_vis                       0.94*/
            parameter1[14]=0.8;		      /* albedo_new_snow_nir                       0.80*/
            parameter1[15]=100.0;	      /* density_now_snow                          100.0*/
            parameter1[16]=0.04;	      /*  z00 roughness length for heat    ../10   */
            parameter1[17]=2700;	      /* specific_heat_overstory                   2700.0*/
            parameter1[18]=2.5;		      /* mass_overstory(kg/m^2)                    25*/
            parameter1[19]=2700;	      /* specific_heat_understory                  2700.0*/
            parameter1[20]=0.1;		      /* mass_understory(kg/m^2)                   10*/
            parameter1[21]=0.3;		      /* root_depth(m)                             0.8*/
            parameter1[22]=0.055;   	  /* albedo_canopy_vis                         0.09*/
            parameter1[23]=0.3;	    	  /* albedo_c_nir_observed                     0.19*/
            parameter1[24]=0.10;	      /* albedo_saturated_soil                     0.10*/
            parameter1[25]=0.35;	      /* albedo_dry_soil                           0.35*/
            parameter1[26]=0.5;		      /* drainage_class                            1.0*/
            parameter1[27]=0.95;	      /* decay_rate_of_root_distribution           0.97*/
            parameter1[28]=200;		      /* minimum_stomatal_resistance(s/m)          100*/
            parameter1[29]=4;		      /* canopy_height                             23*/
            parameter1[30]=0.1;		      /* height of understory                      3*/
            parameter1[31]=30;		      /* the_height_to_measure_wind_speed          37*/
            parameter1[32]=0.05;	      /* the_depth_litter                          0.05*/

            // the following parameters were added by MEM for the Ball and Berry calculations
            parameter1[33]=8;		      /* k Ball			*/
            parameter1[34]=0.0175;		  /* intercept_for_H2O_ball_berry		*/
            parameter1[35]=0.0175/1.6;	  /* intercept_for_C2O_ball_berry		*/
            parameter1[36]=78.2*0.5 + 100.7*0.5; /* maximum capacity of Rubisco at 25C-Vcmax	*/
            parameter1[37]=2.39*parameter1[36] - 14.2;		/* Jmax	at 25 C	*/
            parameter1[38]=28.0;	      /* mb / coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants*/

            parameter1[39]=0.001/RTIMES;  /* leaf resp co.  kg C-1 d-1 kg-1    */
            parameter1[40]=0.002/RTIMES;  /*  stem resp co.   kg C-1 d-1 kg-1      */
            parameter1[41]=0.0015/RTIMES; /* root resp co.   kg C-1 d-1 kg-1   */
            parameter1[44]=0.003/RTIMES;  /*  fine root resp co.   kg C-1 d-1 kg-1	*/
            parameter1[45]=2.3;           /* Q10=2.3 constant for exp. resp.	*/

            // the parameter 46 47 were added by G.Mo for the Vcmax - Nitrogen calculations
            parameter1[46]=(1.75+0.76+1.62+0.61)*0.5;		/* leaf Nitrogen content	mean value + 1 SD g/m2 */
            parameter1[47]=(45.29/78.2 + 62.75/100.7)*0.5;	/* 	slope of Vcmax-N curve  		*/
    }

    return;
}
