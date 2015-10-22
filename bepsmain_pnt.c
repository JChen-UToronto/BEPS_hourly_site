/*************************************************************************
  bepsmain_pnt.c 
   
  BEPS 4.01 for a ponit, to sumilate canbon fluxes, energe fluxes and soil water ...

 ****************************************************************************/

#include "beps.h"
#include "soil.h"

int main()
{
	int jday,rstep,dd,tt,flag,i;
	int landcover,soil_type;			
	float lai_p[366];		
	float lon,lat,nppyr,lai_yr;	
	float lcv,stypev,stv,swv,sdpv;
	float ccdv,cssdv,csmdv,cfsdv,cfmdv,csmv,cmv,csv,cpv;
	float rv,tv,hv,pv,wv,civ;
	float m_rad[366][25],m_tem[366][25],m_hum[366][25],m_pre[366][25],m_wind[366][25];
	double tem,hum,st,sw,snowdepth,temp_soil1,temp_soil[layer+1];
	double Ccd[5],Cssd[5],Csmd[5],Cfsd[5],Cfmd[5],Csm[5],Cm[5],Cs[5],Cp[5];
	double lai,clumping;

	char inp_dir[70],site[33],lc_fn[70],cp_fn[70],lai_fn[60],me_fn[60],outp_fn[80];
	FILE *lc_ptr,*cp_ptr,*laif_ptr,*me_ptr,*outp_ptr;

	double es,esd;
	double theta_vfc[layer+1],theta_vwp[layer+1],thermal_s[layer+1];
	double psi_sat[layer+1],bb[layer+1],fei[layer+1];
	double Ksat[layer+1],theta[layer+1];
	double coef[100]; 
	double CosZs;  
	double parameter[50];

	double var_o[40],var_n[40];
	double v2last[40];
	double outp[10],total[10];

	struct climatedata *meteo;
	struct results *mid_res;
	struct Soil *p_soil;

	meteo=(struct climatedata*)malloc(366*sizeof(struct climatedata));
	mid_res=(struct results*)malloc(366*sizeof(struct results));
	p_soil = (struct Soil*)malloc(sizeof(struct Soil));

	/*****  setup input/output file directory here  *****/
	sprintf(inp_dir, "input");
	sprintf(site, "p1");
	sprintf(lc_fn,"%s/%s_data1.txt",inp_dir,site);
	sprintf(cp_fn,"%s/p1_data2.txt",inp_dir);
	sprintf(lai_fn,"%s/%s_lai.txt",inp_dir,site);
	sprintf(me_fn,"%s/%s_meteo.txt",inp_dir,site);

	sprintf(outp_fn, "output/%s_01.txt",site);

	/*  ---  setting all accu. to 0   --- */
	for (i=0;i<=10;i++)   total[i]=0;
	
	/* Open landcover file */
  	if ((lc_ptr=fopen(lc_fn, "r")) == NULL) {
          printf("\n Unable to open data1 file,  exitting program ...\n\n");
	  exit(0);
    	}
	/* Read longitude, latitude  */
		fscanf(lc_ptr,"%f  %f ",&lon,&lat);
	/* Read landcover/clumping index data for each pix */
		fscanf(lc_ptr,"%f  %f ",&lcv,&civ);
			landcover=(int)lcv;
			clumping=civ;
	/* Read soiltxt data for each pix */
		fscanf(lc_ptr,"%f ",&stypev);
			soil_type=(int)stypev;
	/* Read soiltemp data for each pix */
		fscanf(lc_ptr,"%f ",&stv);
			st=stv;
	/* Read soilwater data for each pix */
		fscanf(lc_ptr,"%f ",&swv);
			sw=swv;
	/* Read snow depth data for each pix */
		fscanf(lc_ptr,"%f ",&sdpv);
			snowdepth=sdpv;

	fclose(lc_ptr);

	/* Open carban pools file */
        if ((cp_ptr=fopen(cp_fn, "r")) == NULL) {
          printf("\n Unable to open data2 file,  exitting program ...\n\n");
          exit(0);
        }

	/* Read ann_LAI and ann_npp for each pix */
		fscanf(cp_ptr,"%f %f ",&lai_yr,&nppyr);

	/* Read carbon pools data for each pix */
		fscanf(cp_ptr,"%f %f %f %f %f %f %f %f %f ",&ccdv,&cssdv,&csmdv,&cfsdv,&cfmdv,&csmv,&cmv,&csv,&cpv);
			Ccd[0]=ccdv*1000;
			Cssd[0]=cssdv*1000;
			Csmd[0]=csmdv*1000;
			Cfsd[0]=cfsdv*1000;
			Cfmd[0]=cfmdv*1000;
			Csm[0]=csmv*1000;
			Cm[0]=cmv*1000;
			Cs[0]=csv*1000;
			Cp[0]=cpv*1000;

	fclose(cp_ptr);

	/* open lai  climate files */
  	if ((laif_ptr=fopen(lai_fn, "r")) == NULL) 
	{
         printf("\n Unable to open laif file,  exitting program ...\n\n");
	  exit(0);
    	}
	
	if ((me_ptr=fopen(me_fn, "r")) == NULL) 
	{
       printf("\n Unable to open meteo. file,  exitting program ...\n\n");
       exit(0);
    	}

	/*  Read lai  data */
	for (jday=1; jday<=365; jday++) 
	{
	fscanf(laif_ptr,"%f ",&lai_p[jday]);

		for (rstep=1;rstep<=24;rstep++) 
		{
		  /* read climate data */ 
		  fscanf(me_ptr,"%d %d %f %f %f %f %f \n",&dd,&tt,&rv,&tv,&hv,&pv,&wv);
	
		  	m_rad[jday-1][rstep-1]=(float)rv;
			m_tem[jday-1][rstep-1]=(float)tv;
			m_hum[jday-1][rstep-1]=(float)hv;
			m_pre[jday-1][rstep-1]=(float)pv/1000;
			m_wind[jday-1][rstep-1]=(float)wv;

		}  // end of time loop
	}  // end of day loop
	fclose(laif_ptr);
	fclose(me_ptr);

	/* Read parameters according to land cover type */
	readparam(landcover,parameter);
			   
	readcoef(landcover,soil_type,coef);  // for soil repiration and NEP calculation 

	/* open output file for   */
        if ((outp_ptr=fopen(outp_fn, "w")) == NULL)
        {   printf("\n Unable to open file <%s>,  exitting ...\n\n",outp_fn);
                exit(0);
    	}


	printf("simulating ...  \n");

	/* DDDDD day loop bgn  DDDDD	*/
   	for (jday=1; jday<=365; jday++) 
   	{

		/* re-calc LAI & renew clump index */
		lai=lai_p[jday]*parameter[2]/clumping; 


		for (rstep=0;rstep<24;rstep++) 
		{
	 		if (jday==1 && rstep==0)  flag=0;
	 		else 	flag=1;

		meteo->Srad = m_rad[jday-1][rstep];
		meteo->temp = m_tem[jday-1][rstep];
		meteo->rain = m_pre[jday-1][rstep];
		meteo->wind = m_wind[jday-1][rstep];
		meteo->LR = -200.0;   //  -200.0 means no measured longwave radiation, the value will be calculated later 
	
		tem = m_tem[jday-1][rstep];
		hum = m_hum[jday-1][rstep];

		/* vapour pressure in mbar */ 
		es =0.46*hum*(tem+273.16)/100;	
		esd  =  6.1078 * exp((17.269*tem)/(237.3 + tem));
		/* Calculating Relative humidity */
		if (es/esd>=1)  meteo->rh = 100;
		else  meteo->rh = 100*es/esd; 

		//  meteo->rh = hum;  // when reading in relative humidity

		if(flag == 0)  // for 1st time step, to initialize var.
		{
		 Init_Soil_Parameters(landcover, soil_type, parameter[27], p_soil); // read soil parameters and set depth. LHE
		 p_soil->r_drainage = parameter[26];
		 Init_Soil_Status(p_soil, st, tem, sw, snowdepth); // LHE

		 for (i=0;i<=40;i++)   var_o[i] = 0;
		 for (i=3;i<=8;i++)   var_o[i] = tem;

		  for(i=9;i<=14;i++) var_o[i] = p_soil->temp_soil_p[i-9];
		  for(i=21;i<=26;i++) var_o[i] = p_soil->thetam_prev[i-21];
		  for(i=27;i<=32;i++) var_o[i] = p_soil->ice_ratio[i-27];

		}
		else   //  for other time step
		  for (i=0;i<=40;i++)   var_o[i] = v2last[i];
		
		/* ---  Calculating cos_solar zenith angle Z  --- */
		s_coszs(jday,rstep,lat,lon,&CosZs); 


		/* --- Starting simulation modules   --- */
		//printf("%d, %d, %f\n", jday, rstep, p_soil->thetam_prev[0]);
		inter_prg(jday,rstep,lai,clumping,parameter,meteo,CosZs,var_o,var_n,p_soil,mid_res);

//		printf("%d, %d, %f, %f\n", jday, rstep, p_soil->thetam_prev[0], p_soil->f_soilwater);

		for (i=0;i<=40;i++)  v2last[i]=var_n[i];

					
		/*  ---  plant respiration /NPP module  ---  */
		temp_soil1=p_soil->temp_soil_c[1]; 
		plantresp(landcover,mid_res,lai_yr,lai,tem,temp_soil1,CosZs); 
			   
		/*  ---  soil respiration module  ---  */

		soilresp(Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,Cp,nppyr,coef,soil_type,p_soil,mid_res); 

  		/*  ---   writing output   ---  */
		outp[1]=mid_res->GPP;
		outp[2]=mid_res->Trans+mid_res->Evap; 
		outp[3]=mid_res->NEP;   
		outp[4]=mid_res->npp_o + mid_res->npp_u;   


//	fprintf(outp_ptr,"%d %d gpp= %f tr= %f Ev= %f \n",jday,rstep,outp[1],outp[2],outp[3]);

		total[1]=total[1]+outp[1];
		total[2]=total[2]+outp[2]; 
		total[3]=total[3]+outp[3];   

		}/* End of hourly loop */
   	}/* End of daily loop  */		    	  
	
   	fprintf(outp_ptr,"total GPP: %f \t ET: %f \tNEP: %f \n",total[1],total[2],total[3]);

	printf("total GPP: %f \t ET: %f \tNEP: %f \n",total[1],total[2],total[3]);


  	fclose(outp_ptr);

  	free(meteo);
  	free(mid_res);
  	free(p_soil);
  
return 1;
}	/* end of main */

