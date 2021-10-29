/// @file s_coszs.c
/// @brief calculate cos_solar zenith angle Z
/// @author W. Ju
/// @date July, 2004

# include "beps.h"

/// @brief Function to calculate cosine solar zenith angle
/// @param  jday  date of year
/// @param  j     local time/UTC time
///               code needs to be edited according to time format
/// @param  lat   latitude of site
/// @param  lon   longitude of site
/// @param  CosZs cosine solar zenith angle
/// @return void
void s_coszs(short jday,short j,float lat,float lon,double* CosZs)
{
	double hr, Hsolar1,Delta,Lat_arc;
	
	Delta=0.006918-0.399912*cos(jday*2.0*3.1415926/365.0)+0.070257*sin(jday*2.0*3.1415926/365.0)
		-0.006758*cos(jday*4.0*3.1415926/365.0)+0.000907*sin(jday*4.0*3.1415926/365.0);
	/* delta is the declination angle of sun.*/
	
	hr =j*24.0/RTIMES+lon/15.0; // UTC time
    //hr =j*24.0/RTIMES; // local time
    if (hr>24) hr=24-hr;
	if (hr<0) hr=24+hr;
	
    Lat_arc=3.1415926*lat/180.0;
    Hsolar1 = (hr - 12.0) * 2.0 *3.1415926 / 24.0;         /* local hour angle in arc. */
	
	*CosZs = cos(Delta) * cos(Lat_arc) * cos(Hsolar1) + sin(Delta) * sin(Lat_arc);
}   
