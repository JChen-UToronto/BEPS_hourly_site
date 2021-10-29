@mainpage BEPS hourly version (v4.10)

The user guide for BEPS hourly version for site (v4.10)

-----------------------------------------
-----------------------------------------

This model was initially developed for boreal ecosystems and has been adapted for all ecosystems over the globe.\n
BEPS mechanistically includes the impacts of various drivers on gross primary productivity (GPP) (climate, CO2 concentration, and nitrogen deposition) and assimilates vegetation structure (LAI) data. \n
BEPS also simulates the dynamics of carbon pools beyond GPP and uses a spin-up procedure to prescribe soil carbon pools for estimating autotrophic respiration (AR) and heterotrophic respiration (HR).

The BEPS hourly version for site (v4.10) can be used in two ways:

1) Dependency import \n
Please copy the header file and source file into traditional IDEs (i.e. Code::block, https://www.codeblocks.org/) and directly build and run the model.

2）CMake \n
Please find the "CMakeLists.txt" file. The BEPS v4.10 model requires minimum 3.17 CMake version and is based on C99 standard. \n
It is recommended to use CLion (https://www.jetbrains.com/clion/) and MingW (https://www.mingw-w64.org/) to compile and run the model.

Make sure the "input" and "output" folders have been created in the current folder of the source codes.

According to users' research interests, the parameters and code structure can be edited. Please remember to make readable comment and git version control after each edition.

Please cite [ARTICLES] for using the BEPS model.

Please see "Modules_variables4BEPS.docx" for detailed parameter descriptions.

-----------------------------------------
-----------------------------------------

The BEPS model requires four input files: 1) Basic information; 2) Carbon pool data; 3) Leaf area index; 4) Meteorological data. \n
Users can find input data examples in the 'input' folder.

1) Basic information (data1 in the input data example)

long, lat, LC, CI, soiltxt, soiltemp, soilwater, snowdp [WITH TAB SPACE]

long --  the longitude of site

lat  --  the latitude of site

LC -- land cover type of site \n
1-ENF 2-DNF 6-DBF 9-EBF 13-shrub 40-C4 plants default-others

CI -- clumping index

soiltxt -- soil texture \n
1-land 2-loamy sand 3-sandy loam 4-loam 5-silty loam 6-sandy clay loam 7-clay loam 8-silty clay loam 9-sandy clay 10-silty clay 11-clay default-Others

soiltemp -- soil temperature

soilwater -- soil water content

snowdp -- snow depth

2) Carbon pool data

LAI_yr, ann_NPP,  ccd,  cssd,  csmd,  cfsd,  cfmd,  csm,  cm,  cs,  cp [WITH TAB SPACE]

3) Leaf area index

Daily float number LAI [WITH TAB SPACE]

4) Meteorological data

DOY, H, SW, TA, VPD/RH, P, WS [WITH TAB SPACE] [LINEBREAK EACH HOUR]

DOY -- day of year (1-365)

H -- hour of day (1-24)

SW -- shortwave radiation

TA -- air temperature

VPD/RH -- vapor pressure deficit OR humidity

P -- precipitation

WS -- wind speed

-----------------------------------------
-----------------------------------------

References for algorithms in this model：

He, L.; Wang, R.; Mostovoy, G.; Liu, J.; Chen, J.M.; Shang, J.; Liu, J.; McNairn, H.; Powers, J. Crop Biomass Mapping Based on Ecosystem Modeling at Regional Scale Using High Resolution Sentinel-2 Data. Remote Sens. 2021, 13, 806. https://doi.org/10.3390/rs13040806

Luo, X., Keenan, T.F. Global evidence for the acclimation of ecosystem photosynthesis to light. Nat Ecol Evol (2020). https://doi.org/10.1038/s41559-020-1258-7

Chen, J.M., Ju, W., Ciais, P. et al. Vegetation structural change since 1981 significantly enhanced the terrestrial carbon sink. Nat Commun 10, 4259 (2019). https://doi.org/10.1038/s41467-019-12257-8 (for daily version)

He, L., et al. (2019). "Diverse photosynthetic capacity of global ecosystems mapped by satellite chlorophyll fluorescence measurements." Remote Sensing of Environment 232: 111344.

He, L.; Mostovoy, G. Cotton Yield Estimate Using Sentinel-2 Data and an Ecosystem Model over the Southern US. Remote Sens. 2019, 11, 2000.

Luo, X, Croft, H, Chen, JM, He, L, Keenan, TF. Improved estimates of global terrestrial photosynthesis using information on leaf chlorophyll content. Glob Change Biol. 2019; 25: 2499– 2514. https://doi.org/10.1111/gcb.14624

Luo, X. Z., Chen, J. M., Liu, J. E., Black, T. A., Croft, H., Staebler, R., . . . McCaughey, H. (2018). Comparison of Big-Leaf, Two-Big-Leaf, and Two-Leaf Upscaling Schemes for Evapotranspiration Estimation Using Coupled Carbon-Water Modeling. Journal of Geophysical Research-Biogeosciences, 123(1), 207-225.

He, L. M., Chen, J. M., Gonsamo, A., Luo, X. Z., Wang, R., Liu, Y., & Liu, R. G. (2018). Changes in the Shadow: The Shifting Role of Shaded Leaves in Global Carbon and Water Cycles Under Climate Change. Geophysical Research Letters, 45(10), 5052-5061.

He, L. M., Chen, J. M., Croft, H., Gonsamo, A., Luo, X. Z., Liu, J. N., . . . Liu, Y. (2017). Nitrogen Availability Dampens the Positive Impacts of CO2 Fertilization on Terrestrial Ecosystem Carbon and Water Cycles. Geophysical Research Letters, 44(22), 11590-11600. doi:10.1002/2017gl075981

He, L., Chen, J. M., Liu, J., Bélair, S., & Luo, X. (2017). Assessment of SMAP soil moisture for global simulation of gross primary production. Journal of Geophysical Research: Biogeosciences, 122, doi:10.1002/2016JG003603. doi:10.1002/2016JG003603

Chen, B., Liu, J., Chen, J. M., Croft, H., Gonsamo, A., He, L., & Luo, X. (2016). Assessment of foliage clumping effects on evapotranspiration estimates in forested ecosystems. Agricultural and Forest Meteorology, 216, 82-92. doi:http://dx.doi.org/10.1016/j.agrformet.2015.09.017

He, L., Chen, J. M., Liu, J., Mo, G., Bélair, S., Zheng, T., . . . Barr, A. G. (2014). Optimization of water uptake and photosynthetic parameters in an ecosystem model using tower flux data. Ecological Modelling, 294(0), 94-104. doi:http://dx.doi.org/10.1016/j.ecolmodel.2014.09.019

Chen, J. M., Mo, G., Pisek, J., Liu, J., Deng, F., Ishizawa, M., & Chan, D. (2012). Effects of foliage clumping on the estimation of global terrestrial gross primary productivity. Global Biogeochemical Cycles, 26. doi:Artn Gb1019 Doi 10.1029/2010gb003996

Chen, B. Z., Chen, J. M., & Ju, W. M. (2007). Remote sensing-based ecosystem-atmosphere simulation scheme (EASS) - Model formulation and test with multiple-year data. Ecological Modelling, 209(2-4), 277-300. doi:DOI 10.1016/j.ecolmodel.2007.06.032

Ju, W., Chen, J. M., Black, T. A., Barr, A. G., Liu, J., & Chen, B. (2006). Modelling multi-year coupled carbon and water fluxes in a boreal aspen forest. Agricultural and Forest Meteorology, 140(1-4), 136-151. doi:10.1016/j.agrformet.2006.08.008

Chen, J. M., Liu, J., Cihlar, J., & Goulden, M. L. (1999). Daily canopy photosynthesis model through temporal and spatial scaling for remote sensing applications. Ecological Modelling, 124(2-3), 99-119. doi:Doi 10.1016/S0304-3800(99)00156-8

Liu, J., Chen, J. M., Cihlar, J., & Park, W. M. (1997). A process-based boreal ecosystem productivity simulator using remote sensing inputs. Remote Sensing of Environment, 62(2), 158-175. 

-----------------------------------------
-----------------------------------------

Compiled into doxygen format by Jiye Leng @UofT

Oct., 2021