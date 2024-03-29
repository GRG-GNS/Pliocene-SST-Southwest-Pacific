# Pliocene-SST-Southwest-Pacific
This repository holds the R Script and R Data associated with  Grant et al., 2023. 'Regional amplified warming in the Southwest Pacific during the mid-Pliocene (3.3–3.0 Ma) and future implications'.  https://doi.org/10.5194/egusphere-2023-108

This script uses raw and processed data to undertake analyses to produce the figures and results presented in the manuscript. Data files referenced in the manuscript and used here are found (https://doi.org/10.5281/zenodo.7935217)

Most recent file upload _15May2023

########### Data files in Grantetal2023_SWPacificPlioceneSST_RDATA.R.RData
#### Biomarkers (Table S2 of manuscript) contains alkenone and TEX index and calibrations

#### SST_mPWP is compiled mid-Pliocene SST for all sites (including references,age, calibrations and uncertainty)
### SST_mPWP_Tex # This includes TEX data (from Biomarkers) formatted separately for plotting purposes 

#### SST_2095.df includes model values and normalised to historic value (.hist) for SSP1, SSP2 and SSP3 scenarios for 2090-2100AD,
#with Southern Hemisphere seasonal range (JJA - winter, ann - annual mean, DJF- summer)

### MIS5e is single site SST data for various scenarios including MIS5e Cortese et al., 2013 (described in text) for plotting 

### PlioMIP.gridmeans is the grid means for PlioCore multi-model mean (Haywood et al., 2020)
### PlioMIP.sites is the site extraction of PlioCore multi-model mean (Haywood et al., 2020)
### PlioMIP.latmeans is the latitudinal means of PlioCore multi-model mean (Haywood et al., 2020) between 140E - 160W and 0.5N - 79.5S

### CMIP6_model is the SSP1,2,3 site SST extracted from CESM2 (Danabasoglu et al., 2020) and INM (Volodin et al., 2018) with respect to HadISST 

###########Terms and units
## Sea Surface Temperatures (SSTs) are in degrees Celsius. 
## Latitude are in degrees north. Longitude in degrees east.
## NZESM (New Zealand Earth system model; Williasm et al., 2016), UKESM (United Kingdom Earth System model; Sellar et al., 2019)
## mPWP (mid-Pliocene Warm Period 3.3 - 3.0 Ma)
## UK'37 - alkenone biomarker SST proxy 
## TEX - TEX86 biomarker SST proxy

########## License 
# Creative Commons Zero v1.0 Universal


