# SREP2017
Analysis scripts for "Uncertainty Quantification of Extratropical Forest Biomass in CMIP5 Models over the Northern Hemisphere"


- <b>Titile</b>  
Uncertainty Quantification of Extratropical Forest Biomass in CMIP5 Models over the Northern Hemisphere"
- <b>Authors</b>  
Cheng-En Yang, Jiafu Mao, Forrest M. Hoffman, Daniel M. Ricciuto, Joshua S. Fu, Chris D. Jones, and Martin Thurner

- Manuscript under reviewed in <I>Scientific Reports</I>


- <b>Contact information</b>  
Cheng-En Yang (<a href="mailto:cyang10@vols.utk.edu">cyang10@vols.utk.edu</a>)  
Department of Civil and Environmental Engineering,  
The University of Tennessee at Knoxville,  
Knoxville, TN, 37996-2313, USA  

- <b>Last update</b>: March 27, 2018  
---

- Folder hierarchy 

scripts/  
├── step1-timeslice/   
├── step2-climatemean/  
├── step3-GSWP3/  
├── step4-MPIBGIv3.txt  
├── step5-create_gridarea.ncl  
├── step5-GLC2000.txt  
├── step6-upscale/  
│   ├── Bbio/  
│   │   ├── max/  
│   │   └── min/  
│   ├── cLeaf/  
│   │   ├── max/  
│   │   └── min/  
│   ├── cWood/  
│   │   ├── max/  
│   │   └── min/  
│   ├── Tbio/  
│   │   ├── max/  
│   │   └── min/  
├── step7-compare_upscaled/  
│   ├── 1-forest_types/  
│   ├── 2-individual_type/  
│   ├── 3-meteorology/  
│   ├── 4-combine/  
│   ├── 5-stacked/  
│   └── 6-relative_error/  
├── step8-correlation/  
│   ├── all_forest_types/  
│   │   ├── modern/  
│   │   └── PI/  
├── step9-Taylor_diagram/  
├── step10-detail_pft/  
└── step11-Whittaker/  

30 directories, 225 files.  
Checksum file: SHA256SUM.  

- Descriptions:  
    * step1 and step2: Preparations for modeled outputs from ESMs.  
    * step3: Preparations for observational climate data from GSWP3.  
    * step4: Preparations for observational forest biomass data from MPI-BGIv3.  
    * step5: Preparations for observational land cover types from GLC2000.  
    * step6: Upscaling fine-resolution observational data sets for each forest compartment according to each ESM's grid cell resolution. Min/Max folders are for calculating the range of forest biomass considering the uncertainty.  
    * step7: Analysis of forest biomass and climate conditions during the contemporary and the pre-industrail periods for lumped and individual forest types.  
    * step8: Correltations between total forest biomass and climate conditions (precipitation and surface temperature).  
    * step9: Evaluating correlations and ratios between modeled and observed biomass for Taylor diagram.  
    * step10: Evaluations of HadGEM2-ES model outputs with detailed PFT-level forest biomass information provided.  
    * step11: Evaluation of forest biomass on climate space (Whittaker diagram).  
