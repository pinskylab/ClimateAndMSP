# ClimateAndMSP
Scripts and data to test whether proactively preparing for climate change-driven shifts in species habitats produces effective ocean plans and whether it imposes tradeoffs.

These scripts are provided in the interests of open science. If you have questions or find errors, please let us know.

Contact:<br/>
Malin Pinsky<br/>
Rutgers University<br/>
[malin.pinsky@rutgers.edu](malin.pinsky@rutgers.edu)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3991884.svg)](https://doi.org/10.5281/zenodo.3991884)


## Software versions
- Software was executed on a laptop running Mac OS and on a Linux workstation running CentOS.
- [R version 3.6.1](https://www.r-project.org/). Main scripting language, with the following packages:
	- maps 3.3.0
	- mapdata 2.3.0
	- maptools 0.9-5
	- rgdal 1.4-4
	- rgeos 0.5-1
	- raster 3.0-2
	- sf 0.7-7
	- data.table 1.12.8
	- beanplot 1.2
	- RColorBrewer 1.1-2
	- lme4 1.1-21
- [R version 3.5.3](https://www.r-project.org/). This version of R was needed to run gurobi, part of the conservation planning algorithm. We used these packages:
	- prioritizr 4.1.4
	- gurobi 8.1-1
- [InVEST-3.7.0](http://releases.naturalcapitalproject.org/?prefix=invest/3.7.0/) and [InVEST-3.7.0.post212+hd9612c04a8dc](https://storage.googleapis.com/natcap-dev-build-artifacts/invest/jdouglass/3.7.0.post212%2Bhd9612c04a8dc/InVEST-3.7.0.post212%2Bhd9612c04a8dc-mac.zip). Used for generating wind and wave energy estimates. The latter development version does not appear to be available anymore, but presumably, the bugs have been fixed in the newer releases.
- [Gurobi 8.1.1](http://gurobi.com/) with an academic license. Used for solving conservation planning problems.

## Directory structure
- [code](code/): scripts that are run using data, data_dl, and output, produce temp, output, figures, and tables
- [data/natcapt](data/natcap/): has data and parameter files generated for this project. Only used for InVEST runs.
	- Note that the [species range projections](https://www.bco-dmo.org/dataset/753124)) from [Morley et al. 2018 PLOS ONE](https://doi.org/10.1371/journal.pone.0196127) are also needed. See below for how they are integrated.
- data_dl: has data that is easily downloaded (not tracked by git)
	- marineregions/World_EEZ_v10_20180221/: Flanders Marine Institute (2018). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 10. Available online at http://www.marineregions.org/ https://doi.org/10.14284/312
	- natcap/WaveEnergy/: Natural Capital project InVEST wave energy data 3.6.0 from [here](http://data.naturalcapitalproject.org/invest-data/3.6.0/WaveEnergy.zip)
	- natcap/WindEnergy/: Nature Capital project InVEST wind energy data 3.6.0 from [here](http://data.naturalcapitalproject.org/invest-data/3.6.0/WindEnergy.zip)
	- natcap/Marine/: Natural Capital project InVEST marine base data 3.6.0 from [here](http://data.naturalcapitalproject.org/invest-data/3.6.0/Marine.zip)
	- oceanadapt/: all-regions-full.rds from [OceanAdapt](https://github.com/pinskylab/OceanAdapt) version [update2019](https://zenodo.org/record/3890214)
	- sau/: Catch data v47-1 from [Sea Around Us](http://www.seaaroundus.org) for Large Marine Ecosystems 1-3, 5-9, 65
	- statcan/: [Canadian census population centres](https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/bound-limit-2016-eng.cfm), specifically ldpl000b16a_e
	- usgs/: [US National Map Cities and Towns](https://nationalmap.gov/small_scale/atlasftp.html#citiesx), specifically citiesx020_nt00007
	- WDPA/: World Database on Protected Areas, Aug 2019, from [Protected Planet](https://www.protectedplanet.net/marine)
- [figures](figures/): plotted figures
- [output](output/): files output by scripts and also used as input for some scripts (tracked by git)
	- [prioritizr_runs](output/prioritizr_runs/): the conservation plans produced by prioritizr
- [tables](tables/): output tables
- temp: temporary files output by scripts (not tracked by git)
- temp_figures: temporary figures output by scripts (not tracked by git)

## Analysis workflow (roughly... see below for more details)
- 0.2_averageSBT.r: to set up a climatology
- 1.0_processWDPA.r: grid the World Database on Protected Areas (WDPA)
- 1.1_summarizeWDPAbygrid.r: calculate the fraction of protected area in each grid cell
- 1.2_evalWDPA.r: compare protected areas to projected shifts in species habitats
- 1.3_WDPAstats.r: summarize the shifts
- 1.4_compare_WDPA_projections_data.r: Compare the projections to ecological survey observations
- 3.0_make_NatCap_files.r: make input files for the wind and wave energy calculations by InVEST
- 3.1_processNatCap.r: process the output from InVEST
- 5.0_define_CMSP.r: set up the ocean planning analysis
- 5.0.1_process_species_proj.r: summarize species habitat suitability on the ocean planning analysis grid
- 5.1_prioritizr.r: run the ocean planning analysis with prioritizr
- 5.2_evalprioritizr_ensemble.R: compare the ocean plans to the ensemble species projections
- 5.2_evalprioritizr.r: compare the ocean plans to the individual species projections
- 5.3_map_species_proj.r: optional maps of species habitat suitability
- 5.3_plot_and_stats_prioritizr_evaluation.r: some optional plots and analysis examining the ocean plans
- 6.0_WDPAnetworkstatsbyregion.R: statistics on the protected area networks
- 6.1_evalRandomNetworks.r: evaluate randomly selected protected areas against shifts in species distribution
- 6.2_turnover_by_CMSPgrid.R: community change in each ocean planning grid cell
- prioritizr_frontier.r: ocean planning runs with a limited total area
- [plot_figures.r](scripts/plot_figures.r) to make figures and some tables

## Scripts and files behind figures, tables, and other calculations
All figures and Table S1 are produced by [plot_figures.r](code/plot_figures.r). The nested lists below indicate which files are needed for each figure or table, as well as which scripts produce those files.
<br/>

1. General statistics from [plot_figures.r](code/plot_figures.r).
	1. temp/presmap\_\*\_rcp\*\_\*.csv.gz and temp/biomassmap\_\*\_rcp\*\_\*.csv.gz
		1. [5.0.1_process_species_proj.r](code/5.0.1_process_species_proj.r)
			1. Species distribution projections from [Morley et al. 2018 PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0196127), available from [BCO DMO](https://www.bco-dmo.org/dataset/753124)
	1. temp/wdpaturnbyMPAbymod.csv.gz
		1. [1.2_evalWDPA.r](code/1.2_evalWDPA.r)
			1. [output/wdpa_cov_by_grid0.05.csv.gz](output/wdpa_cov_by_grid0.05.csv.gz)
				1. [1.1_summarizeWDPAbygrid.r](code/1.1_summarizeWDPAbygrid.r)
					1. temp/wdpa_by_grid0.05_intersect.rds and temp/SPsf2.rds
						1. [1.0_processWDPA.r](code/1.0_processWDPA.r)
							1. Species distribution projections from Morley et al. 2018 PLOS ONE: see above
							1. dataDL/WDPA/WDPA_Aug2019_marine-shapefile/WDPA_Aug2019_marine-shapefile
			1. Species distribution projections from Morley et al. 2018 PLOS ONE: see above
1. [Figure 1](figures/Fig1_study_regions.png)
	1. [output/turnover_by_CMSPgrid.csv](output/turnover_by_CMSPgrid.csv)
		1. [6.2_turnover_by_CMSPgrid.R](code/6.2_turnover_by_CMSPgrid.R)
			1. temp/presmap\_\*\_rcp\*\_\*.csv.gz: see above
	1. [output/region_grid.csv.gz](output/region_grid.csv.gz)
		1. [5.0_define_CMSP.r](code/5.0_define_CMSP.r)
			1. temp/SPsf2.rds: see above
			1. dataDL/marineregions/World_EEZ_v10_20180221/eez_v10.shp
			1. dataDL/natcap/Marine/Land/global_polyline.shp
1. [Figure 2](figures/Fig2_prioritizr_goalsmetbymod.png)
	1. [output/goalsmetbymod_hist_all.csv](output/goalsmetbymod_hist_all.csv) and [output/goalsmetbymod_2per_all.csv](output/goalsmetbymod_2per_all.csv)
		1. [5.2_evalprioritizr.r](code/5.2_evalprioritizr.r)
			1. output/prioritizr_runs/solution\_\*.csv and output/prioritizr\_runs/spp\_\*.csv
				1. [5.1_prioritizr.r](code/5.1_prioritizr.r)
					1. temp/presmap\_\*\_rcp\*\_\*.csv.gz and temp/biomassmap\_\*\_rcp\*\_\*.csv.gz: see above
					1. Kappa thresholds from Morley et al. 2018 PLOS ONE available [here](https://raw.githubusercontent.com/pinskylab/project_velocity/master/output/modeldiag_Nov2017_fitallreg_2017.csv)
					1. [output/wind_npv.csv.gz](output/wind_npv.csv.gz) and [output/wave_npv.csv.gz](output/wave_npv.csv.gz)
						1. [3.1_processNatCap.r](code/3.1_processNatCap.r)
							1. ../NatCap_temp/westcoastwind/output/npv_US_millions.tif and ../NatCap_temp/eastcoastwind/output/npv_US_millions.tif
								1. Produced with InVEST 3.7.0 (see below)
									1. dataDL/natcap/WindEnergy/\*.\*
									1. dataDL/natcap/Marine/\*.\*
									1. temp/AOI_east.shp and temp/AOI_west.shp, made by [3.0_make_NatCap_files.r](code/3.0_make_NatCap_files.r), which needs temp/SPsf2.rds (see above)
							1. temp/SPsf2.rds: see above
							1. ../NatCap_temp/westcoastwave/output/npv_usd.tif and ../NatCap_temp/eastcoastwave/output/npv_usd.tif
								1. Produced with InVEST 3.7.0 (see below)
									1. temp/AOI_east.shp and temp/AOI_west.shp: see above
									1. dataDL/natcap/WaveEnergy/\*.\*
									1. dataDL/natcap/Marine/\*.\*
									1. [output/landgridpts_northamerica.csv](output/landgridpts_northamerica.csv), made by [3.0_make_NatCap_files.r](code/3.0_make_NatCap_files.r), which uses [data/natcap/NAmainland_lines.shp](data/natcap/NAmainland_lines.shp), dataDL/usgs/citiesx020_nt00007/citiesx020.shp, and dataDL/statcan/lpc_000b16a_e/lpc_000b16a_e.shp
									1. [data/natcap/Machine_Pelamis_Economic_Kimetal2012r0.05.csv](data/natcap/Machine_Pelamis_Economic_Kimetal2012r0.05.csv)
					1. [output/fishery_spps.csv](output/fishery_spps.csv)
						1. [5.0_define_CMSP.r](code/5.0_define_CMSP.r)
							1. Species distribution projections from Morley et al. 2018 PLOS ONE: see above
							1. dataDL/sau/\*.\*
					1. [output/region_grid.csv.gz](output/region_grid.csv.gz): see above
			1. Kappa thresholds from Morley et al. 2018 PLOS ONE: see above
			1. [output/region_grid.csv.gz](output/region_grid.csv.gz): see above
			1. temp/presmap\_\*\_rcp\*\_\*.csv.gz and temp/biomassmap\_\*\_rcp\*\_\*.csv.gz: see above
1. [Figure 3](figures/Fig3_planareas.png)
	1. output/prioritizr_runs/solution_*.csv: see above
1. [Figure 4](figures/Fig4_prioritizr_frontiers.png)
	1. temp/frontierall_2019-12-31_075440.csv
		1. [prioritizr_frontier.r](code/prioritizr_frontier.r)
			1. temp/presmap\_\*\_rcp\*\_\*.csv.gz and temp/biomassmap\_\*\_rcp\*\_\*.csv.gz: see above
			1. Kappa thresholds from Morley et al. 2018 PLOS ONE: see above
			1. [output/wind_npv.csv.gz](output/wind_npv.csv.gz) and [output/wave_npv.csv.gz](output/wave_npv.csv.gz): see above
			1. [output/fishery_spps.csv](output/fishery_spps.csv): see above
			1. [output/region_grid.csv.gz](output/region_grid.csv.gz): see above
	1. output/prioritizr_runs/spp_*.csv: see above
1. [Figure 5](figures/Fig5_MPA_turnover.png)
	1. temp/wdpaturnbyMPAbymod.csv.gz and temp/wdpaturnbynetbymod.csv.gz
		1. [1.2_evalWDPA.r](code/1.2_evalWDPA.r)
			1. Species distribution projections from Morley et al. 2018 PLOS ONE: see above
			1. [output/wdpa_cov_by_grid0.05.csv.gz](output/wdpa_cov_by_grid0.05.csv.gz): see above
1. [Figure S1](figures/FigS1_prioritizr_frontiers.png)
	1. temp/frontierall_2019-12-22_071607.csv and temp/frontierall_2019-12-31_075440.csv
		1. [prioritizr_frontier.r](code/prioritizr_frontier.r): see above
		1. output/prioritizr_runs/spp_*.csv: see above
1. [Figure S2](code/FigS2_randMPAs.png)
	1. [output/randMPAs_byBT.csv](output/randMPAs_byBT.csv)
		1. [6.1_evalRandomNetworks](code/6.1_evalRandomNetworks)
			1. [output/climatology.csv.gz](output/climatology.csv.gz)
				1. [0.2_averageSBT.r](code/0.2_averageSBT.r)
					1. Temperature projections from Morley et al. 2018 PLOS ONE: see above
			1. [output/region_grid.csv.gz](output/region_grid.csv.gz)
			1. Species distribution projections from Morley et al. 2018 PLOS ONE: see above
	1. [output/MPA_network_stats.csv](output/MPA_network_stats.csv)
		1. [6.0_WDPAnetworkstatsbyregion.R](code/6.0_WDPAnetworkstatsbyregion.R)
			1. [output/climatology.csv.gz](output/climatology.csv.gz): see above
			1. [output/wdpa_cov_by_grid0.05.csv.gz](output/wdpa_cov_by_grid0.05.csv.gz): see above
			1. temp/wdpaturnbyMPAbymod.csv.gz: see above
			1. [output/region_grid.csv.gz](output/region_grid.csv.gz): see above
1. [Figure S3](figures/FigS3_MPA_trawl.png)
	1. temp/wdpa_trawlsppobs_byreg.csv.gz and [output/MPAvstrawl_NPV_PPV.csv](output/MPAvstrawl_NPV_PPV.csv)
		1. [1.4_compare_WDPA_projections_data.r](code/1.4_compare_WDPA_projections_data.r)
			1. dataDL/oceanadapt/all-regions-full.rds
			1. temp/wdpa_by_spp_projnow.rds
				1. [1.2_evalWDPA.r](code/1.2_evalWDPA.r): see above
			1. [output/name_conversions_obs_to_proj.csv](output/name_conversions_obs_to_proj.csv)
				1. generated by early part of [1.4_compare_WDPA_projections_data.r](code/1.4_compare_WDPA_projections_data.r) in an interactive loop
			1. temp/wdpa_by_grid0.05_intersect.rds: see above
			1. dataDL/WDPA/WDPA_Aug2019_marine-shapefile/WDPA_Aug2019_marine-shapefile-polygons.shp: see above
			1. Kappa thresholds from Morley et al. 2018 PLOS ONE: see above
1. [Table S1](tables/tableS1.csv)
	1. output/prioritizr_runs/solution\_\*.csv and output/prioritizr_runs/spp\_\*.csv: see above
1. Table S4
	1. [output/fishery_spps.csv](output/fishery_spps.csv): see above

## InVEST runs
These are the detailed program settings used when InVEST was run.
- East coast wind
	- run with InVEST-3.7.0.post212+hd9612c04a8dc version
	- wind data points: dataDL/natcap/WindEnergy/input/Global_EEZ_WEBPAR_90pct_100ms.csv
	- area of interest: temp/AOI_east.shp: see above
	- Bathymetric DEM: dataDL/natcap/Marine/DEMs/global_dem/w001001.adf
	- Land polygon: dataDL/natcap/Marine/Land/global_polygon.shp
	- Global wind energy parameters: dataDL/natcap/WindEnergy/input/global_wind_energy_parameters.csv
	- turbine type: dataDL/natcap/WindEnergy/input/5_0_turbine.csv
	- number of turbines: 16 
	- min depth: 3m
	- max depth: 60 m
	- min distance: 0 m
	- max distance: 200000 m							
	- turn on valuation
		- foundation: $2.5M/turbine
		- discount rate 0.05
		- Grid connections: [none]
		- average shore to grid distance: 5 km
		- no price table
		- price of energy: $0.161/kWhr
		- rate of price change: 0.025
- West coast wind
	- run with InVEST-3.7.0.post212+hd9612c04a8dc version
	- wind data points: dataDL/natcap/WindEnergy/input/Global_EEZ_WEBPAR_90pct_100ms.csv
	- area of interest: temp/AOI_west.shp: see above
	- Bathymetric DEM: dataDL/natcap/Marine/DEMs/global_dem/w001001.adf
	- Land polygon: dataDL/natcap/Marine/Land/global_polygon.shp
	- Global wind energy parameters: dataDL/natcap/WindEnergy/input/global_wind_energy_parameters.csv
	- turbine type: dataDL/natcap/WindEnergy/input/5_0_turbine.csv
	- number of turbines: 16 
	- min depth: 3m
	- max depth: 60 m
	- min distance: 0 m
	- max distance: 200000 m							
	- turn on valuation
		- foundation: $2.5M/turbine
		- discount rate 0.05
		- Grid connections: [none]
		- average shore to grid distance: 5 km
		- no price table
		- price of energy: $0.161/kWhr
		- rate of price change: 0.025
- East coast wave
	- run with InVEST-3.7.0.post212+hd9612c04a8dc version
	- wave data: dataDL/natcap/WaveEnergy/input/WaveData
	- analysis area: Global
	- AOI: temp/AOI_east.shp: see above
	- Machine performance table: dataDL/natcap/WaveEnergy/input/Machine_Pelamis_Performance.csv
	- Machine parameter table: dataDL/natcap/WaveEnergy/input/Machine_Pelamis_Parameter.csv
	- Global DEM: dataDL/natcap/Marine/DEMs/global_dem/w001001.adf
	- Valuation: yes
		- Grid connections: output/landgridpts_northamerica.csv (see above)
		- Economic table: [data/natcap/Machine_Pelamis_Economic_Kimetal2012r0.05.csv](data/natcap/Machine_Pelamis_Economic_Kimetal2012r0.05.csv)
		- \# machines: 100
- West coast wave
	- run with InVEST-3.7.0.post212+hd9612c04a8dc version
	- wave data: dataDL/natcap/WaveEnergy/input/WaveData
	- analysis area: Global
	- AOI: temp/AOI_west.shp: see above
	- Machine performance table: dataDL/natcap/WaveEnergy/input/Machine_Pelamis_Performance.csv
	- Machine parameter table: dataDL/natcap/WaveEnergy/input/Machine_Pelamis_Parameter.csv
	- Global DEM: dataDL/natcap/Marine/DEMs/global_dem/w001001.adf
	- Valuation: yes
		- Grid connections: output/landgridpts_northamerica.csv (see above)
		- Economic table: [data/natcap/Machine_Pelamis_Economic_Kimetal2012r0.05.csv](data/natcap/Machine_Pelamis_Economic_Kimetal2012r0.05.csv)
		- \# machines: 100
