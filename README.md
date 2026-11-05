# SPECTRA: Normalized Biovolume Size Spectra (NBSS) from ISIIS-3 and LISST data

R and Python code to calculate Normalized Biovolume Size Spectra (NBSS) from ISIIS-2,
ISIIS-3 (Cowen and Guigand, 2008) in-situ image data and LISST (Sequoia
Instruments) laser diffraction data, following the approach of Dugenne et
al. (2024) and the SPECTRA project (NCC Oregon coast), which serves as the main
template for the other three workflows. There are four total workflows,two from 
ISIIS-2 (MEZCAL and MBONISIIS-2), and ISIIS-3 (SPECTRA and MBON ISIIS-3). 
These separate workflows are necessary to accomodate, first, distinct transect and tow sampling design:
- MEZCAL tows are separated by day and night, by replicates, by Transect and by Cruise
- MBON tows are separated by day and night, by Transect and by Cruise
- SPECTRA tows are organized by Transect and Cruise.

Second, ISIIS-2 imagery is collected by only one camera, whereas ISIIS-3 imagery is collected by
two cameras, which implies extra processing to spatiotemporally merge the data from the cameras.
The tree structure below describes the general organization of the projects, with custom or specific
files to each workflow clarified througout the description.


## Repository Structure

```

├─ LICENSE                        GPL-3.0
├─ README.md                      This file
├─ requirements.txt               List of python dependencies needed
├─ R_requirements.R               R script for installing necessary dependencies
├─ ISIIS_NBSS_workflow_*
   ├── Makefile                       Workflow recipes (run via `make`)
   ├── ancillary/                     Supporting reference data (CSV)
   │   ├── shelf_break_200*.csv                       Shelf break longitudes
   │   ├── transect_lats*.csv                         Transect latitudes
   │   ├── transect_shelf_land*.csv                   Transect land boundaries
   │   ├── MEZCAL_FT_AND_GROUPINGv2.csv               Prob thresholds ISIIS-2
   │   ├── MEZCAL_GROUPS_AND_CFv2.csv                 Correction factors ISIIS-2 
   │   ├── MEZCAL_training_lib_sizes.csv              ISIIS-2 sizes from training library
   │   ├── contour_custom_ISIIS-2_size_classes.csv    Size class bins for ISIIS-2
   │   ├── SPECTRA_90_prob_filtering_manyDetritus_python.csv  Prob thresholds ISIIS-3
   │   ├── classification_metrics_90.0_manyDetritus_python.csv  to calculate ISIIS-3 correction factors
   │   ├── SPECTRA_training_lib_sizes.csv             ISIIS-3 sizes from training library
   │   ├── contour_custom_ISIIS-3_size_classes.csv    Size class bins for ISIIS-3
   │   └── LISST_size_classes_um.csv                  LISST size classes
   └── SCRIPTS/
       ├── functions/
       │   ├── ISIIS_PSS_binning_funcs.R             Core binning/NBSS functions
       │   ├── LISST_funcs.R                         LISST data functions
       │   └── model_functions.R                     Statistical model functions
       ├── *.R                                       Workflow scripts (see below)
       ├── Process_TarBall_PixelCount.py              Pixel counting from tarballs
       └── calculate_dark_pixels.py                   Dark pixel analysis
```

## Workflow

All scripts assume they are run from each workflow base directory (`ISIIS_NBSS_workflow_*`).

### 1. Image preprocessing (Python)

`scripts/Process_TarBall_PixelCount.py` extracts image tarballs and calls
`scripts/calculate_dark_pixels.py` to count light/dark pixels, compute contour
area and fitted ellipse parameters for each image crop. Output goes to `CSV/`.

### 2. Merge size & metadata (R)

| Script | Purpose |
|--------|---------|
| `0.1.0_get_tarfile_list.R` | Generate list of tar files to process |
| `1_combine_merge_size.R` | Combine ISIIS merged data (CTD/environmental) with per-particle size data |
| `2.0.1_combineCam_merged_data.R` | ISIIS-3:Merge CamTop and CamBottom into single combined data set |
| `2.0.2_get_singleCam_data.R` | ISIIS-3:Split combined camera data into non-overlapping per-camera records |
| `2_NEW_create_file_groups.R` | inspects time ranges and depth profiles of the combined data to generate a cast id index file that opens the combined data by chunks, each chunk representing a single cast |

### 3. Size binning and Normalized Biovolume Size Spectra (R)

`2.1_ISIIS_binned_abundance_data_7hrCorr.R` reads the `($CRUISE)_($TRANSECT_ID)_cast_id_depth_file_list_UPDATED.csv`, creates intermediate, clean summaries for the size spectra estimates, and calculates plankton abundance for each 1 m depth bin along a cast, both in terms of total individuals and biovolume per cubic meter 

`scripts/3_NEW_ISIIS_bin_and_NBSS.R` bins all detected particles by size
class, calculates biovolume per bin, and normalizes by bin width to produce
NBSS.

### 4. Kriging interpolation (R)

`scripts/4_Kriging_new_size_spectra.R` spatially interpolates (kriges) the
size spectra across each transect to produce continuous fields of particle
size distributions.

### 5. Normalized abundance plots (R)

`scripts/4.2_Normalized_abundance_plots.R` generates diagnostic plots of the
kriged size spectra.

### 6. LISST comparison (R)

`scripts/5_LISST_ISIIS_NBSS.R` merges LISST-derived NBSS with
ISIIS-derived NBSS for instrument comparison.

### 7. Supplementary scripts

| Script | Purpose |
|--------|---------|
| `supp_create_ISIIS-3_sizeClasses.R` | Generate ISIIS-2 and ISIIS-3 size class bin definitions |
| `supp_generate_bathy_polygon.R` | Create bathymetry polygon for spatial filtering |
| `supp_LISST_size_spectra_processing.R` | Calculate NBSS from LISST raw data |

### 8. Manuscript figures and tables

These are the scripts used for the figures and tables included in PROOCE-D-26-00232R1 
"Seston and plankton size spectra variability in a highly productive Eastern Boundary Current System"
currently under review in Progress in Oceanography  

| Script | Output |
|--------|--------|
| **Mapping and spatial visualization** | |
| `ManuscriptVar_Map_transects.R` | Manuscript transect map |
| `Map_transects_all.R` | Map of all sampling transects |
| **Statistical analyses and models** | |
| `ManuscriptVar_PSS_mixed_models_cast_ID_total_biovolume_only_exploratory_LMM.R` | Exploratory LMM of total biovolume with cast-level effects |
| `ManuscriptVar_PSS_mixed_models_environmental_R2_pooled_shelf_2x3.R` | Environmental LMM analysis and 2 × 3 panel figure |
| `ManuscriptVar_biovolume_vs_slope_16_LMM.R` | LMM analysis of biovolume–slope relationships |
| `ManuscriptVar_season_shelf_NBSS_three_responses_LMM.R` | LMM analysis of seasonal and shelf effects on three NBSS responses |
| **Model summaries and diagnostics** | |
| `ManuscriptVar_summarize_biovolume_vs_slope_16_LMM_effects.R` | Summary of biovolume–slope LMM effects |
| `ManuscriptVar_summarize_cast_model_residual_diagnostics.R` | Cast-level model residual diagnostics |
| `ManuscriptVar_summarize_environmental_LMM_effects.R` | Summary of environmental LMM effects |
| `ManuscriptVar_summarize_season_shelf_NBSS_three_responses_LMM.R` | Summary of seasonal and shelf LMM effects |
| **Kriging and spatial analyses** | |
| `ManuscriptVar_compile_kriging_cv_summary.R` | Compilation of kriging cross-validation results |
| `ManuscriptVar_multi_panel_kriged_figures.R` | Multi-panel kriged figures |
| `ManuscriptVar_multi_panel_kriged_forPres.R` | Kriged figures formatted for presentations |
| `ManuscriptVar_figTowYo.R` | Tow-level visualization and spatial plots |
| **Comparative analyses and supplementary tables** | |
| `ManuscriptVar_multi_panel_LISST_ISIIS_comparison.R` | Multi-panel LISST–ISIIS comparison |
| `ManuscriptVar_PSS_paired_source_comparison.R` | Paired source comparison |
| `ManuscriptVar_TableS5.R` | Supplementary Table S5 |



## Dependencies

### R

Run Rscript `R_requirements.R` 

### Python

See `requirements.txt`. Key packages include: `opencv-python`, `numpy`,
`pandas`, `pathlib`.

can be used as python3 -m pip install -r requirements.txt

## Usage

Each R script accepts command-line arguments documented via `-h`:

```bash
Rscript scripts/3_NEW_ISIIS_bin_and_NBSS.R -h
```

The `Makefile` documents the full workflow pipeline as named recipes (`make
<recipe>`). Edit the CRUISE and TRANSECT variables at the top of relevant
targets to match your data, and `TOW_ID` for all the projects except SPECTRA that 
contain multiple tows per transect per cruise. NOTE: these Makefiles are extensive, 
and test/legacy code still needs to be removed. The repostitory owner will do so in time.

## Citation

If you use this code, please cite the associated publication (to be added).

## License

GNU General Public License v3.0. See `LICENSE`.
