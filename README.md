# Pathways_PGS_Cardiometabolic_Risk_CVD

## Background

This repository includes the code used for the analysis in "Polygenic risk for Cardiometabolic and Cardiovascular Disease in a Multi-Ethnic Cohort of Breast Cancer Survivors." (2025)

The repository is largely divided into the initial submission analyses and the code required for the revisions.

### `initial_submission`

The code in this folder was written and prepared by Haiyang Sheng, Alexandra Zimbalist, and Peter Fiorica. It includes some of the initial calculations of polygenic scores, Fine-Grey subdistribution hazard ratios, and figures. Following the request for revisions, we made substantial changes in the next folder.

### `revisions1`

This folder includes code primarily written by Peter Fiorica and Ally Zimbalist. It is stark change from the code above because the PGS were calculated using [PGSC_CALC](https://github.com/PGScatalog/pgsc_calc) and combined with [PRSMix+](https://github.com/buutrg/PRSmix). Further stratification and additional requests from the reviewers were followed after these initial steps.

An important consider for this folder is that PGS calculation took place on the University at Buffalo's CCR (HPC). The remaining analyses were performed on Kaiser Permante's ACI, which is not an HPC. Because of this there is a change in code structure across the two platforms.

The main notes for these processes are available at `PRS_Processing_KPNC_ACI.Rmd`

## Contact

For questions about this repo, please contact Peter Fiorica at pnfioric\@buffalo.edu.
