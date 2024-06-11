# Pathways_PGS_Cardiometabolic_Risk_CVD
This repository includes the code used for the analysis in "Assessing genetic risk for cardiometabolic risk factors and cardiovascular disease in a multi-ethnic cohort of breast cancer survivors."

The code in this repository was written and prepared by Haiyang Sheng, Alexandra Zimbalist, and Peter Fiorica.

## `00_PGS_Calculation_Script.R`
This script was initially prepared by Haiyang Sheng. It calculates PGS for each individual from a list of input PGS from the PGS Catalog.

**Input:** 

  - Post-QC Genotypes in .gds format from the Pathways Study that are divided by ancestral population.
  - PGS format polygenic score file where the columns are `rsID`, `chr_name`, `effect_allele`, `other_allele`, `effect_weight`, `locus_name`, `OR`

**Output:**
  
  - Raw polygenic scores files for each chromosome and ancestral population, where each folder represents a polygenic score.
      - Folders/files follow a pattern such as `/PGS000123/EUR/chr14.scores.txt`
  - A file for each chromosome containing the number of variants the PGS uses from that chromosome. (`chr.14scores.txt.nvar`)
  - A file for each chromosome that reports the number of overlapping variants between the genotyped individual and the PGS. (`chr.14scores.txt.valid.nvar`)
  
## `01_merge.R` 
This script merges PGS information across chromosomes to calculate whole-genotype PGS for individuals.

**Input:** 

  - A  series of folders where each folder represents a given polygenic score. Within the given PGS folder, there are subfolders corresponding to genetic ancestry populations. Files for each chromosome for the genotyped individual within their subpopulation are located within these ancestry folders.

**Output:** 

  - A compiled file (.txt/xlsx) containing all PGS for all individuals where rows are individuals and columns correspond to different PGS.

## `02_merge.valid.nvar`
This script merges variant information across chromosomes to calculate PGS summary statistics for missingness.
**Input:** 

  - The `chr.*scores.txt.valid.nvar` and `chr.*scores.txt.nvar` files from `01_merge.R` that are needed to assure that >99% of the PGS variants are present in the genotyped individuals.

**Output**

  - A summary fie (.txt/xlsx) that provides summary level information for all PGS reporting the number of present, missing, and flipped variants that are present in the PGS and genotyped group of individuals.

## `03_incident_tables_git.Rmd`
This script calculates Fine-Gray subdistribution hazard models to calculate the HR of an incident event within the SIRE subgroups and overall.

**Input:** 

  - `.rds` files for PGS Z-scores and event information within each SIRE group.

**Output:**

  - `.rds` file containing HR, P-values, and 95% CIs for each PGS-event association within and across SIRE subgroups

## `04_Figures_Tables.Rmd`
This script generates figures and tables for the manuscript.

**Input:** 

-  `.rds` output file from `03_incident_tables_git.Rmd`

**Output:**

- Figures(.png/pdf) and tables (.xlsx/csv) to be included in the manuscript/supplement

# Contact

For questions about the contents of this repository, please contact Peter Fiorica (peter.fiorica@roswellpark.org).
