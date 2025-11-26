#!/bin/bash
##SBATCH --partition=general-compute --qos=general-compute
##SBATCH --time=23:59:00
##SBATCH --nodes=1
##SBATCH --constraint=IB
##SBATCH --mem=32000
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=1
##SBATCH --job-name="04_PRSMIx+"
##SBATCH --output=logs_prsmix/R-%x_%j.out
##SBATCH --error=logs_prsmix/R-%x_%j.err

#module load gcc openmpi r-bundle-bioconductor/3.15-R-4.2.0

##set -x
cd /bd-fs-mnt/TenantShare/exthome/home/l578311	

LOG_DIR="logs_prsmix"
mkdir -p "$LOG_DIR"

phenos=$(cat sig_phenotypes.txt)
ancestries=("ALL" "EUR" "AFR" "EAS" "AMR")


for pheno_name in $phenos; do
  for ancestry in "${ancestries[@]}"; do
    pheno="${pheno_name}_phenotype.txt"
    
    echo "Running PRSMix for $pheno_name and $ancestry ..."
    
    log_file="${LOG_DIR}/no_r2_${ancestry}_${pheno_name}.log"

    Rscript 05_KPNC_run_PRSMix.R --pheno_file /bd-fs-mnt/TenantShare/exthome/home/l578311/pgs_sandbox/phenotypes/${ancestry}_${pheno} \
    	--covariate_file /bd-fs-mnt/TenantShare/exthome/home/l578311/pgs_sandbox/phenotypes/${ancestry}_${pheno_name}_covariates.txt \
    	--pgs_file /bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_lists/${pheno_name}_file_paths.txt \
    	--pgs_list /bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_lists/${pheno_name}_PGS_list.txt \
    	--pheno_name ${pheno_name} --isbinary TRUE \
    	--original_beta_files_list /bd-fs-mnt/TenantShare/exthome/home/l578311/janise_share/weights.txt \
    	--out /bd-fs-mnt/TenantShare/exthome/home/l578311/no_r2_prsmix_results/${ancestry}_${pheno_name} \
  	  > "$log_file" 2>&1

  done
done

