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
treatments=("rad_tx_yn" "anthra_yn" "tras_yn" "horm_yn")



for treat in "${treatments[@]}"; do
  if [ ! -d "/bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_results/${treat}" ]; then
    mkdir "/bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_results/${treat}"
    echo "created /bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_results/${treat}"
  else
    echo "/bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_results/${treat} already exists!"
  fi
  for pheno_name in $phenos; do  
    for i in 0 1; do 
    pheno="${pheno_name}_${treat}_${i}"
    
    echo "Running PRSMix for $pheno_name and $treat ..."
    
    log_file="${LOG_DIR}/${pheno}.log"

    Rscript 05_KPNC_run_PRSMix.R --pheno_file /bd-fs-mnt/TenantShare/exthome/home/l578311/pgs_sandbox/phenotypes/ALL_${pheno}_phenotype.txt \
    	--covariate_file /bd-fs-mnt/TenantShare/exthome/home/l578311/pgs_sandbox/phenotypes/ALL_${pheno}_covariates.txt \
    	--pgs_file /bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_lists/${pheno_name}_file_paths.txt \
    	--pgs_list /bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_lists/${pheno_name}_PGS_list.txt \
    	--pheno_name ${pheno_name} \
    	--original_beta_files_list /bd-fs-mnt/TenantShare/exthome/home/l578311/janise_share/weights.txt \
    	--treatment ${treat} \
    	--out /bd-fs-mnt/TenantShare/exthome/home/l578311/prsmix_results/${treat}/ALL_${pheno} \
  	  > "$log_file" 2>&1
	  
	  done
  done
done

