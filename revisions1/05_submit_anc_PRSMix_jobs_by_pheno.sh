#!/bin/bash
# Directory containing the list of phenotypes that have significantly associated PGS files
OUTPUT_DIR="prsmix_lists"
# Path to the script to call
CALL_SCRIPT="05_call_ancestry_PRSMix.sh"


# List of possible ancestries
ancestries=("AFR" "EUR" "AMR" "EAS")

# Loop through all *_file_paths.txt files in the directory
for file in "$OUTPUT_DIR"/*_file_paths.txt; do
  # Extract the phenotype (outcome) name from the filename
  phenotype=$(basename "$file" | sed 's/_file_paths.txt//')

  # Loop through each ancestry type
  for ancestry in "${ancestries[@]}"; do
    # Submit a job using both phenotype and ancestry as arguments
    sbatch --job-name="${phenotype}_${ancestry}" \
           --export=PHENO="$phenotype",ANCESTRY="$ancestry" \
           "$CALL_SCRIPT" "$phenotype" "$ancestry"    
    # Optionally, print a confirmation message to the console
    echo "Submitted job for phenotype: $phenotype with ancestry: $ancestry"
  done
done

