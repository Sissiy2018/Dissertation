#!/bin/bash

# Loop through sample sizes
for sample_size in {50..500..50}
do
    # Define input VCF file and output file names based on the sample size
    vcf_file="chr22_Ne10000_seed123_sample${sample_size}.vcf"
    output_prefix="gap_constant_${sample_size}_Ne10000"

    # Run the plink2 command with the specified input and output files
    plink2 --vcf "$vcf_file" --max-alleles 2 --make-pgen --out "$output_prefix"

    # Print a message indicating the completion of the current sample size
    echo "Completed PLINK processing for sample size ${sample_size}"
done

echo "All PLINK processing runs completed!"

