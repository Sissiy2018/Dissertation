#!/bin/bash

# Loop through sample sizes from 50 to 500 in steps of 50
for sample_size in {25..500..25}
do
    # Define input pgen file and output file names based on the sample size
    pgen_file="chr22_${sample_size}_Ne10000.pgen"
    infer_output_file="fakemap_inferARG_${sample_size}_Ne10000.threads"
    argn_output_file="fakemap_inferARG_${sample_size}_Ne10000.argn"
    nodes_output_file="fakemap_inferARG_nodes_${sample_size}_Ne10000"

    # Run the threads infer command with the specified input and output files
    threads infer --pgen "$pgen_file" --map_gz fake20Mb_Ne10000.map --demography constant_pop.demo --out "$infer_output_file" --mutation_rate 1e-8

    # Convert the threads output to ARGN format
    threads convert --threads "$infer_output_file" --argn "$argn_output_file"

    # Run the Python script to get recombination nodes
    python3 get_recombination_nodes.py -a "$argn_output_file" -o "$nodes_output_file"

    # Print a message indicating the completion of the current sample size
    echo "Completed processing for sample size ${sample_size}"
done

echo "All processing runs completed!"
