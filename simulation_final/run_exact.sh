#!/bin/bash

# Loop through sample sizes from 50 to 500 in steps of 50
for sample_size in {25..500..25}
do
    # Define input file and output file names based on the sample size
    argn_output_file="chr22_${sample_size}_Ne10000.argn"
    nodes_output_file="exactARG_nodes_${sample_size}_Ne10000"

    # Run the Python script to get recombination nodes
    python3 get_recombination_nodes.py -a "$argn_output_file" -o "$nodes_output_file"

    # Print a message indicating the completion of the current sample size
    echo "Completed processing for sample size ${sample_size}"
done

echo "All processing runs completed!"