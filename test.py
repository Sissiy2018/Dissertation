import msprime
import numpy as np
import pandas as pd
import argparse
import gzip
import arg_needle_lib
import csv
def create_rate_map(map_file):
    # Initialize lists to store positions and rates
    positions = [0]
    gmaps = [0]

    # Read positions from map.txt file
    with open(map_file, 'r') as f:
        for line in f:
            if line.strip():  # Check if line is not empty
                parts = line.split()
                position = float(parts[3])  # Extract position (assuming it's the fourth column)
                positions.append(position)
                gmap = float(parts[2])
                gmaps.append(gmap)

    # Calculate combination rates and record them
    combined_rates = []
    for i in range(1, len(positions)):
        pos1 = positions[i - 1]
        pos2 = positions[i]
        gmap1 = gmaps[i - 1]
        gmap2 = gmaps[i]
        combination_rate = ((gmap2 - gmap1) / (pos2 - pos1)) * 1_000_000 * 1e-8
        combined_rates.append(combination_rate)
    
    rate_map = msprime.RateMap(
        position=positions,
        rate=combined_rates
    )

    return rate_map, positions, combined_rates


def export_output(mts, name, positions, combined_rates):
    a = arg_needle_lib.tskit_to_arg(mts)
    arg_needle_lib.serialize_arg(a, name + ".argn")

    # Output Tree Sequence to a file
    mts.dump(name + ".trees")

    # Output VCF (Variant Call Format) file
    vcf_file = name + ".vcf"
    with open(vcf_file, "w") as vcf_out:
        mts.write_vcf(vcf_out)
    
    with open(name + ".csv", "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
    
        # Write header
        writer.writerow(["Position", "Rate"])
    
        # Write data rows
        for i in range(len(positions) - 1):
            writer.writerow([positions[i], combined_rates[i]])


def sim_one_const(pop_size, sample_size, rate_map, mu_rate):
    # Create a demographic model for a constant-size population
    demographic_model = msprime.Demography()
    demographic_model.add_population(initial_size=pop_size)
    
    # Simulate the tree sequence
    ts = msprime.sim_ancestry(samples=sample_size, demography=demographic_model, recombination_rate=rate_map,
                              ploidy=2)
    
    # Simulate mutations
    mts = msprime.sim_mutations(ts, rate= mu_rate)
    
    return mts

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate trees')
    parser.add_argument("-s", '--sample', dest="sample_size", required=True, help="Number of haploid samples, which is not the same as the diploid sample size required by msprime")
    #parser.add_argument("-l", '--len', dest="len", required=True, help="Length of the sequence")
    parser.add_argument("-p", '--pop', dest="pop_size", required=True, help="Effective population size")
    parser.add_argument("-mu", '--mutation', dest="mu_rate", required=True, help="Mutation rate")
    parser.add_argument("-m", '--map', dest="map_file", required=True, help="Path to the genetic map file")
    parser.add_argument("-o", '--output', dest="output", required=True, help="Output file prefix")
    
    args = vars(parser.parse_args())
    print(args)
    
    pop_size = int(args['pop_size'])
    sample_size = int(args['sample_size'])
    #length = int(args['len'])
    mu_rate = float(args['mu_rate'])
    map_file = args['map_file']
    
    rate_map, positions, combined_rates = create_rate_map(map_file)
    mts = sim_one_const(pop_size, sample_size, rate_map, mu_rate)
    export_output(mts, args['output'],positions, combined_rates)

# bash script
#!/bin/bash
#SBATCH --job-name=ARG_job
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G


#source myenv/bin/activate
# python test.py -s 100 -p 1000 -mu 1e-8 -m plink.chr1.GRCh38.map -o output_test