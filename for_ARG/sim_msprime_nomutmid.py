import msprime
import numpy as np
import pandas as pd
import argparse
import gzip
import arg_needle_lib
import csv
import random
#import timeancestry as tac

def create_rate_map(map_file):
    # Initialize lists to store positions and rates
    positions = [0]
    gmaps = [0]

    # Read positions from map.txt file
    with open(map_file, 'r') as f:
        for line in f:
            if line.strip():  # Check if line is not empty
                parts = line.split()
                position = float(parts[3])  # Extract position (assuming it's the third column)
                positions.append(position)
                gmap = float(parts[2])
                gmaps.append(gmap)

    ## random 500kb segement
    # Ensure there is enough data for a 500kb segment
    if positions[-1] - positions[0] < 500_000:
        raise ValueError("Not enough data for a 500kb segment.")

    # Randomly select a start position
    max_start_position = positions[-1] - 500_000
    start_position = random.uniform(positions[0], max_start_position)
    end_position = start_position + 500_000

    # Filter positions and rates to only include the middle 500kb
    filtered_positions = []
    filtered_gmaps = []
    for pos, gmap in zip(positions, gmaps):
        if start_position <= pos <= end_position:
            filtered_positions.append(pos)
            filtered_gmaps.append(gmap)

    # Normalize positions to ensure the first position is zero
    normalized_positions = [pos - filtered_positions[0] for pos in filtered_positions]

    # Calculate combination rates and record them
    combined_rates = []
    for i in range(1, len(normalized_positions)):
        pos1 = normalized_positions[i - 1]
        pos2 = normalized_positions[i]
        gmap1 = filtered_gmaps[i - 1]
        gmap2 = filtered_gmaps[i]
        combination_rate = ((gmap2 - gmap1) / (pos2 - pos1)) * 1_000_000 * 1e-8
        combined_rates.append(combination_rate)

    rate_map = msprime.RateMap(
        position=normalized_positions,
        rate=combined_rates
    )

    with open("recombination_rates_chr1mid.csv", "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
    
        # Write header
        writer.writerow(["Position", "Rate"])
    
        # Write data rows
        for i in range(len(normalized_positions) - 1):
            writer.writerow([normalized_positions[i], combined_rates[i]])
    
    return rate_map



def export_output(ts, name):
    a = arg_needle_lib.tskit_to_arg(ts)
    arg_needle_lib.serialize_arg(a, name + ".argn")

    # Output Tree Sequence to a file
    ts.dump(name + ".trees")

    # Output VCF (Variant Call Format) file
    vcf_file = name + ".vcf"
    with open(vcf_file, "w") as vcf_out:
        ts.write_vcf(vcf_out)

def sim_one_const(pop_size, sample_size, rate_map):
    # Create a demographic model for a constant-size population
    demographic_model = msprime.Demography()
    demographic_model.add_population(initial_size=pop_size)
    
    # Simulate the tree sequence
    ts = msprime.sim_ancestry(samples=sample_size, demography=demographic_model, recombination_rate=rate_map,
                              ploidy=2)
    
    # Simulate mutations
    #mts = msprime.sim_mutations(ts, rate= mu_rate)
    
    return ts

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate trees')
    parser.add_argument("-s", '--sample', dest="sample_size", required=True, help="Number of haploid samples, which is not the same as the diploid sample size required by msprime")
    #parser.add_argument("-l", '--len', dest="len", required=True, help="Length of the sequence")
    parser.add_argument("-p", '--pop', dest="pop_size", required=True, help="Effective population size")
    #parser.add_argument("-mu", '--mutation', dest="mu_rate", required=True, help="Mutation rate")
    parser.add_argument("-m", '--map', dest="map_file", required=True, help="Path to the genetic map file")
    parser.add_argument("-o", '--output', dest="output", required=True, help="Output file prefix")
    
    args = vars(parser.parse_args())
    print(args)
    
    pop_size = int(args['pop_size'])
    sample_size = int(args['sample_size'])
    #length = int(args['len'])
    #mu_rate = float(args['mu_rate'])
    map_file = args['map_file']
    
    rate_map = create_rate_map(map_file)
    mts = sim_one_const(pop_size, sample_size, rate_map)
    export_output(mts, args['output'])

