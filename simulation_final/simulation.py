import msprime
import numpy as np
import pandas as pd
import argparse
import gzip
import csv
import random
import arg_needle_lib

def create_rate_map10Mb(map_file, seed=None):
    if seed is not None:
        random.seed(seed)
    
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
    
    # Random 10Mb segment
    # Ensure there is enough data for a segment
    if positions[-1] - positions[0] < 10_000_000:
        raise ValueError("Not enough data for a 10Mb segment.")

    # Randomly select a start position
    max_start_position = positions[-1] - 10_000_000
    start_position = random.uniform(positions[0], max_start_position)
    end_position = start_position + 10_000_000

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
        combination_rate = ((gmap2 - gmap1) / (pos2 - pos1))
        combined_rates.append(combination_rate)

    rate_map = msprime.RateMap(
        position=normalized_positions,
        rate=combined_rates
    )

    with open("10Mbrates_chr22_500_Ne10000.csv", "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
    
        # Write header
        writer.writerow(["Position", "Rate"])
    
        # Write data rows
        for i in range(len(normalized_positions) - 1):
            writer.writerow([normalized_positions[i], combined_rates[i]])
    
    # Calculate genetic positions based on rates
    genetic_positions = [0]  # The first genetic position is set to 0
    cumulative_distance = 0
    for i in range(1, len(normalized_positions)):
        physical_distance = normalized_positions[i] - normalized_positions[i - 1]
        rate = combined_rates[i - 1]
        genetic_distance = rate * physical_distance
        cumulative_distance += genetic_distance
        genetic_positions.append(round(cumulative_distance, 6))
    
    with open("10Mbmap_chr22_500_Ne10000.map", "w") as mapfile:
        for index, (pos, gmap) in enumerate(zip(filtered_positions, filtered_gmaps)):
            chromosome = 22
            snp_id = f"."
            genetic_distance = gmap
            physical_position = pos
            mapfile.write(f"{chromosome}\t{snp_id}\t{genetic_distance}\t{physical_position}\n")
    
    # Write the map file
    with open("10Mbmap_chr22_500_Ne10000_forthread.map", 'w') as mapfile:
        for pos, gmap in zip(normalized_positions, genetic_positions):
            chromosome = 22
            snp_id = f"."
            genetic_distance = gmap
            physical_position = pos
            mapfile.write(f"{chromosome}\t{snp_id}\t{genetic_distance}\t{physical_position}\n")
    
    return rate_map

def export_vcf(mts, name):
    # Output VCF (Variant Call Format) file
    vcf_file = name + ".vcf"
    with open(vcf_file, "w") as vcf_out:
        mts.write_vcf(vcf_out)

    arg = arg_needle_lib.tskit_to_arg(mts)
    arg_needle_lib.serialize_arg(arg, name + ".argn")


def sim_one_const(pop_size, sample_size, rate_map, mu_rate, seed=None):
    # Create a demographic model for a constant-size population
    demographic_model = msprime.Demography()
    demographic_model.add_population(initial_size=pop_size)
    
    # Simulate the tree sequence
    ts = msprime.sim_ancestry(
        samples=sample_size, 
        demography=demographic_model, 
        recombination_rate=rate_map,
        ploidy=2,
        random_seed=seed
    )
    
    # Simulate mutations
    mts = msprime.sim_mutations(ts, rate=mu_rate, random_seed=seed)
    
    return mts

if __name__ == '__main__':
    seed = 42
    pop_size = 10000
    sample_size = 500
    mu_rate = 1e-8  # Mutation rate per base pair per generation
    map_file = 'plink.chr22.GRCh38.map'
    output_prefix = 'chr22_10Mb_500_Ne10000'

    print("Creating rate map...")
    rate_map = create_rate_map10Mb(map_file, seed)
    
    print("Running the simulation...")
    mts = sim_one_const(pop_size, sample_size, rate_map, mu_rate, seed)
    
    print("Exporting VCF...")
    export_vcf(mts, output_prefix)
    print("Done!")
