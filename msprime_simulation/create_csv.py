import msprime
import numpy as np
import pandas as pd
import argparse
import gzip
import csv
import random
import os
import arg_needle_lib

def create_recombination_rate_csv(output_csv, map_file, step=2000):
    # Define the regions and gaps
    regions = [
        (0, 3_000_000, 1e-7),
        (3_500_000, 5_500_000, 0.5e-7),
        (6_000_000, 8_000_000, 1e-8),
        (8_500_000, 10_500_000, 0.7e-8),
        (11_000_000, 13_000_000, 0.4e-8),
        (13_500_000, 15_500_000, 1e-9),
        (16_000_000, 18_000_000, 0.7e-9),
        (18_500_000, 20_000_000, 0.4e-9),
    ]
    gaps = [(3_000_000, 3_500_000), (5_500_000, 6_000_000), (8_000_000, 8_500_000), (10_500_000, 1_100_000),
            (13_000_000, 13_500_000), (15_500_000, 16_000_000), (18_000_000, 18_500_000), ]

    positions_and_rates = []

    # Write regions with constant recombination rates and noise
    for start, end, rate in regions:
        current_position = start
        while current_position < end:
            # Introduce noise: ±5% of the rate
            noise = rate * random.uniform(-0.05, 0.05)
            noisy_rate = rate + noise
            positions_and_rates.append((current_position, noisy_rate))
            current_position += step  # Increment by step size

    # Write gaps with zero recombination rate
    for start, end in gaps:
        current_position = start
        while current_position < end:
            positions_and_rates.append((current_position, 0))
            current_position += step  # Increment by step size

    # Sort positions and rates by position
    positions_and_rates.sort()

    # Write sorted positions and rates to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Position", "Rate"])
        writer.writerows(positions_and_rates)

    positions, rates = zip(*positions_and_rates)

    rate_map = msprime.RateMap(
        position=positions,
        rate=rates[:-1],
    )

    # Calculate genetic positions based on rates
    genetic_positions = [0]  # The first genetic position is set to 0
    cumulative_distance = 0
    for i in range(1, len(positions)):
        physical_distance = positions[i] - positions[i - 1]
        rate = rates[i - 1]
        genetic_distance = rate * physical_distance
        cumulative_distance += genetic_distance
        genetic_positions.append(round(cumulative_distance,6))
    
    # Write the map file
    with open(map_file, 'w') as mapfile:
        for pos, gmap in zip(positions, genetic_positions):
            chromosome = 22
            snp_id = "."
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
    pop_size = 10000
    sample_size = 500
    mu_rate = 1e-8  # Mutation rate per base pair per generation
    map_file = 'test_manual.map'
    output_prefix = 'test'

    # Define the output paths
    recombination_rate_csv = "test_manual.csv"

    # Create recombination rate CSV and genetic map file
    rate_map = create_recombination_rate_csv(recombination_rate_csv, map_file)
    
    # Run the simulation
    mts = sim_one_const(pop_size, sample_size, rate_map, mu_rate)
    
    # Export VCF
    export_vcf(mts, output_prefix)