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
        (0, 3_000_000, 0.8e-8),
        (3_500_001, 5_500_000, 0.6e-8),
        (6_000_001, 8_000_000, 0.4e-8),
        (8_500_001, 10_500_000, 0.2e-8),
        (11_000_001, 13_000_000, 1e-9),
        (13_500_001, 15_500_000, 0.8e-9),
        (16_000_001, 18_000_000, 0.6e-9),
        (18_500_001, 20_000_000, 0.4e-9),
    ]
    gaps = [(3_000_001, 3_500_000), (5_500_001, 6_000_000), (8_000_001, 8_500_000), (10_500_001, 1_100_000),
            (13_000_001, 13_500_000), (15_500_001, 16_000_000), (18_000_001, 18_500_000), ]

    # Initialize arrays to store x and y values
    x = []
    y = []

# Generate y values based on normal distributions in each region
    for start, end, peak in regions:
        mu = peak
        sigma = peak * 0.5
        y_values = np.linspace(0, 2*mu, 1000)
        x_values = np.linspace(start, end, 1000)
        pdf_values = (1/(sigma * np.sqrt(2 * np.pi))) * np.exp(-((y_values - mu)**2) / (2 * sigma**2))/1e16
        x.extend(x_values.tolist())
        y.extend(pdf_values.tolist())

    # Set y values to zero in the gaps
    for gap_start, gap_end in gaps:
        x_values = np.linspace(gap_start, gap_end, 1000)
        y_values = np.zeros_like(x_values)
        x.extend(x_values.tolist())
        y.extend(y_values.tolist())

    # Combine x and y into positions_and_rates list
    positions_and_rates = [(x[i], y[i]) for i in range(len(x))]

    # Sort positions_and_rates by position (x values)
    positions_and_rates.sort()

    positions, rates = zip(*positions_and_rates)

    # Write sorted positions and rates to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Position", "Rate"])
        writer.writerows(positions_and_rates)

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
    map_file = 'normal_manual.map'
    output_prefix = 'normal'

    # Define the output paths
    recombination_rate_csv = "normal_manual.csv"

    # Create recombination rate CSV and genetic map file
    rate_map = create_recombination_rate_csv(recombination_rate_csv, map_file)
    
    # Run the simulation
    mts = sim_one_const(pop_size, sample_size, rate_map, mu_rate)
    
    # Export VCF
    export_vcf(mts, output_prefix)