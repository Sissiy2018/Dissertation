import msprime
import numpy as np
import pandas as pd
import argparse
import gzip
import csv
import random
import arg_needle_lib

def create_recombination_rate_csv(output_prefix, stage_number, step=2000, seed = None):
    if seed is not None:
        random.seed(seed)
    
    # Define the regions and gaps
    regions = [
        (0, 20_000_000, 2e-08),
    ]
    gaps = [ ]

    positions_and_rates = []

    # Write regions with constant recombination rates and noise
    for start, end, rate in regions:
        current_position = start
        while current_position <= end:
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
    csv_filename = f"{output_prefix}_stage{stage_number}.csv"
    with open(csv_filename, 'w', newline='') as csvfile:
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
        genetic_distance = rate * physical_distance * 100
        cumulative_distance += genetic_distance
        genetic_positions.append(round(cumulative_distance,6))
    
    # Write the map file
    map_filename = f'{output_prefix}_stage{stage_number}.map'
    with open(map_filename, 'w') as mapfile:
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

def sim_epoch(pop_size, sample_size, mu_rate, start_time, end_time, initial_state=None, rate_map=None, seed=None):
    # Create a demographic model for a constant-size population
    demographic_model = msprime.Demography()
    demographic_model.add_population(initial_size=pop_size)
    
    if initial_state is not None:
        # Continue the simulation from the provided initial state
        ts = msprime.sim_ancestry(
            sequence_length = 20_000_000,
            initial_state=initial_state,
            demography=demographic_model, 
            #recombination_rate=2e-08,
            recombination_rate = rate_map,
            ploidy=2,
            random_seed=seed,
            start_time=start_time,
            end_time=end_time
        )
    else:
        # Simulate the tree sequence for the epoch from scratch
        ts = msprime.sim_ancestry(
            sequence_length = 20_000_000,
            samples=sample_size, 
            demography=demographic_model, 
            #recombination_rate=2e-08,
            recombination_rate = rate_map,
            ploidy=2,
            random_seed=seed,
            start_time=start_time,
            end_time=end_time
        )
    
    # Simulate mutations
    mts = msprime.sim_mutations(ts, rate=mu_rate, start_time=start_time,
            end_time=end_time, random_seed=seed)

    max_time = max(mts.tables.nodes.time)
    total_generations = int(max_time)

    print(f"The total number of generations is {total_generations}")
    
    return mts

if __name__ == '__main__':
    seed = 123
    pop_size = 10000
    sample_size = 500
    mu_rate = 1e-8  # Mutation rate per base pair per generation
    map_file = 'plink.chr22.GRCh38.map'
    convert_map_file = 'chr22_transform.map'
    output_prefix = f'gap_constant_500_Ne10000'

    print("Creating rate map for first epoch...")
    rate_map1 = create_recombination_rate_csv(output_prefix, 1, seed=seed)
    
    print("Running the simulation for generations 1 to 5000...")
    mts1 = sim_epoch(pop_size, sample_size, mu_rate, start_time=0, end_time=5000, 
                     initial_state=None, rate_map=None, seed=seed)

    print("Creating rate map for second epoch...")
    rate_map2 = create_recombination_rate_csv(output_prefix, 2, seed=seed + 1)
    
    print("Running the simulation for generations 5000 to 10000...")
    mts2 = sim_epoch(pop_size, sample_size, mu_rate, start_time=5000, end_time=10000, 
                     initial_state=mts1, rate_map = None, seed=seed+1)

    print("Creating rate map for third epoch...")
    rate_map3 = create_recombination_rate_csv(output_prefix, 3, seed=seed + 2)

    print("Running the simulation for generations 10000 to ...")
    mts3 = sim_epoch(pop_size, sample_size, mu_rate, start_time=10000, end_time=None, 
                     initial_state=mts2, rate_map= rate_map3, seed=seed+2)
    
    print("Exporting VCF...")
    export_vcf(mts3, output_prefix)
    print("Done!")