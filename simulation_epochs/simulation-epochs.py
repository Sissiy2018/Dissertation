import msprime
import numpy as np
import pandas as pd
import argparse
import gzip
import csv
import random
import arg_needle_lib

def create_rate_map20Mb(map_file, convert_map_file, output_prefix, stage_number, seed=None):
    if seed is not None:
        random.seed(seed)
    
    # Read the HapMap format recombination map
    rate_map = msprime.RateMap.read_hapmap(convert_map_file,rate_col=None, )
    
    # Randomly select a 10Mb segment
    total_length = rate_map.right[-1] - rate_map.right[0]
    if total_length < 20_000_000:
        raise ValueError("Not enough data for a 20Mb segment.")
    
    #print(total_length)
    max_start_position = rate_map.right[-1] - 20_000_000
    
    #print(max_start_position)
    start_position = random.choice([x for x in rate_map.left if rate_map.right[0] < x < max_start_position])
    end_position = max((x for x in rate_map.right if x < start_position + 20_000_000), default=None)

    # Slice the map to get the 20Mb segment
    sliced_rate_map = rate_map.slice(start_position, end_position, trim = True)

    new_left = [sliced_rate_map.right[-1], 20_000_000]
    new_position = np.append(sliced_rate_map.left, new_left)

    new_rate = np.append(sliced_rate_map.rate, 2e-08) 

    new_rate_map = msprime.RateMap(position = new_position, rate = new_rate)

    csv_filename = f"{output_prefix}_stage{stage_number}.csv"
    with open(csv_filename, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(["Position", "Rate"])
        # Write data rows
        for pos, rate in zip(new_rate_map.left, new_rate_map.rate):
            writer.writerow([pos, rate])
    
    # Calculate genetic positions based on rates
    genetic_positions = [0]
    cumulative_distance = 0
    for i in range(0, len(new_rate_map.left)):
        physical_distance = new_rate_map.span[i]
        rate = new_rate_map.rate[i]
        genetic_distance = rate * physical_distance * 100
        cumulative_distance += genetic_distance
        genetic_positions.append(round(cumulative_distance, 6))
    
    filtered_rows = []
    with open(map_file, 'r') as infile:
        # Read the file line by line
        for line in infile:
            # Split the line into columns based on tab or space
            columns = line.split()
        
            if len(columns) >= 4:
                value = float(columns[3])
                        
                # Check if the value falls within the specified range
                if start_position <= value <= end_position:
                    filtered_rows.append(line)
    
    map_filename = f'{output_prefix}_stage{stage_number}.map'
    # Write the filtered rows to the output file
    with open(map_filename, 'w') as outfile:
        outfile.writelines(filtered_rows) 
    
    thread_map_filename = f'{output_prefix}_stage{stage_number}_forthread.map'
    # Write the map file for thread
    with open(thread_map_filename, 'w') as mapfile:
        for pos, gmap in zip(new_rate_map.left, genetic_positions):
            chromosome = 22
            snp_id = f"."
            genetic_distance = gmap
            physical_position = pos
            mapfile.write(f"{chromosome}\t{snp_id}\t{genetic_distance}\t{physical_position}\n")
    
    return new_rate_map

def export_vcf(mts, name):
    # Output VCF (Variant Call Format) file
    vcf_file = name + ".vcf"
    with open(vcf_file, "w") as vcf_out:
        mts.write_vcf(vcf_out)

    arg = arg_needle_lib.tskit_to_arg(mts)
    arg_needle_lib.serialize_arg(arg, name + ".argn")

def sim_epoch(pop_size, sample_size, rate_map, mu_rate, start_time, end_time, initial_state=None, seed=None):
    # Create a demographic model for a constant-size population
    demographic_model = msprime.Demography()
    demographic_model.add_population(initial_size=pop_size)
    
    if initial_state is not None:
        # Continue the simulation from the provided initial state
        ts = msprime.sim_ancestry(
            #sequence_length = 20_000_000,
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
            #sequence_length = 20_000_000,
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
    mts = msprime.sim_mutations(ts, rate=mu_rate, random_seed=seed)
    
    return mts

if __name__ == '__main__':
    seed = 123
    pop_size = 10000
    sample_size = 500
    mu_rate = 1e-8  # Mutation rate per base pair per generation
    map_file = 'plink.chr22.GRCh38.map'
    convert_map_file = 'chr22_transform.map'
    output_prefix = f'constant_chr22_500_Ne10000'

    print("Creating rate map for first epoch...")
    rate_map1 = create_rate_map20Mb(map_file, convert_map_file, output_prefix, 1, seed=seed)
    
    print("Running the simulation for generations 1 to 5000...")
    mts1 = sim_epoch(pop_size, sample_size, rate_map1, mu_rate, start_time=0, end_time=5000, initial_state=None, seed=seed)

    print("Creating rate map for second epoch...")
    rate_map2 = create_rate_map20Mb(map_file, convert_map_file, output_prefix, 2, seed=seed+1)
    
    print("Running the simulation for generations 5000 to 10000...")
    mts2 = sim_epoch(pop_size, sample_size, rate_map2, mu_rate, start_time=5000, end_time=10000, initial_state=mts1, seed=seed+1)

    print("Creating rate map for third epoch...")
    rate_map3 = create_rate_map20Mb(map_file, convert_map_file, output_prefix, 3, seed=seed+2)

    print("Running the simulation for generations 10000 to ...")
    mts3 = sim_epoch(pop_size, sample_size, rate_map3, mu_rate, start_time=10000, end_time=None, initial_state=mts2, seed=seed+2)
    
    print("Exporting VCF...")
    export_vcf(mts3, output_prefix)
    print("Done!")