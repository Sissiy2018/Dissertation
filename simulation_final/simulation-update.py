import msprime
import numpy as np
import pandas as pd
import argparse
import gzip
import csv
import random
import arg_needle_lib

def convert_map_to_hapmap(input_map_file, output_hapmap_file):
    with open(input_map_file, 'r') as infile, open(output_hapmap_file, 'w') as outfile:
        # Write the header
        outfile.write("Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n")

        # Read and process each line
        for line in infile:
            parts = line.strip().split()
            if len(parts) == 4:
                chromosome = parts[0]
                position_bp = parts[3]
                map_cM = parts[2]
                rate_cM_per_Mb = "0.000000"  # Assuming a constant rate between positions
                outfile.write(f"{chromosome}\t{position_bp}\t{rate_cM_per_Mb}\t{map_cM}\n")

def create_rate_map20Mb(map_file, convert_map_file, seed=None):
    if seed is not None:
        random.seed(seed)
    
    # Read the HapMap format recombination map
    rate_map = msprime.RateMap.read_hapmap(convert_map_file,rate_col=None, )
    
    print(rate_map.right)
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

    with open("20Mbrates_chr22_250_Ne5000.csv", "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(["Position", "Rate"])
        # Write data rows
        for pos, rate in zip(sliced_rate_map.left, sliced_rate_map.rate):
            writer.writerow([pos, rate])
    
    # Calculate genetic positions based on rates
    genetic_positions = [0]
    cumulative_distance = 0
    for i in range(0, len(sliced_rate_map.left)):
        physical_distance = sliced_rate_map.span[i]
        rate = sliced_rate_map.rate[i]
        genetic_distance = rate * physical_distance*100
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
        
        # Write the filtered rows to the output file
    with open('20Mbmap_chr22_250_Ne5000.map', 'w') as outfile:
        outfile.writelines(filtered_rows) 
    
    # Write the map file for thread
    with open("20Mbmap_chr22_250_Ne5000_forthread.map", 'w') as mapfile:
        for pos, gmap in zip(sliced_rate_map.left, genetic_positions):
            chromosome = 22
            snp_id = f"."
            genetic_distance = gmap
            physical_position = pos
            mapfile.write(f"{chromosome}\t{snp_id}\t{genetic_distance}\t{physical_position}\n")
    
    return sliced_rate_map

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
    pop_size = 5000
    sample_size = 250
    mu_rate = 1e-8  # Mutation rate per base pair per generation
    map_file = 'plink.chr22.GRCh38.map'
    convert_map_file = 'chr22_transform.map'
    output_prefix = 'chr22_250_Ne5000'

    convert_map_file = convert_map_to_hapmap(map_file, convert_map_file)

    print("Creating rate map...")
    rate_map = create_rate_map20Mb(map_file, 'chr22_transform.map', seed)
    
    print("Running the simulation...")
    mts = sim_one_const(pop_size, sample_size, rate_map, mu_rate, seed)
    
    print("Exporting VCF...")
    export_vcf(mts, output_prefix)
    print("Done!")
