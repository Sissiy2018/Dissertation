import msprime
import numpy as np
import pandas as pd
import argparse
import gzip
import csv
import random
import arg_needle_lib

def create_rate_map(map_file, convert_map_file, output_prefix):
    
    # Read the HapMap format recombination map
    rate_map = msprime.RateMap.read_hapmap(convert_map_file,rate_col=None, )

    csv_filename = f"{output_prefix}.csv"
    with open(csv_filename, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(["Position", "Rate"])
        # Write data rows
        for pos, rate in zip(rate_map.left, rate_map.rate):
            writer.writerow([pos, rate])
    
    # Calculate genetic positions based on rates
    genetic_positions = [0]
    cumulative_distance = 0
    for i in range(0, len(rate_map.left)):
        physical_distance = rate_map.span[i]
        rate = rate_map.rate[i]
        genetic_distance = rate * physical_distance * 100
        cumulative_distance += genetic_distance
        genetic_positions.append(round(cumulative_distance, 6))
    
    thread_map_filename = f'{output_prefix}_forthread.map'
    # Write the map file for thread
    with open(thread_map_filename, 'w') as mapfile:
        for pos, gmap in zip(rate_map.left, genetic_positions):
            chromosome = 22
            snp_id = f"."
            genetic_distance = gmap
            physical_position = pos
            mapfile.write(f"{chromosome}\t{snp_id}\t{genetic_distance}\t{physical_position}\n")
    
    return rate_map

if __name__ == '__main__':
    map_file = 'plink.chr22.GRCh38.map'
    convert_map_file = 'chr22_transform.map'
    output_prefix = f'chr22_500_Ne10000'
    
    rate_map1 = create_rate_map(map_file, convert_map_file, output_prefix,)