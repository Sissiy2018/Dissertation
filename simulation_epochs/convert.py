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

if __name__ == '__main__':
    
    map_file = 'plink.chr22.GRCh38.map'
    convert_map_file = 'chr22_transform.map'

    convert_map_file = convert_map_to_hapmap(map_file, convert_map_file)