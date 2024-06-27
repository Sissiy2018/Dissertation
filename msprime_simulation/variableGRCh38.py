map_file = "plink.chr1.GRCh38.map"
output_file = "combined_chr1.txt"

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

# Calculate combination rates and record them
combined_rates = []
for i in range(1, len(positions)):
    pos1 = positions[i - 1]
    pos2 = positions[i]
    gmap1 = gmaps[i - 1]
    gmap2 = gmaps[i]
    combination_rate = ((gmap2 - gmap1)/(pos2-pos1))*(10^6)
    combined_rates.append(combination_rate)

import msprime
from IPython.display import SVG, display
import sys
import tskit

rate_map = msprime.RateMap(
    position= positions,
    rate= combined_rates
)

pop_size = 20000
num_dip = 1000

# Simulate ancestry with a constant population size
ts = msprime.sim_ancestry(samples=num_dip, population_size=pop_size,
                          random_seed=42, recombination_rate=rate_map,
                          model=msprime.StandardCoalescent())
ts

ts.dump("variableGRCh38_20000pop_1000dip.trees")

vcf_file = "variableGRCh38_20000pop_1000dip.vcf"
with open(vcf_file, "w") as vcf_out:
    ts.write_vcf(vcf_out)