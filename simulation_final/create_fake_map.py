import csv

def create_map_from_csv(csv_file, map_file):
    # Initialize lists to store positions and rates
    positions = []
    rates = []

    # Read positions and rates from the CSV file
    with open(csv_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header
        for row in reader:
            position = float(row[0])
            rate = float(row[1])
            positions.append(position)
            rates.append(1e-8)
    
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

if __name__ == "__main__":
    csv_file = "manualrates_20Mb_500_Ne10000.csv"
    map_file = "fake20Mb_500_Ne10000.map"
    create_map_from_csv(csv_file, map_file)