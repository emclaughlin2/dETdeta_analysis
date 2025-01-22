import re
import numpy as np

# Initialize variables
energies = []

# Read the file
with open("run_data_log.txt", "r") as file:
    for line in file:
        # Check for lines with "high chi2 channel"
        if "OHCal high chi2 channel" in line:
            # Extract the energy value using a regular expression
            match = re.search(r"energy\s+([-+]?[0-9]*\.?[0-9]+)", line)
            if match:
                energies.append(float(match.group(1)))

# Calculate average and standard deviation
if energies:
    avg_energy = np.mean(energies)
    std_dev_energy = np.std(energies)
    print(len(energies))
    print(f"Average Energy: {avg_energy:.6f}")
    print(f"Standard Deviation of Energy: {std_dev_energy:.6f}")
else:
    print("No energy values found.")

