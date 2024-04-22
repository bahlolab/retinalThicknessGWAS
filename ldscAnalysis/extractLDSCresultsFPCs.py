#!/usr/bin/env python

import os

# Function to extract h2 and se from log file
def extract_h2_from_log(log_file):
    with open(log_file, 'r') as f:
        lines = f.readlines()
    
    h2 = None
    se = None
    for line in lines:
        if "Total Observed scale h2:" in line:
            parts = line.split("(")
            h2 = float(parts[0].split(":")[1].strip())
            se = float(parts[1].split(")")[0].strip())
            break
    
    return h2, se

# Function to extract intercept and se from log file
def extract_intercept_from_log(log_file):
    with open(log_file, 'r') as f:
        lines = f.readlines()
    
    intercept = None
    se = None
    for line in lines:
        if "Intercept:" in line:
            parts = line.split("(")
            intercept = float(parts[0].split(":")[1].strip())
            se = float(parts[1].split(")")[0].strip())
            break
    
    return intercept, se

# univariate
# Output files
output_h2 = open("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_univariate_h2.txt", "w")
output_intercept = open("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_univariate_ldscIntercept.txt", "w")

# Write headers to output files
output_h2.write("fPC\th2\tse\n")
output_intercept.write("fPC\tintercept\tse\n")

# Iterate over phenotypes 1:6
for phenotype in  ["fPC1", "fPC2", "fPC3", "fPC4", "fPC5", "fPC6"]:
    
    # Extract h2 and se from log file
    log_file= f"/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs/{phenotype}.univariate_h2.log"
    h2, se_h2 = extract_h2_from_log(log_file)
    
    # Write h2 to output file
    output_h2.write(f"{phenotype}\t{h2}\t{se_h2}\n")
    
    # Extract intercept and se from log file
    intercept, se_intercept = extract_intercept_from_log(log_file)
    
    # Write intercept to output file
    output_intercept.write(f"{phenotype}\t{intercept}\t{se_intercept}\n")

# Close output files
output_h2.close()
output_intercept.close()

# same for stratified
# Output files
output_h2_stratified = open("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_stratified_h2.txt", "w")
output_intercept_stratified = open("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs_stratified_ldscIntercept.txt", "w")

# Write headers to output files
output_h2_stratified.write("fPC\th2\tse\n")
output_intercept_stratified.write("fPC\tintercept\tse\n")


for phenotype in  ["fPC1", "fPC2", "fPC3", "fPC4", "fPC5", "fPC6"]:

    # Extract h2 and se from log file
    log_file= f"/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs/{phenotype}.stratified_h2.log"
    h2, se_h2 = extract_h2_from_log(log_file)

    # Write h2 to output file
    output_h2_stratified.write(f"{phenotype}\t{h2}\t{se_h2}\n")

    # Extract intercept and se from log file
    intercept, se_intercept = extract_intercept_from_log(log_file)
    
    # Write intercept to output file
    output_intercept_stratified.write(f"{phenotype}\t{intercept}\t{se_intercept}\n")

# Close output files
output_h2_stratified.close()
output_intercept_stratified.close()

print("Files generated successfully!")



