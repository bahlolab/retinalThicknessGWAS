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

# Read phenotypes and slice from pixels.txt
phenotypes = []
slices = []
with open('pixels.txt', 'r') as f:
    for line in f:
        parts = line.strip().split("\t")
        phenotypes.append(parts[0])
        slices.append(parts[1])
        
# Output files
output_h2 = open("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/univariate_h2.txt", "w")
output_intercept = open("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/univariate_ldscIntercept.txt", "w")
output_h2_stratified = open("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/stratified_h2.txt", "w")
output_intercept_stratified = open("/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/stratified_ldscIntercept.txt", "w")

# Write headers to output files
output_h2.write("pixel\th2\tse\n")
output_intercept.write("pixel\tintercept\tse\n")
output_h2_stratified.write("pixel\th2\tse\n")
output_intercept_stratified.write("pixel\tintercept\tse\n")

# Iterate over phenotypes
for phenotype, slice in zip(phenotypes, slices):
    
    # Extract h2 and se from log file
    log_file= f"/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/{slice}/{phenotype}.univariate_h2.log"
    h2, se_h2 = extract_h2_from_log(log_file)
    
    # Write h2 to output file
    output_h2.write(f"{phenotype}\t{h2}\t{se_h2}\n")
    
    # Extract intercept and se from log file
    intercept, se_intercept = extract_intercept_from_log(log_file)
    
    # Write intercept to output file
    output_intercept.write(f"{phenotype}\t{intercept}\t{se_intercept}\n")

    ## same for stratified
    # Extract h2 and se from log file
    log_file= f"/vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/{slice}/{phenotype}.stratified_h2.log"
    h2, se_h2 = extract_h2_from_log(log_file)

    # Write h2 to output file
    output_h2_stratified.write(f"{phenotype}\t{h2}\t{se_h2}\n")
    
    # Extract intercept and se from log file
    intercept, se_intercept = extract_intercept_from_log(log_file)

    # Write intercept to output file
    output_intercept_stratified.write(f"{phenotype}\t{intercept}\t{se_intercept}\n")

# Close output files
output_h2.close()
output_intercept.close()
output_h2_stratified.close()
output_intercept_stratified.close()

print("Files generated successfully!")



