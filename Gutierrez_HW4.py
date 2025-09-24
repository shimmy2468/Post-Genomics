import pandas as pd

# Read the CSV file
df = pd.read_csv("Final_Normal.csv")
    
# Extract required columns
results = []
for _, row in df.iterrows():
    # Extract chromosome (remove 'chr' prefix)
    chromosome = str(row.iloc[2]).replace('chr', '')
    # Extract position left coordinate
    position = row.iloc[3]
    # Extract reference and alternative alleles
    ref_allele = row.iloc[5]  # ref_seq
    if row.iloc[6] == row.iloc[5]:  # because there are two var_seq columns
        alt_allele = row.iloc[7]  # var_seq2
    else:
        alt_allele = row.iloc[6]  # var_seq1
    
    # Format: chromosome position ref_allele alt_allele 1
    snp_line = f"chromosome {chromosome} {position} {ref_allele} {alt_allele} 1"
    results.append(snp_line)

# Write to output file
with open("Final_Normal_forNexus.txt", 'w') as f:
    for line in results:
        f.write(line + '\n')


