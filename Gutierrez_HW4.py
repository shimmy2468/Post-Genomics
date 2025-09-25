import pandas as pd

# Read the CSV file
df = pd.read_csv("Final_Tumor.csv")

# Extract required columns and create lines for both strands (1 and -1)
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

    # Create one line for the positive strand (1) and one for the negative strand (-1)
    line_pos = f"chromosome {chromosome} {position} {ref_allele} {alt_allele} 1"
    line_neg = f"chromosome {chromosome} {position} {ref_allele} {alt_allele} -1"
    results.append(line_pos)
    results.append(line_neg)

# Write to output file
with open("Final_Tumor_forNexus.txt", 'w') as f:
    for line in results:
        f.write(line + '\n')


