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

    # Create one line for the positive strand (1)
    line_pos = f"chromosome {chromosome} {position} {ref_allele} {alt_allele} 1"
    results.append(line_pos)

# Write to output file
with open("Final_Tumor_forNexus.txt", 'w') as f:
    for line in results:
        f.write(line + '\n')

# For clean merging of the original and pervariant files

def normalize_chrom(x):
    if pd.isna(x):
        return x
    x = str(x)
    # make lowercase 'chr' prefix consistent
    if x.startswith('chr'):
        return x
    return 'chr' + x


final_path = 'Final_Tumor.csv'
pervariant_path = 'tumor_pervariant.tsv'
fathmm_path = 'tumor_predictions.txt'

final = pd.read_csv(final_path, dtype=str, low_memory=False)
pervar = pd.read_csv(pervariant_path, sep='\t', dtype=str, low_memory=False)

# preserve original order
final['_orig_index'] = range(len(final))


final['chrom_norm'] = final['chrom'].astype(str).apply(normalize_chrom)

pervar['chrom_norm'] = pervar['chrom'].astype(str).apply(normalize_chrom)

if 'id' in pervar.columns:
        pervar['chrom_norm'] = pervar['id'].astype(str).str.split(':').str[0]

# Position columns
if 'left' in final.columns:
    final['pos_norm'] = final['left'].astype(int)
elif 'position' in final.columns:
    final['pos_norm'] = final['position'].astype(int)

if 'position' in pervar.columns:
    pervar['pos_norm'] = pervar['position'].astype(int)
else:
    # try to extract from id
    pervar['pos_norm'] = pervar['id'].astype(str).str.split(':').str[1].astype(int)

# Alleles: normalize ref/alt
if 'ref_seq' in final.columns:
    final['ref_norm'] = final['ref_seq'].astype(str)
elif 'ref' in final.columns:
    final['ref_norm'] = final['ref'].astype(str)
else:
    final['ref_norm'] = None

# pick var_seq1 if available for alt
if 'var_seq1' in final.columns:
    final['alt_norm'] = final['var_seq1'].astype(str)
elif 'Allele' in final.columns:
    # Allele might be single-letter alt
    final['alt_norm'] = final['Allele'].astype(str)
else:
    final['alt_norm'] = None

if 'ref' in pervar.columns:
    pervar['ref_norm'] = pervar['ref'].astype(str)
else:
    pervar['ref_norm'] = None

if 'alt' in pervar.columns:
    pervar['alt_norm'] = pervar['alt'].astype(str)
else:
    # try to parse from id (chr:pos:ref/alt:...)
    pervar['alt_norm'] = pervar['id'].astype(str).str.split(':').str[2].str.split('/').str[1]

# Decide merge keys - default to chrom+pos+ref+alt when alleles present
keys_final = ['chrom_norm', 'pos_norm']
keys_pervar = ['chrom_norm', 'pos_norm']

if final['ref_norm'].notna().any() and pervar['ref_norm'].notna().any():
    final['ref_norm'] = final['ref_norm'].replace({'None':'', 'nan':''}).fillna('')
    pervar['ref_norm'] = pervar['ref_norm'].replace({'None':'', 'nan':''}).fillna('')
    if final['ref_norm'].str.len().sum() > 0 and pervar['ref_norm'].str.len().sum() > 0:
        keys_final.append('ref_norm')
        keys_pervar.append('ref_norm')

if final['alt_norm'].notna().any() and pervar['alt_norm'].notna().any():
    final['alt_norm'] = final['alt_norm'].replace({'None':'', 'nan':''}).fillna('')
    pervar['alt_norm'] = pervar['alt_norm'].replace({'None':'', 'nan':''}).fillna('')
    if final['alt_norm'].str.len().sum() > 0 and pervar['alt_norm'].str.len().sum() > 0:
        keys_final.append('alt_norm')
        keys_pervar.append('alt_norm')

# Most of the code above is to make sure I get the right keys for merging with the data being so inconsistent

print('Merge keys chosen:', keys_final)

# Select columns from pervar to merge - exclude columns used for joining
pervar_cols = [c for c in pervar.columns if c not in (keys_pervar + ['_orig_index'])]

merged = final.merge(pervar[pervar_cols + keys_pervar], left_on=keys_final, right_on=keys_pervar, how='left')


out_path = 'Final_Tumor_merged.csv'
merged.to_csv(out_path, index=False)

# Quick verification stats
matches = merged[keys_final].notna().all(axis=1).sum()
total = len(merged)
print(f'Matched rows (non-null join keys): approx {matches} / {total}')

# Most of the code references 'norm', but I just change the paths for the normal and tumor files respectively with the script remaining the same. I just run it twice

# Last step with FATHMM-XF predictions
print('\nProcessing FATHMM-XF predictions...')

# Read FATHMM file - skip comment lines starting with #
fathmm = pd.read_csv(fathmm_path, sep='\t', comment='#', dtype=str, low_memory=False)

# Rename columns to avoid confusion 
fathmm.columns = ['fathmm_chrom', 'fathmm_position', 'fathmm_ref', 'fathmm_alt', 
                   'fathmm_coding_score', 'fathmm_noncoding_score', 'fathmm_warning']

# Normalize FATHMM chromosome (add 'chr' prefix if not present)
fathmm['chrom_norm'] = fathmm['fathmm_chrom'].astype(str).apply(lambda x: 'chr' + x if not x.startswith('chr') else x)

# Convert position to integer
fathmm['pos_norm'] = fathmm['fathmm_position'].astype(int)

# Normalize ref and alt alleles
fathmm['ref_norm'] = fathmm['fathmm_ref']
fathmm['alt_norm'] = fathmm['fathmm_alt']

# Also uppercase normalize the merged dataframe's ref/alt for matching
merged['ref_norm'] = merged['ref_norm'].astype(str).str.upper()
merged['alt_norm'] = merged['alt_norm'].astype(str).str.upper()

# Determine merge keys for FATHMM
keys_merged_fathmm = ['chrom_norm', 'pos_norm']
keys_fathmm = ['chrom_norm', 'pos_norm']

# Add ref and alt to merge keys if they exist
if 'ref_norm' in keys_final and 'alt_norm' in keys_final:
    keys_merged_fathmm.extend(['ref_norm', 'alt_norm'])
    keys_fathmm.extend(['ref_norm', 'alt_norm'])

print('Merge keys chosen for FATHMM:', keys_merged_fathmm)

# Select FATHMM columns to merge 
fathmm_cols_to_merge = ['fathmm_coding_score', 'fathmm_noncoding_score', 'fathmm_warning']

# Perform the merge with FATHMM predictions
merged = merged.merge(fathmm[fathmm_cols_to_merge + keys_fathmm], 
                      left_on=keys_merged_fathmm, 
                      right_on=keys_fathmm, 
                      how='left')

# Restore original order
merged = merged.sort_values('_orig_index').drop(columns=['_orig_index'])

# Clean up temporary normalization columns
columns_to_drop = ['chrom_norm', 'pos_norm', 'ref_norm', 'alt_norm']
columns_to_drop = [col for col in columns_to_drop if col in merged.columns]
merged = merged.drop(columns=columns_to_drop)

out_path = 'Fathmm_Tumor_merged.csv'
merged.to_csv(out_path, index=False)

# Quick verification stats
print('\n=== Merge Statistics ===')
print(f'Total rows: {len(merged)}')

# Check pervariant merge
if any(col in merged.columns for col in pervar_cols):
    pervar_matched = merged[pervar_cols[0]].notna().sum() if len(pervar_cols) > 0 else 0
    print(f'Rows with pervariant data: {pervar_matched}')

# Check FATHMM merge
fathmm_matched = merged['fathmm_noncoding_score'].notna().sum()
print(f'Rows with FATHMM predictions: {fathmm_matched}')
print(f'Rows with FATHMM coding scores: {merged["fathmm_coding_score"].notna().sum()}')
print(f'Rows with FATHMM non-coding scores: {merged["fathmm_noncoding_score"].notna().sum()}')

print(f'\nOutput saved to: {out_path}')