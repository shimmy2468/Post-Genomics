import pandas as pd
from pathlib import Path


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

print('Merge keys chosen:', keys_final)

# Select columns from pervar to merge - exclude columns used for joining
pervar_cols = [c for c in pervar.columns if c not in (keys_pervar + ['_orig_index'])]

merged = final.merge(pervar[pervar_cols + keys_pervar], left_on=keys_final, right_on=keys_pervar, how='left')

# Restore original order
merged = merged.sort_values('_orig_index').drop(columns=['_orig_index'])

out_path = 'Final_Tumor_merged.csv'
merged.to_csv(out_path, index=False)

# Quick verification stats
matches = merged[keys_final].notna().all(axis=1).sum()
total = len(merged)
print(f'Matched rows (non-null join keys): approx {matches} / {total}')

# Most of the code references 'norm', but I just change the paths for the normal and tumor files respectively with the script remaining the same