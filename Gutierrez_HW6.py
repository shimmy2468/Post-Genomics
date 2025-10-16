import pandas as pd

# First part
Normal_subset = pd.read_csv('Final_Normal.csv')
Normal_subset = Normal_subset[["chrom", "left", "ref_seq", "alt_seq", "Patient_ID", "VCF_ID"]]
print(Normal_subset['Patient_ID'].unique()) # 5

Tumor_subset = pd.read_csv('Final_Tumor.csv')
Tumor_subset = Tumor_subset[["chrom", "left", "ref_seq", "alt_seq", "Patient_ID", "VCF_ID"]]
print(Tumor_subset['Patient_ID'].unique()) # 4

# 1.1.3
Normal_group = Normal_subset.groupby(["chrom", "left", "ref_seq", "alt_seq"]).agg(list).reset_index()
Tumor_group = Tumor_subset.groupby(["chrom", "left", "ref_seq", "alt_seq"]).agg(list).reset_index()


# 1.1.4
Normal_group['N#'] = Normal_subset.groupby(["chrom", "left", "ref_seq", "alt_seq"])['Patient_ID'].transform('nunique')
Tumor_group['T#'] = Tumor_subset.groupby(["chrom", "left", "ref_seq", "alt_seq"])['Patient_ID'].transform('nunique')


#1.1.5 renaming columns patient_id and vcf_id to have _normal or _tumor at the end
Normal_group.rename(columns={"Patient_ID": "Patient_ID_normal", "VCF_ID": "VCF_ID_normal"}, inplace=True)
Tumor_group.rename(columns={"Patient_ID": "Patient_ID_tumor", "VCF_ID": "VCF_ID_tumor"}, inplace=True)

print(Normal_group)
print(Tumor_group)

# 1.2 Merging to CSV
AML = pd.merge(Normal_group, Tumor_group, on=["chrom", "left", "ref_seq", "alt_seq"], how='outer')
AML.to_csv('AML.csv', index=False)

AML = pd.read_csv('AML.csv')

# get unique normal and tumor variants
Normal_variants = AML[AML['N#'].notna() & AML['T#'].isna()]
print(Normal_variants) # 0
Tumor_variants = AML[AML['T#'].notna() & AML['N#'].isna()]
print(Tumor_variants) # 1408
# shared variants
Shared_variants = AML[AML['N#'].notna() & AML['T#'].notna()]
print(Shared_variants) # 165

# 1.3 concatenate Final_Normal and Final_Tumor along the axis = 0, with this Expand/Explode the rows based on the CSQ columns and save this file as AML_Expand.csv. Remove duplicate rows.
Final_Normal = pd.read_csv('Final_Normal.csv')
Final_Tumor = pd.read_csv('Final_Tumor.csv')
AML_Expand = pd.concat([Final_Normal, Final_Tumor], axis=0)
AML_Expand = AML_Expand.explode('CSQ')
AML_Expand = AML_Expand.drop_duplicates()

print(AML_Expand)
