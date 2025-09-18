import pandas as pd
import os

#normal
normal_folder = "Normal_CSV"
normal_files = [file for file in os.listdir(normal_folder)]

normal_dfs = []
for file in normal_files:
    path = os.path.join(normal_folder, file)
    normal_dfs.append(pd.read_csv(path))
merged_normal_df = pd.concat(normal_dfs, ignore_index=True)

#tumor
tumor_folder = "Tumor_CSV"
tumor_files = [file for file in os.listdir(tumor_folder)]

tumor_dfs = []
for file in tumor_files:
    path = os.path.join(tumor_folder, file)
    tumor_dfs.append(pd.read_csv(path))

merged_tumor_df = pd.concat(tumor_dfs, ignore_index=True)

'''
Homework 3.A
(i) Merge the 5 normal CSV files together and the 5 tumor CSV files, 
result should 2 separate Dataframes, one with Normal variants and another with Tumor variants.

SUGGESTION: Prior to coding, create 2 empty folders, Normal_CSV and Tumor_CSV, 
manually move the 5 normal CSVs into the Normal_CSV folder, and then move
the 5 tumor CSVs into the Tumor_CSV folder, this can be done 
by using the search bar in the Finder(Mac) or Folder(Windows) app. 
The script can then point to the directory (similar to HW1) to read and merge the files within, 
using a function within the pandas (pd) package.

Reading in a CSV file Example:
DataFrame1 = pd.read_csv("DataFrame1.csv")

Merging Example:
newDataFrame = pd.concat(DataFrame1, DataFrame2, axis=0) 
'''


# Function adds in alt_seq column to, input is a dataframe and function returns a dataframe
def addALT_Seq(csv):
    alt = []
    for row in range(csv.shape[0]):
        ref_seq = csv["ref_seq"][row]
        if ref_seq == csv["var_seq1"][row]:
            alt.append(csv["var_seq2"][row])
        else:
            alt.append(csv["var_seq1"][row])
    csv.insert(csv.shape[1], "alt_seq", alt)
    return csv

'''
Homework 3.A
(ii) Using the output from A(i), run the addALT_Seq() function:
Example:
newDataFrame_withALTseq = addALT_Seq(NewDataFrame)
'''
normal_merged_withALTseq = addALT_Seq(merged_normal_df) 

tumor_merged_withALTseq = addALT_Seq(merged_tumor_df)

'''
Homework 3.A
(iii) Using the output from A(ii), remove duplicates based on the given columns:
[“chrom”, “left”, “ref_seq”, “alt_seq”, “Patient_ID”]
Save the two DataFrames as: Final_Normal and Final_Tumor

Remove Duplicates Example:
Final = newDataFrame_withALTseq.drop_duplicates(columns)
'''

Final_Normal = normal_merged_withALTseq.drop_duplicates(subset=["chrom", "left", "ref_seq", "alt_seq", "Patient_ID"])
Final_Normal.to_csv("Final_Normal.csv", index=False)

Final_Tumor = tumor_merged_withALTseq.drop_duplicates(subset=["chrom", "left", "ref_seq", "alt_seq", "Patient_ID"])
Final_Tumor.to_csv("Final_Tumor.csv", index=False)

'''
OUTPUT CHECK
Homework 3.A
(iv) Run the lines below:
'''
print("The number of (Rows, Columns) in Tumor:")
print(Final_Tumor.shape)
print("The number of (Rows, Columns) in Normal:")
print(Final_Normal.shape)






