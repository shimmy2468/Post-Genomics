import pandas as pd

df = pd.read_csv('AML_gene.csv')

# get unique values of the gene column
unique_genes = df['Gene'].dropna().astype(str).unique()

# unique genes into a txt file, one per line
with open('unique_genes.txt', 'w', encoding='utf-8') as f:
    for gene in unique_genes:
        f.write(gene + '\n')


