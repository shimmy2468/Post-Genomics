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

# 1.3 concatenate 
Final_Normal = pd.read_csv('Final_Normal.csv')
Final_Tumor = pd.read_csv('Final_Tumor.csv')

AML_Expand = pd.concat([Final_Normal, Final_Tumor], axis=0)
AML_Expand = AML_Expand.reset_index(drop=True)

# Define CSQ columns
CSQ_columns = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "RefSeq", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET", "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF", "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF", "ExAC_EAS_MAF", "ExAC_FIN_MAF", "ExAC_NFE_MAF", "ExAC_OTH_MAF", "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "ENTREZ", "EVIDENCE"]

# Clean and process the CSQ columns
for col in CSQ_columns:
    # Convert string representations of lists to actual lists
    AML_Expand[col] = AML_Expand[col].apply(lambda x: eval(x) if isinstance(x, str) and x.startswith('[') else x)

# Create a list of records for the expanded dataframe
expanded_records = []
non_csq_columns = [col for col in AML_Expand.columns if col not in CSQ_columns]

# Iterate through each row and create expanded records
for _, row in AML_Expand.iterrows():
    # Get the number of elements in the first non-empty CSQ column
    csq_lengths = [len(row[col]) for col in CSQ_columns if isinstance(row[col], list)]
    if not csq_lengths:
        continue
    n_elements = csq_lengths[0]
    
    # Create a record for each element
    for i in range(n_elements):
        new_record = {}
        # Copy non-CSQ columns
        for col in non_csq_columns:
            new_record[col] = row[col]
        # Extract i-th element from each CSQ column
        for col in CSQ_columns:
            if isinstance(row[col], list) and len(row[col]) > i:
                # Clean the string by removing quotes and extra spaces
                value = row[col][i].strip().strip("'").strip()
                new_record[col] = value
            else:
                new_record[col] = ''
        expanded_records.append(new_record)

# Create new dataframe from expanded records
AML_Expand = pd.DataFrame(expanded_records)

# Remove duplicate rows
AML_Expand = AML_Expand.drop_duplicates()

# Reset index after removing duplicates
AML_Expand = AML_Expand.reset_index(drop=True)

# Save to CSV
AML_Expand.to_csv('AML_Expand.csv', index=False)
print(AML_Expand)

# 1.3.2 Create AML_gene
AML_gene = AML_Expand[['SYMBOL','Gene','Feature']].drop_duplicates().reset_index(drop=True)
AML_gene.to_csv('AML_gene.csv', index=False)
print(AML_gene)

# 1.3.3 AML_tx
AML_tx = AML_Expand[["chrom", "left", "right", "ref_seq", "alt_seq", "Feature", "cDNA_position", "BIOTYPE"]].drop_duplicates().reset_index(drop=True)
AML_tx.to_csv('AML_tx.csv', index=False)
print(AML_tx)


# Part 2
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import load_iris

iris = load_iris()
df = pd.DataFrame(data=iris.data, columns=iris.feature_names)
df['target'] = iris.target
print(df.head())

X = df.iloc[:, :-1].values
y = df.iloc[:, -1].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# For 2.3 I reduced number of trees from 100 to 10 and limited max depth to 2
classifier = RandomForestClassifier(n_estimators=10, max_depth=2, random_state=42)
print(f'RandomForest parameters: n_estimators={classifier.n_estimators}, max_depth={classifier.max_depth}, random_state={classifier.random_state}')
classifier.fit(X_train, y_train)
y_pred = classifier.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)
print(f'Accuracy: {accuracy * 100:.2f}%')

conf_matrix = confusion_matrix(y_test, y_pred)

plt.figure(figsize=(8, 6))
sns.heatmap(conf_matrix, annot=True, fmt='g', cmap='Blues', cbar=False, 
            xticklabels=iris.target_names, yticklabels=iris.target_names)

plt.title('Confusion Matrix Heatmap')
plt.xlabel('Predicted Labels')
plt.ylabel('True Labels')
plt.show()

feature_importances = classifier.feature_importances_

plt.barh(iris.feature_names, feature_importances)
plt.xlabel('Feature Importance')
plt.title('Feature Importance in Random Forest Classifier')
plt.show()

# Part 3
import numpy as np
from sklearn.datasets import load_digits

data, labels = load_digits(return_X_y=True)
(n_samples, n_features), n_digits = data.shape, np.unique(labels).size

print(f"# digits: {n_digits}; # samples: {n_samples}; # features {n_features}")

from time import time

from sklearn import metrics
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def bench_k_means(kmeans, name, data, labels):
    """Benchmark to evaluate the KMeans initialization methods.

    Parameters
    ----------
    kmeans : KMeans instance
        A :class:`~sklearn.cluster.KMeans` instance with the initialization
        already set.
    name : str
        Name given to the strategy. It will be used to show the results in a
        table.
    data : ndarray of shape (n_samples, n_features)
        The data to cluster.
    labels : ndarray of shape (n_samples,)
        The labels used to compute the clustering metrics which requires some
        supervision.
    """
    t0 = time()
    estimator = make_pipeline(StandardScaler(), kmeans).fit(data)
    fit_time = time() - t0
    results = [name, fit_time, estimator[-1].inertia_]

    # Define the metrics which require only the true labels and estimator
    # labels
    clustering_metrics = [
        metrics.homogeneity_score,
        metrics.completeness_score,
        metrics.v_measure_score,
        metrics.adjusted_rand_score,
        metrics.adjusted_mutual_info_score,
    ]
    results += [m(labels, estimator[-1].labels_) for m in clustering_metrics]

    # The silhouette score requires the full dataset
    results += [
        metrics.silhouette_score(
            data,
            estimator[-1].labels_,
            metric="euclidean",
            sample_size=300,
        )
    ]

    # Show the results
    formatter_result = (
        "{:9s}\t{:.3f}s\t{:.0f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}"
    )
    print(formatter_result.format(*results))

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

print(82 * "_")
print("init\t\ttime\tinertia\thomo\tcompl\tv-meas\tARI\tAMI\tsilhouette")

kmeans = KMeans(init="k-means++", n_clusters=n_digits, n_init=4, random_state=0)
bench_k_means(kmeans=kmeans, name="k-means++", data=data, labels=labels)

kmeans = KMeans(init="random", n_clusters=n_digits, n_init=4, random_state=0)
bench_k_means(kmeans=kmeans, name="random", data=data, labels=labels)

pca = PCA(n_components=n_digits).fit(data)
kmeans = KMeans(init=pca.components_, n_clusters=n_digits, n_init=1)
bench_k_means(kmeans=kmeans, name="PCA-based", data=data, labels=labels)

print(82 * "_")

reduced_data = PCA(n_components=2).fit_transform(data)
kmeans = KMeans(init="k-means++", n_clusters=n_digits, n_init=4)
kmeans.fit(reduced_data)

# Step size of the mesh. Decrease to increase the quality of the VQ.
h = 0.02  # point in the mesh [x_min, x_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = reduced_data[:, 0].min() - 1, reduced_data[:, 0].max() + 1
y_min, y_max = reduced_data[:, 1].min() - 1, reduced_data[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

# Obtain labels for each point in mesh. Use last trained model.
Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
plt.figure(1)
plt.clf()
plt.imshow(
    Z,
    interpolation="nearest",
    extent=(xx.min(), xx.max(), yy.min(), yy.max()),
    cmap=plt.cm.Paired,
    aspect="auto",
    origin="lower",
)

plt.plot(reduced_data[:, 0], reduced_data[:, 1], "k.", markersize=2)
# Plot the centroids as a white X
centroids = kmeans.cluster_centers_
plt.scatter(
    centroids[:, 0],
    centroids[:, 1],
    marker="x",
    s=169,
    linewidths=3,
    color="w",
    zorder=10,
)
plt.title(
    "K-means clustering on the digits dataset (PCA-reduced data)\n"
    "Centroids are marked with white cross"
)
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xticks(())
plt.yticks(())
plt.show()