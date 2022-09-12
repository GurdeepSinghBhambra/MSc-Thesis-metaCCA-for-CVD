import pandas as pd 


# preprocessing multi gene annotation to single gene
annot_df = pd.read_csv('./data/processed/annotated_lookup.annot', delimiter=" ")

filtered_annot_df = annot_df[annot_df['ANNOT'].apply(lambda x: "|" not in x)]

filtered_annot_df.to_csv('./data/processed/annotated_lookup_singlehg19genes.txt', sep="\t", index=False)

print("Single SNP filter INFO:")
print("File:", "Annotated Single HG19 genes lookup table")
print("SNPs:", filtered_annot_df.shape[0], "/", annot_df.shape[0])
print("Genes:", filtered_annot_df['ANNOT'].nunique(), annot_df['ANNOT'].nunique())
print("File saved at:", "./data/processed/annotated_lookup_singlehg19genes.txt")
print("\n\n")


# filtering SXY_study1
df = pd.read_csv('./data/processed/sxy_study_1_maf1_pruned.txt', delimiter="\t")

filtered_df = df[df['SNPID_UKB'].isin(filtered_annot_df['SNP'])]

filtered_df.to_csv('./data/processed/sxy_study_1_maf1_pruned_singlehg19genes.txt', sep="\t", index=False)

print("Single SNP filter INFO:")
print("File:", "SXY_STUDY_1 HG19 filtered genes")
print("SNPs:", filtered_df.shape[0], "/", df.shape[0])
# print("Genes:", filtered_df['ANNOT'].nunique(), filtered_df['ANNOT'].nunique())
print("File saved at:", "./data/processed/sxy_study_1_maf1_pruned_singlehg19genes.txt")
print("\n\n")


# filtering SXX
df = pd.read_csv('./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted_nohead_ld02_pruned.txt', delimiter="\t", header=None)

filtered_df = df[df[2].isin(filtered_annot_df['SNP'])]

filtered_df.to_csv('./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted_nohead_ld02_pruned_singlehg19genes.txt', sep="\t", index=False, header=False)

print("Single SNP filter INFO:")
print("File:", "SXX_STUDY_1 HG19 filtered genes")
print("SNPs:", filtered_df.shape[0], "/", annot_df.shape[0])
# print("Genes:", filtered_df['ANNOT'].nunique(), filtered_df['ANNOT'].nunique())
print("File saved at:", "./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted_nohead_ld02_pruned_singlehg19genes.txt")
print("\n\n")

