##
#
# References: 
# 1. Osborne A. data_preparation_summary_statistics.sh. https://github.com/AmyJaneOsborne/CCA_scripts/blob/main/data_preparation_summary_statistics.sh
# 2. Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MA, Bender D, Maller J, Sklar P, de Bakker PI, Daly MJ, Sham PC. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet. 2007 Sep;81(3):559-75. doi: 10.1086/519795. Epub 2007 Jul 25. PMID: 17701901; PMCID: PMC1950838.
#
##


# creating 1000g VCF file 
./plink --bfile './data/g1000p3_EUR/g1000p3_EUR' --recode vcf --out './data/g1000p3_EUR/g1000p3_EUR' --memory 13000 

# getting rsids into a file
awk '{print $3}' ./data/processed/maf1_data.vcf > ./data/processed/rsids.txt

# getting snps for 1000gp accorsing to the merged snps
bcftools view --include ID==@./data/processed/rsids.txt ./data/g1000p3_EUR/g1000p3_EUR.vcf > ./data/g1000p3_EUR/g1000p3_EUR_overlap.vcf

# counting cols
# awk '{print NF}' ./data/g1000p3_EUR/g1000p3_EUR_overlap.vcf | sort -nu | tail -n 1
# head -31 ./data/g1000p3_EUR/g1000p3_EUR_overlap.vcf | tail -1 | awk '{print NF}'
grep -m 1 -w '^#CHROM' ./data/g1000p3_EUR/g1000p3_EUR_overlap.vcf | awk '{print NF}'
# cols: 512

# counting snps, 7106549 SNPs
grep -v '#' ./data/g1000p3_EUR/g1000p3_EUR_overlap.vcf | wc -l

# literally data is put in g1000 data
awk 'NR==FNR{a[$3]=" "$0; next}{ print $0 (a[$3]?a[$3]:" missing")}' OFS='\t' ./data/processed/maf1_data.vcf ./data/g1000p3_EUR/g1000p3_EUR_overlap.vcf > ./data/g1000p3_EUR/g1000p3_EUR_common.vcf

# 7106580 SNPs
wc -l ./data/g1000p3_EUR/g1000p3_EUR_common.vcf 

# ref alt fixing
# COLUMN INDEX: REF @ 516, ALT @ 517
awk '{if ($4==toupper($517)) { gsub("0\\/0","2"); gsub("[12345]\\/[12345]","0"); gsub("0\\/[12345]","1"); gsub("[12345]\\/0","1"); print } }' ./data/g1000p3_EUR/g1000p3_EUR_common.vcf > ./data/g1000p3_EUR/g1000p3_EUR_common_refaltfixed.vcf

awk '{if (substr($5,1,1)==toupper($517)) {gsub("[02345]\\/[02345]","0"); gsub("1\\/1","2"); gsub("[02345]\\/1","1"); gsub("1\\/[02345]","1"); print } }' ./data/g1000p3_EUR/g1000p3_EUR_common.vcf > ./data/g1000p3_EUR/g1000p3_EUR_common_refaltOK_1.txt

awk '{if (substr($5,3,1)==toupper($517)) {gsub("[02345]\\/[02345]","0"); gsub("1\\/1","2"); gsub("[02345]\\/1","1"); gsub("1\\/[02345]","1"); print } }' ./data/g1000p3_EUR/g1000p3_EUR_common.vcf > ./data/g1000p3_EUR/g1000p3_EUR_common_refaltOK_2.txt

awk '{if (substr($5,5,1)==toupper($517)) {gsub("[02345]\\/[02345]","0"); gsub("1\\/1","2"); gsub("[02345]\\/1","1"); gsub("1\\/[02345]","1"); print } }' ./data/g1000p3_EUR/g1000p3_EUR_common.vcf > ./data/g1000p3_EUR/g1000p3_EUR_common_refaltOK_3.txt

# concatenating files
# removing headers and vsf data from g1000p3_EUR_common_refaltOK_1.txt
awk 'NR>28 {print $0}' ./data/g1000p3_EUR/g1000p3_EUR_common_refaltOK_1.txt > ./data/g1000p3_EUR/g1000p3_EUR_common_refaltOK_1_nohead.txt

# merging the files which were not empty
cat ./data/g1000p3_EUR/g1000p3_EUR_common_refaltfixed.vcf ./data/g1000p3_EUR/g1000p3_EUR_common_refaltOK_1_nohead.txt >> ./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted.txt

# 0 lines with (a number [0-5]/ a number [0-5])
grep "[0-5]/[0-5]" ./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted.txt | wc -l

# checking for lost snps, result: 44
awk 'NR==FNR{a[$3]=" "$0; next}{ print $0 (a[$3]?a[$3]:" absent")}' OFS='\t' ./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted.txt ./data/g1000p3_EUR/g1000p3_EUR_common.vcf | grep absent > ./data/missing_snps.txt


## LD Pruning
# Creating VCF file with header (just taking a clean header)
head -31 ./data/g1000p3_EUR/g1000p3_EUR_overlap.vcf > ./data/g1000p3_EUR/g1000p3_EUR_header.vcf

# cleaning data for pruning (removing header)
grep -v '#' ./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted.txt > ./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted_nohead.txt

# merging header vcf file and noheader 1000g file
cat ./data/g1000p3_EUR/g1000p3_EUR_header.vcf ./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted_nohead.txt >> ./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted_clean_head.vcf

mkdir ./data/g1000p3_EUR/prune

# creating binary files from clean head vcf
plink --vcf ./data/g1000p3_EUR/g1000p3_EUR_common_GTs_sorted_clean_head.vcf --memory 14000 --keep-allele-order --double-id --make-bed --out ./data/g1000p3_EUR/prune/g1000p3_EUR_common_GTs_sorted_clean_head_ld_bed

# pruning using binary files 
plink2 --bfile ./data/g1000p3_EUR/prune/g1000p3_EUR_common_GTs_sorted_clean_head_ld_bed --indep-pairwise 1000000 5 0.2 --threads 24 --memory 14000 --out ./data/g1000p3_EUR/prune/g1000p3_EUR_common_GTs_sorted_clean_head_ld_bed_1000k_02

# pruned snps: 389852 
wc -l ./data/g1000p3_EUR/prune/g1000p3_EUR_common_GTs_sorted_clean_head_ld_bed_1000k_02.prune.in 

# filtering pruned snps for GT file, sxy file and maf01_data vcf file
python sxy_sxx_vcf_filtering.py # Time taken: 4 mins, 30 secs | all filtering: OK


## Gene Annotation using plink
# creating file with .assoc format for annotating with plink and also which could be used later
# Format: ['CHR' 'SNP' 'BP' 'REF' 'ALT'], space delimited file

# getting pruned data only
awk '{ print $1 " " $3 " " $2 " " $4 " " $5}' ./data/processed/maf1_data_pruned.vcf > ./data/processed/for_gene_annot_noheader.txt

# merging the header along with data
echo 'CHR SNP BP REF ALT' > ./data/processed/for_gene_annot_header.txt

# merging header and data to the same data file
cat ./data/processed/for_gene_annot_header.txt ./data/processed/for_gene_annot_noheader.txt >> ./data/processed/for_gene_annot.txt

# plink annotation (this will create a file with multiple gene annotations)
plink --annotate ./data/processed/for_gene_annot.txt ranges=./data/plink_data/glist-hg19 prune minimal --out ./data/processed/annotated_lookup

# filter annotated snps in sxy and GT to reduce their size and filtering files according to snps with only single gene
python gene_filtering.py

# Final gene count: 15198
# Final SNP count: 166130

