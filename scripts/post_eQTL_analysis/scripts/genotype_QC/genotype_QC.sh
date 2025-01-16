#!/bin/bash
# Genotype Data Quality Control (QC) Pipeline


# Define the name of the study
study="ALL"

##############################################################
#                 STEP 1: INITIAL PREPARATION                #
##############################################################

# Convert input SNP data to PLINK binary format (.bed, .bim, .fam)
# --allow-extra-chr: Allows non-standard chromosomes
plink --file archivos_axiom/${study} --make-bed --out ${study} --allow-extra-chr

# Check missingness per individual and per SNP
echo "Checking missingness..."
plink --bfile ${study} --missing --allow-extra-chr --out miss

# Check sex discrepancies based on genotype data
echo "Checking sex discrepancies..."
plink --bfile ${study} --check-sex --allow-extra-chr --out sex

# Identify problematic sex discrepancies
grep PROBLEM sex.sexcheck > sex.drop

##############################################################
#                    STEP 2: VARIANT QC                     #
##############################################################

# Investigate missingness again (post-sex check)
plink --bfile ${study} --missing --allow-extra-chr --out miss

# Generate missingness histograms in R
echo "Generating missingness histograms..."
Rscript --no-save hist_miss.R

# Remove SNPs with high genotype missingness (>5%)
echo "Filtering SNPs with >5% missing genotypes..."
plink --bfile ${study} --geno 0.05 --make-bed --allow-extra-chr --out ${study}_geno

# Perform case-control rate difference QC
echo "Running case-control missingness test..."
plink --bfile ${study}_geno --test-missing --allow-extra-chr --out case-control

# Filter case-control SNPs with p-value < 1e-5
echo "Filtering case-control SNPs..."
Rscript --no-save case_control.R
plink --bfile ${study}_geno --exclude case_control.drop --make-bed --allow-extra-chr --out ${study}_case_control

# Minor Allele Frequency (MAF) check
echo "Checking MAF distribution..."
plink --bfile ${study}_case_control --freq --allow-extra-chr --out MAF_check
Rscript --no-save MAF_check.R

# Remove SNPs with MAF < 0.05
echo "Filtering SNPs with MAF < 0.05..."
plink --bfile ${study}_case_control --maf 0.05 --make-bed --allow-extra-chr --out ${study}_MAF

# Hardy-Weinberg Equilibrium (HWE) check
echo "Running Hardy-Weinberg Equilibrium check..."
plink --bfile ${study}_MAF --hardy midp --allow-extra-chr --out hwe
Rscript --no-save hwe_zoom.R
Rscript --no-save hwe.R

# Filter controls for HWE (p < 1e-6)
plink --bfile ${study}_MAF --hwe 1e-6 --make-bed --allow-extra-chr --out ${study}_hwe_controls

# Filter cases for extreme HWE deviations (p < 1e-10)
plink --bfile ${study}_hwe_controls --hwe 1e-10 --hwe-all --make-bed --allow-extra-chr --out ${study}_hwe

##############################################################
#                    STEP 3: SAMPLE QC                      #
##############################################################

# Check missingness per individual
echo "Checking missingness per individual..."
plink --bfile ${study}_hwe --missing --allow-extra-chr --out miss_post

# Remove individuals with missingness > 5%
plink --bfile ${study}_hwe --mind 0.05 --make-bed --allow-extra-chr --out ${study}_mind

# Check heterozygosity rates
echo "Checking heterozygosity rates..."
plink --bfile ${study}_mind --het --allow-extra-chr --out het
Rscript --no-save check_heterozygosity_rate.R
Rscript --no-save heterozygosity_outliers_list.R

# Remove individuals with heterozygosity rate >3SD from the mean
plink --bfile ${study}_mind --remove fail-het-qc_remove.txt --make-bed --allow-extra-chr --out ${study}_het

# Generate LD-pruned SNP set for relatedness and population analysis
plink --bfile ${study}_het --indep-pairwise 1000 5 0.2 --allow-extra-chr --out indep

##############################################################
#                 STEP 4: RELATEDNESS CHECK                 #
##############################################################

# Identity-by-Descent (IBD) to identify related individuals
echo "Calculating relatedness (IBD)..."
plink --bfile ${study}_het --extract indep.prune.in --genome --min 0.2 --out pihat_min0.2 --allow-extra-chr

# Remove related individuals
Rscript --no-save cryptic_relatedness_deletion.R
plink --bfile ${study}_het --remove pihat.drop --make-bed --allow-extra-chr --out ${study}_pihat

##############################################################
#            STEP 5: POPULATION STRATIFICATION              #
##############################################################

# Merge dataset with 1000 Genomes reference data for PCA
# Files 1kG_MDS5* were made by the bioinformatics platform in a previous analysis
# Extract shared SNPs
awk '{print $2}' ${study}_pihat.bim > ${study}_SNPs.txt
plink --bfile 1kG_MDS5 --extract ${study}_SNPs.txt --make-bed --out 1kG_MDS6

# Align builds and reference alleles
awk '{print $2,$4}' ${study}_MDS.map > ${study}.txt
plink --bfile 1kG_MDS6 --update-map ${study}.txt --make-bed --out 1kG_MDS7
awk '{print$2,$5}' 1kG_MDS7.bim > 1kg_ref-list.txt
plink --bfile ${study}_MDS --reference-allele 1kg_ref-list.txt --make-bed --out ${study}-adj

# Resolve strand issues and flip SNPs
awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp
awk '{print$2,$5,$6}' ${study}-adj.bim > ${study}-adj_tmp
sort 1kGMDS7_tmp ${study}-adj_tmp |uniq -u > all_differences.txt
awk '{print $1}' all_differences.txt | sort -u > flip_list.txt
plink --bfile ${study}-adj --flip flip_list.txt --make-bed --allow-extra-chr --out ${study}_flipped

# Exclude problematic SNPs
awk '{print$2,$5,$6}' ${study}_flipped.bim > ${study}_flipped_tmp
sort 1kGMDS7_tmp ${study}_flipped_tmp |uniq -u  > uncorresponding_SNPs.txt
awk '{print $1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt
plink --bfile ${study}_flipped --exclude SNPs_for_exclusion.txt --make-bed --out ${study}_MDS2
plink --bfile 1kG_MDS7 --exclude SNPs_for_exclusion.txt --make-bed --out 1kG_MDS8

# Merge datasets: study and 1000 Genomes Data
plink --bfile ${study}_MDS2 --bmerge 1kG_MDS8 --allow-no-sex --make-bed --out MDS_merge2

# Perform MDS analysis
plink --bfile MDS_merge2 --extract indep.prune.in --genome --out MDS_merge2
plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2

# Download the file with population information of the 1000 genomes dataset
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel

# Convert population codes into superpopulation codes (i.e., AFR,AMR,ASN, and EUR)
awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt
sed 's/JPT/ASN/g' race_1kG.txt>race_1kG2.txt
sed 's/ASW/AFR/g' race_1kG2.txt>race_1kG3.txt
sed 's/CEU/EUR/g' race_1kG3.txt>race_1kG4.txt
sed 's/CHB/ASN/g' race_1kG4.txt>race_1kG5.txt
sed 's/CHD/ASN/g' race_1kG5.txt>race_1kG6.txt
sed 's/YRI/AFR/g' race_1kG6.txt>race_1kG7.txt
sed 's/LWK/AFR/g' race_1kG7.txt>race_1kG8.txt
sed 's/TSI/EUR/g' race_1kG8.txt>race_1kG9.txt
sed 's/MXL/AMR/g' race_1kG9.txt>race_1kG10.txt
sed 's/GBR/EUR/g' race_1kG10.txt>race_1kG11.txt
sed 's/FIN/EUR/g' race_1kG11.txt>race_1kG12.txt
sed 's/CHS/ASN/g' race_1kG12.txt>race_1kG13.txt
sed 's/PUR/AMR/g' race_1kG13.txt>race_1kG14.txt

# Create a racefile of your own data
awk '{print$1,$2,"OWN"}' ${study}_MDS.fam>racefile_own.txt

# Concatenate racefiles
cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

# Generate MDS plot
Rscript MDS_merged.R

# Exclude ethnic outliers based on MDS
Rscript --no-save EUR_MDS_merge2.R
plink --bfile ${study}_pihat --keep EUR_MDS_merge2 --make-bed --out ${study}_popStrat

# Convert file to vcf for downstream analysis 
plink --bfile ${study}_popStrat --recode vcf-iid --out ${study}_passed_QC

echo "Pipeline completed successfully!"
