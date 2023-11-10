#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=vcf
#SBATCH --output=vcfout_%j
#SBATCH --error=vcferr_%j
#SBATCH --partition=ceti
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mlondono@unm.edu

module load miniconda3
source activate stacks-env

#Stacks1
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks1/populations.snps.vcf --freq2 --out erna1_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks1/populations.snps.vcf --depth --out erna1_depth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks1/populations.snps.vcf --site-mean-depth --out erna1_sitedepth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks1/populations.snps.vcf --missing-indv --out erna1_miss_indiv
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks1/populations.snps.vcf --missing-site --out erna_miss_site

#Stacks2
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks2/populations.snps.vcf --freq2 --out erna2_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks2/populations.snps.vcf --depth --out erna2_depth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks2/populations.snps.vcf --site-mean-depth --out erna2_sitedepth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks2/populations.snps.vcf --missing-indv --out erna2_miss_indiv
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks2/populations.snps.vcf --missing-site --out erna2_miss_site

#Stacks3
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks3/populations.snps.vcf --freq2 --out erna3_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks3/populations.snps.vcf --depth --out erna3_depth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks3/populations.snps.vcf --site-mean-depth --out erna3_sitedepth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks3/populations.snps.vcf --missing-indv --out erna3_miss_indiv
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks3/populations.snps.vcf --missing-site --out erna3_miss_site

#Stacks4
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks4/populations.snps.vcf --freq2 --out erna4_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks4/populations.snps.vcf --depth --out erna4_depth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks4/populations.snps.vcf --site-mean-depth --out erna4_sitedepth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks4/populations.snps.vcf --missing-indv --out erna4_miss_indiv
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks4/populations.snps.vcf --missing-site --out erna4_miss_site

#Stacks5
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks5/populations.snps.vcf --freq2 --out erna5_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks5/populations.snps.vcf --depth --out erna5_depth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks5/populations.snps.vcf --site-mean-depth --out erna5_sitedepth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks5/populations.snps.vcf --missing-indv --out erna5_miss_indiv
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks5/populations.snps.vcf --missing-site --out erna5_miss_site

#Stacks6
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6/populations.snps.vcf --freq2 --out erna6_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6/populations.snps.vcf --depth --out erna6_depth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6/populations.snps.vcf --site-mean-depth --out erna6_sitedepth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6/populations.snps.vcf --missing-indv --out erna6_miss_indiv
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6/populations.snps.vcf --missing-site --out erna6_miss_site

#Stacks6R
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6R/populations.snps.vcf --freq2 --out erna6R_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6R/populations.snps.vcf --depth --out erna6R_depth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6R/populations.snps.vcf --site-mean-depth --out erna6R_sitedepth
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6R/populations.snps.vcf --missing-indv --out erna6R_miss_indiv
vcftools --vcf /export/martinsons/manuela/erna/alone_erna/stacks6R/populations.snps.vcf --missing-site --out erna6R_miss_site
