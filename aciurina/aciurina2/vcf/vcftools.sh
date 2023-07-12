#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=vcf_test1
#SBATCH --output=out_%j
#SBATCH --error=err_%j
#SBATCH --partition=ceti
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mlondono@unm.edu

module load miniconda3
source activate stacks-env

#denovo with new popmap
vcftools --vcf /export/martinsons/manuela/aciurina/wthnewpopmap/stacks/populations.snps.vcf --freq2 --out aciu2_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/aciurina/wthnewpopmap/stacks/populations.snps.vcf --depth --out aciu2_depth
vcftools --vcf /export/martinsons/manuela/aciurina/wthnewpopmap/stacks/populations.snps.vcf --site-mean-depth --out aciu2_sitedepth
vcftools --vcf /export/martinsons/manuela/aciurina/wthnewpopmap/stacks/populations.snps.vcf --missing-indv --out aciu2_miss_indiv
vcftools --vcf /export/martinsons/manuela/aciurina/wthnewpopmap/stacks/populations.snps.vcf --missing-site --out aciu2_miss_site
