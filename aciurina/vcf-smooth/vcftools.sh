#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=vcf
#SBATCH --output=vcf_out_%j
#SBATCH --error=vcf_err_%j
#SBATCH --partition=ceti
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mlondono@unm.edu

cd /export/martinsons/manuela/smooth/vcf

module load miniconda3
source activate stacks-env

#denovo4
vcftools --vcf smooth.snps.vcf --freq2 --out smooth_alleles --max-alleles 2
vcftools --vcf smooth.snps.vcf --depth --out smooth_depth
vcftools --vcf smooth.snps.vcf --site-mean-depth --out smooth_sitedepth
vcftools --vcf smooth.snps.vcf --missing-indv --out smooth_miss_indiv
vcftools --vcf smooth.snps.vcf --missing-site --out smooth_miss_site
~                                                                               
~                         
