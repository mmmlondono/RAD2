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

cd /export/martinsons/manuela/mfc/vcf

module load miniconda3
source activate stacks-env

#denovo4
vcftools --vcf mfc.snps.vcf --freq2 --out mfc_alleles --max-alleles 2
vcftools --vcf mfc.snps.vcf --depth --out mfc_depth
vcftools --vcf mfc.snps.vcf --site-mean-depth --out mfc_sitedepth
vcftools --vcf mfc.snps.vcf --missing-indv --out mfc_miss_indiv
vcftools --vcf mfc.snps.vcf --missing-site --out mfc_miss_site
~                                                                               
~                         
