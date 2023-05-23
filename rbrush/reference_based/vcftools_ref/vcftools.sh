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

#denovo 
vcftools --vcf /export/martinsons/manuela/rbrush/reference/populations_refbased_output/populations.snps.vcf --freq2 --out ref_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/rbrush/reference/populations_refbased_output/populations.snps.vcf --depth --out ref_depth
vcftools --vcf /export/martinsons/manuela/rbrush/reference/populations_refbased_output/populations.snps.vcf --site-mean-depth --out ref_sitedepth
vcftools --vcf /export/martinsons/manuela/rbrush/reference/populations_refbased_output/populations.snps.vcf --missing-indv --out ref_miss_indiv
vcftools --vcf /export/martinsons/manuela/rbrush/reference/populations_refbased_output/populations.snps.vcf --missing-site --out ref_miss_site

