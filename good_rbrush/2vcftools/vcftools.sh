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

vcftools --vcf /export/martinsons/manuela/goods_rbrush/2stacks/populations.snps.vcf --freq2 --out good2_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/goods_rbrush/2stacks/populations.snps.vcf --depth --out good2_depth
vcftools --vcf /export/martinsons/manuela/goods_rbrush/2stacks/populations.snps.vcf --site-mean-depth --out good2_sitedepth
vcftools --vcf /export/martinsons/manuela/goods_rbrush/2stacks/populations.snps.vcf --missing-indv --out good2_miss_indiv
vcftools --vcf /export/martinsons/manuela/goods_rbrush/2stacks/populations.snps.vcf --missing-site --out good2_miss_site
