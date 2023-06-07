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

#denovo1
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks1/populations.snps.vcf --freq2 --out 1_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks1/populations.snps.vcf --depth --out 1_depth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks1/populations.snps.vcf --site-mean-depth --out 1_sitedepth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks1/populations.snps.vcf --missing-indv --out 1_miss_indiv
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks1/populations.snps.vcf --missing-site --out 1_miss_site

#denovo2
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks2/populations.snps.vcf --freq2 --out 2_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks2/populations.snps.vcf --depth --out 2_depth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks2/populations.snps.vcf --site-mean-depth --out 2_sitedepth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks2/populations.snps.vcf --missing-indv --out 2_miss_indiv
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks2/populations.snps.vcf --missing-site --out 2_miss_site

#denovo3
#vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks3/populations.snps.vcf --freq2 --out 3_alleles --max-alleles 2
#vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks3/populations.snps.vcf --depth --out 3_depth
#vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks3/populations.snps.vcf --site-mean-depth --out 3_sitedepth
#vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks3/populations.snps.vcf --missing-indv --out 3_miss_indiv
#vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks3/populations.snps.vcf --missing-site --out 3_miss_site

#denovo4
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks4/populations.snps.vcf --freq2 --out 4_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks4/populations.snps.vcf --depth --out 4_depth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks4/populations.snps.vcf --site-mean-depth --out 4_sitedepth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks4/populations.snps.vcf --missing-indv --out 4_miss_indiv
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks4/populations.snps.vcf --missing-site --out 4_miss_site

#denovo5
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks5/populations.snps.vcf --freq2 --out 5_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks5/populations.snps.vcf --depth --out 5_depth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks5/populations.snps.vcf --site-mean-depth --out 5_sitedepth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks5/populations.snps.vcf --missing-indv --out 5_miss_indiv
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks5/populations.snps.vcf --missing-site --out 5_miss_site

#denovo6
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks6/populations.snps.vcf --freq2 --out 6_alleles --max-alleles 2
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks6/populations.snps.vcf --depth --out 6_depth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks6/populations.snps.vcf --site-mean-depth --out 6_sitedepth
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks6/populations.snps.vcf --missing-indv --out 6_miss_indiv
vcftools --vcf /export/martinsons/manuela/cfs_rbrush/stacks6/populations.snps.vcf --missing-site --out 6_miss_site
