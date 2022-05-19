#!/bin/bash
#SBATCH -t 0-23:30
#SBATCH -p thinnodes
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --mem=80GB

# This script will run WhatsHap for phase set tagging on the selected chromosome and population pair

# Usage:
# pop_chr_vcf_whatshap.sh <population> <chromosome>

module load cesga/2020
module load whatshap/1.1

# input directory
INdir=/mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression

# reference genome
ref=/mnt/netapp1/Store_csebdjgl/reference_genomes/Ref_Genome_LyCa/lc4.fa

# population
pop=($(echo "${1}"))

# population bam names
bam_pops=($(grep "${pop}" ${INdir}/lp_ll_introgression_populations.txt | cut -f1))

# chromosome
chr=($(echo "${2}"))


# bam populations of pop

# WhatsHap for PS
echo "Phase sets of $chr in $sp"
whatshap phase \
 -o ${INdir}/lp_ll_introgression_filtered_${pop}_${chr}_ps.vcf \
 --tag=PS \
 --reference=${ref} \
 ${INdir}/lp_ll_introgression_filtered_${pop}_${chr}.vcf \
 $(ls /mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_bams/*_LyCa_ref_sorted_rg_rmdup_sorted_indelrealigner.bam | 
    grep -f <(echo "${bam_pops[@]}" | tr ' ' '\n'))

