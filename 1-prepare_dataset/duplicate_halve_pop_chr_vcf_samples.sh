#!/bin/bash
# Usage:
# pop_chr_vcf_whatshap.sh <switch> <population> <chromosome>

# switch for duplicate or halve
switch=($(echo "${1}"))

# input directory
INdir=/mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression

# population
pop=($(echo "${2}"))

# chromosome
chr=($(echo "${3}"))

# prefix
PREfix=${INdir}/lp_ll_introgression_filtered_${pop}_${chr}_ps

# number of samples
nsamples=($(grep -m1 "#CHR" ${PREfix}.vcf | tr '\t' '\n' | grep "_" | wc -l))

if [ $switch == "duplicate" ]
 then 
  
  # copy the header except for the table line
  echo "creating duplicated VCF header"
  grep "##" ${PREfix}.vcf > ${PREfix}_duplicate.vcf
  
  # create a new VCF table header (#CHR line)
  # with double the samples (new ones are called the same + "_2")
  # sed command removes trailing \t 
  echo "adding samples names"
  paste <(grep "#CHR" ${PREfix}.vcf) \
   <(paste -d '_' <(grep "#CHR" ${PREfix}.vcf | tr '\t' '\n' | grep "_") <(yes "2" | head -n ${nsamples}) | 
     tr '\n' '\t') | sed 's/^[ \t]*//;s/[ \t]*$//' \
  >> ${PREfix}_duplicate.vcf
  
  # extract extract sample's genotypes
  echo "extracting genotypes of ${nsamples} samples"
  grep -v "#" ${PREfix}.vcf | rev | cut -f1-${nsamples} | rev > tmp.gts
    
  
  # paste the vcf without the header with the genotype and add it to the duplicated vcf
  echo "pasting the GTs and completing duplicated VCF"
  paste <(grep -v "#" ${PREfix}.vcf) <(cat tmp.gts) >> ${PREfix}_duplicate.vcf
  
  # remove temporary file
  echo "removing tmp file"
  rm tmp.gts

elif [ $switch == "halve" ]
 then
 
fi
