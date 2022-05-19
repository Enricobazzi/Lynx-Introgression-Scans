---
title: "1-prepare_dataset"
author: "Enrico"
date: "5/17/2022"
output: html_document
---

The genomic dataset we used for the introgression scans was comprised of all of the individuals we have sampled from the following populations:

-   *Lynx pardinus* :

    -   Sierra Morena (N=19)

-   *Lynx lynx* :

    -   Western clade Kirov + Urals (N=20)
    -   Eastern clade Yakutia + Primorsky krai (N=19)
    -   Southern clade Caucasus (N=12 or N=9)

I have created a table that stores this information (lp_ll_introgression_populations.txt) this way

| **sample** | **population** |
|:----------:|:--------------:|
|     sm     |      lpa       |
|     ki     |      wel       |
|     ur     |      wel       |
|     ya     |      eel       |
|     vl     |      eel       |
|     ca     |      sel       |

------------------------------------------------------------------------

## Variant Calling using GATK

I made a list of the bams to be included in the calling

```{bash}
# On CESGA FT2 server go to Canada lynx reference genome BAM files
cd /mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_bams

# Create a bamlist of samples from populations we want to include in our analysis
ls *er.bam | grep -E "ll|lp" | grep -E "ca|sm|ya|vl|ki|ur" | grep -vE "ca_0249|ca_0253" > lp_ll_introgression.bamlist
```

From our previous work we have already generated BAM files for each of these samples (yay!), generated using the Canada lynx assembly as the reference genome (mLynCan4_v1.p; GCA_007474595.1; Rhie et al., 2021).

To generate a VCF file from these BAMs, we performed variant calling on each scaffold in parallel, while additionally dividing all of the bigger scaffolds (18 autosomic + X) into 80 chunks each. This way more jobs have to be sbatched, but we have no memory or time problems.

We called scaffolds that could be called without being divided into chunks using the script perchr_haplotypecaller.sh

```{bash}
# short non-autosomic (??) scaffolds 
chr_list=($(cat /mnt/netapp1/Store_csebdjgl/reference_genomes/Ref_Genome_LyCa/non-autosomic_scaffolds_list.txt | grep -vE "mLynCan4_MT|Super_Scaffold_10"))
for chr in ${chr_list[@]}
 do
  echo "sbatching perchr_haplotypecaller for $chr"
  sbatch --mem=120GB perchr_haplotypecaller.sh $chr
done
```

Chunk divided scaffolds were then called using the script perchunk_haplotypecaller.sh

```{bash}
# big autosomic scaffolds (only 10 at the time for cluster job submission limitations):
chr_list=($(cat /mnt/netapp1/Store_csebdjgl/reference_genomes/Ref_Genome_LyCa/autosomic_scaffolds_list.txt))
for chr in ${chr_list[@]:0:10}
 do
  for i in {1..10}
   do
    echo "sbatching perchunk_haplotypecaller for $chr w $i"
    sbatch --mem=120GB perchunk_haplotypecaller.sh $chr $i
  done
done

# X chromosome (Super_Scaffold_10)
for i in {1..10}
 do
  echo "sbatching perchunk_haplotypecaller for Super_Scaffold_10 w $i"
  sbatch --mem=120GB perchunk_haplotypecaller.sh Super_Scaffold_10 $i
done
```

Resulting chunk and scaffold VCFs were concatenated using bcf-tools

```{bash}
# Create list of vcf files to concatenate
ls lp_ll_introgression_*caffold*.vcf | tr ' ' '\n' > lp_ll_introgression_vcfs.list

# Concatenate chromosome VCFs into whole genome VCF
bcftools concat -f lp_ll_introgression_vcfs.list \
 --output-type v \
 --output lp_ll_introgression_LyCa_ref.unsorted.vcf
```

Resulting VCF was then sorted and chunk/scaffold VCFs moved to a different folder of intermediary VCFs

```{bash}
# Sort with bcftools sort:
bcftools sort -O v -o lp_ll_introgression_LyCa_ref.sorted.vcf lp_ll_introgression_LyCa_ref.unsorted.vcf

mkdir lp_ll_introgression_intermediate_vcfs

for vcf in $(cat lp_ll_introgression_vcfs.list)
 do
  echo ${vcf}
  mv ${vcf} intermediate_vcfs
  mv ${vcf}.idx intermediate_vcfs
done
```

------------------------------------------------------------------------

## Standard Variant Filtering

The following variants were filtered from the VCF generated during the calling:

(1) Variants from Repetitive/Low mappability regions (defined by the reference genome annotation)
(2) Indels + Non-biallelic variants
(3) Non-variant SNPs (allele frequency = 1)
(4) and (5) Standard quality filters, as GATK standard practices

This was performed running a customized script (lp_ll_introgression_vcf_filters_1-5.sh) that uses a combination of bedtools, bcftools and gatk.

An additional custom script (summary_table_filters_1-5.sh) was then run to extract a table summarizing the filtering process, indicating how many variants were filtered (f_vars) and how many variants are left (e_vars) at each step.

------------------------------------------------------------------------

## Phasing variants

Phasing of variants will be conducted with a pipeline that first uses WhatsAp v.1.1 to create haplogroups (??) from individual read and population data. The output of WhatsAp is then passed to SHAPEIT v.4.2.1 that will infer the haplotypes of each sample for each chromosome.

All this is based on what Lorena already ran for the samples mapped to the *Felix catus* reference genome:
[Lorena - Phasing](https://github.com/lorenalorenzo/Phasing)

### Splitting the VCF

To divide my VCF into single population VCFs and further dividing those into single chromosome VCFs I ran a custom bash script (pop_chr_vcf_split.sh). The populations are defined as at the beginning of this md. The chromosomes I decided to keep are the larger ones, 18 autosomes and the X chromosome.

### Generate genetic map

To run SHAPEIT I also need to provide a genetic map for the SNPs to phase. As we don't have one, we will manually generate a genetic map by multiplying the physical distance in bp between SNPs and genome wide average recombination rate, which is 1.9 cM/Mbp. By cumulatively summing the multiplication of the physical distance from previous the SNP by 0.0000019, we obtain the cM value of each SNP. This approximation is not ideal but it's the only way we can provide a map. To calculate this I wrote a custom script (make_chr_gmap.sh) which will output a gmap table for each chromosome, made of 3 columns: position, chromosome, cM (format useful for SHAPEIT).
