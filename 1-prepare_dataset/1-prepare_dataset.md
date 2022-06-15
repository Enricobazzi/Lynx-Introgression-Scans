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

From our previous work we have already generated BAM files for each of these samples (yay!), generated using the Canada lynx assembly as the reference genome (mLynCan4_v1.p; GCA_007474595.1; [Rhie et al., 2021](http://www.nature.com/articles/s41586-021-03451-0)).

To generate a VCF file from these BAMs, we performed variant calling on each scaffold in parallel, while additionally dividing all of the bigger scaffolds (18 autosomic + X) into 80 chunks each. This way more jobs have to be sbatched, but we have no memory or time problems.

We called scaffolds that could be called without being divided into chunks using the script [perchr_haplotypecaller.sh](./perchr_haplotypecaller.sh)

```{bash}
# short non-autosomic (??) scaffolds 
chr_list=($(cat /mnt/netapp1/Store_csebdjgl/reference_genomes/Ref_Genome_LyCa/non-autosomic_scaffolds_list.txt | grep -vE "mLynCan4_MT|Super_Scaffold_10"))
for chr in ${chr_list[@]}
 do
  echo "sbatching perchr_haplotypecaller for $chr"
  sbatch --mem=120GB perchr_haplotypecaller.sh $chr
done
```

Chunk divided scaffolds were then called using the script [perchunk_haplotypecaller.sh](./perchunk_haplotypecaller.sh)

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

(1) Variants from Repetitive/Low mappability regions (defined by the reference genome annotation), (2) Indels + Non-biallelic variants, (3) Non-variant SNPs (allele frequency = 1), (4) and (5) Standard quality filters, as GATK standard practices

This was performed running a customized [script](./lp_ll_introgression_vcf_filters_1-5.sh) that uses a combination of bedtools, bcftools and gatk.

An additional custom [script](./summary_table_filters_1-5.sh) was then run to extract a table summarizing the filtering process, indicating how many variants were filtered (f_vars) and how many variants are left (e_vars) at each step.

This command was run to generate a BED of SNPs filtered because of low quality

```{bash}
bedtools subtract \
 -a /GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression_LyCa_ref.sorted.filter3.vcf \
 -b /GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression_LyCa_ref.sorted.filter5.vcf |
 awk '{print $1, $2-1, $2}' | tr ' ' '\t' \
 > /GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression/filter_beds/qual_filter.bed
``` 

------------------------------------------------------------------------

## Phasing variants

Phasing of variants will be conducted with a pipeline that first uses [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html) v.1.1 ([Martin et al., 2016](https://www.biorxiv.org/content/10.1101/085050v2)) to create phase sets from individual read and population data. The output of WhatsHap is then passed to [SHAPEIT4](https://odelaneau.github.io/shapeit4/) v.4.2.1 ([Delaneau et al., 2019](https://www.nature.com/articles/s41467-019-13225-y)) that will infer the haplotypes of each sample for each chromosome.

All this is based on what Lorena already ran for the samples mapped to the *Felix catus* reference genome:

[Lorena - Phasing](https://github.com/lorenalorenzo/Phasing)

### Splitting the VCF

To divide my VCF into single population VCFs and further dividing those into single chromosome VCFs I ran a custom bash [script](./pop_chr_vcf_split.sh). The populations are defined as at the beginning of this md. The chromosomes I decided to keep are the larger ones, 18 autosomes and the X chromosome.

### Generate genetic map

To run SHAPEIT4 I also need to provide a genetic map for the SNPs to phase. As we don't have one, we will manually generate a genetic map by multiplying the physical distance in bp between SNPs and genome wide average recombination rate, which is 1.9 cM/Mbp. By cumulatively summing the multiplication of the physical distance from previous the SNP by 0.0000019, we obtain the cM value of each SNP. This approximation is not ideal but it's the only way we can provide a map. To calculate this I wrote a custom [script](./make_chr_gmap.sh) which will output a gmap table for each chromosome, made of 3 columns: position, chromosome, cM (format useful for SHAPEIT4).

### Generate Phase sets with WhatsHap

For more precise phasing, we first run the software WhatsHap using the --tag=PS (see [link](https://whatshap.readthedocs.io/en/latest/guide.html#representation-of-phasing-information-in-vcfs)).

Phase sets were generated from the VCF of each chromosome of each population by running in parallel a custom [script](./pop_chr_vcf_whatshap.sh).

```{bash}
pop_list=($(cat /mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression/lp_ll_introgression_populations.txt | cut -f2 | sort -u))
chr_list=($(cat /mnt/netapp1/Store_csebdjgl/reference_genomes/Ref_Genome_LyCa/big_scaffolds.bed | cut -f1))
for pop in ${pop_list[@]}
 do
  for chr in ${chr_list[@]}
   do
    echo "sbatching whatshap phasing of ${chr} VCF of ${pop}"
    sbatch pop_chr_vcf_whatshap.sh ${pop} ${chr}
  done
done
```

*note for the future* : time to run 19-20 samples ranged from \~1 to \~7 hours - adjust sbatch time for faster queues.

### Phase using SHAPEIT4

Because SHAPEIT requires a minimum of 20 samples in order to phase a VCF, we will have to trick it by duplicating the genotypes of our samples. To do this I wrote a custom [script](./duplicate_pop_chr_vcf.sh) based on what [Lorena](https://github.com/lorenalorenzo/Phasing/blob/main/duplicating_for_phasing.sh) did, that will duplicate the samples of the output of the WhatsHap Phase-Set VCF.

I run it for all of the populations and chromosomes.

```{bash}
pop_list=($(cat /mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression/lp_ll_introgression_populations.txt | cut -f2 | sort -u))
chr_list=($(cat /mnt/netapp1/Store_csebdjgl/reference_genomes/Ref_Genome_LyCa/big_scaffolds.bed | cut -f1))
for pop in ${pop_list[@]}
 do
  for chr in ${chr_list[@]}
   do
    echo "duplicating whatshap ${pop}'s phase set VCF of ${chr}"
    ./duplicate_pop_chr_vcf.sh ${pop} ${chr}
  done
done
```

I then need to zip and index each duplicated phase set VCF.

```{bash}
pop_list=($(cat /mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression/lp_ll_introgression_populations.txt | cut -f2 | sort -u))
chr_list=($(cat /mnt/netapp1/Store_csebdjgl/reference_genomes/Ref_Genome_LyCa/big_scaffolds.bed | cut -f1))
for pop in ${pop_list[@]}
 do
  for chr in ${chr_list[@]}
   do
    echo "zipping ${pop}'s duplicated phase set VCF of ${chr}"
    bgzip lp_ll_introgression_filtered_${pop}_${chr}_ps_duplicate.vcf
    echo "indexing ${pop}'s zipped duplicated phase set VCF of ${chr}"
    bcftools index lp_ll_introgression_filtered_${pop}_${chr}_ps_duplicate.vcf.gz
  done
done
```

The data is now ready to be phased using SHAPEIT4. To do so in parallel, I used a custom made [script](./pop_chr_vcf_shapeit.sh) that runs SHAPEIT4 for each population and chromosome combinations.
MCMC iterations were set to "10b,1p,1b,1p,1b,1p,1b,1p,10m" as suggested by the SHAPEIT4 manual.

```{bash}
pop_list=($(cat /mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression/lp_ll_introgression_populations.txt | cut -f2 | sort -u))
chr_list=($(cat /mnt/netapp1/Store_csebdjgl/reference_genomes/Ref_Genome_LyCa/big_scaffolds.bed | cut -f1))
for pop in ${pop_list[@]}
 do
  for chr in ${chr_list[@]}
   do
   echo "sbatching SHAPEIT4 phasing of ${chr} VCF of ${pop}"
   sbatch pop_chr_vcf_shapeit.sh ${pop} ${chr}
  done
done
```

*note for the future* : run time is extremely short, between 2 and 15 minutes! - adjust sbatch time for faster queues!

To remove the duplicated samples from the phased VCF I use a custom made [script](./gt_masker_pop_chr_vcf.sh) that will also remove any GT that was imputed by SHAPEIT4, that I would like to keep as missing data.

```{bash}
pop_list=($(cat /mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression/lp_ll_introgression_populations.txt | cut -f2 | sort -u))
chr_list=($(cat /mnt/netapp1/Store_csebdjgl/reference_genomes/Ref_Genome_LyCa/big_scaffolds.bed | cut -f1))
for pop in ${pop_list[@]}
 do
  for chr in ${chr_list[@]}
   do
    echo "removing duplicates and un-imputing GTs from ${pop}'s phased VCF of ${chr}"
    ./gt_masker_pop_chr_vcf.sh ${pop} ${chr}
  done
done
```

The different single chromosome population VCFs for each populations were combined into a single VCF with data from all of the samples and all of the chromosomes using the software [bcftools merge](https://vcftools.github.io/htslib.html#merge) ([Li 2011](https://academic.oup.com/bioinformatics/article/27/21/2987/217423?login=true)) and [vcftools concat](https://vcftools.github.io/perl_module.html#vcf-concat) ([Danecek et al. 2011](https://academic.oup.com/bioinformatics/article/27/15/2156/402296)) in a custom [script](./combine_phased_vcfs.sh).

------------------------------------------------------------------------

## Additional Filtering

Two additional filtering are applied to the data.

One is based on missing genotype information and is implemented to avoid analyzing variants that are not informative enough and might introduce noise into our results.

The other is based on read depth and is implemented to avoid including in the analysis possible paralogs whose SNP profiles do not reflect real genetic diversity.

Both these filters are population specific. We are aiming to generate VCFs with data from only a pair of populations, in order to identify windows introgressed from one population into the other. This means the SNPs to be filtered out should be calculated for each population independently and then applied only if the population is included in the VCF for the analysis. Three distinct VCFs will be generated, for each population pair we are aiming for, which are the three Eurasian lynx populations always paired with the Iberian lynx one. 

### Calculating filter based on Missing Data

I calculated the number of missing genotypes in each population for each SNP in order to draw a distribution of data missingness across the entire genome. Using a bash [script](./pop_vcf_missing_gts.sh) we generate a table. The table can be then read by a R [script](./filter_missing_data.R) that will output a resume table and a graph that help visualize how different limits on missing data affect the number of SNPs filtered.

The final decision, based on these results, is to filter out any SNP with **15%** or more missing data, which results in a loss of ~5.4% of SNPs in Lynx pardinus, ~1.5% in Western and Eastern Eurasian lynx and <1% in the Southern Eurasian lynx. Higher missing rate in Lynx pardinus is probably given by the combination of the older sequencing technologies used for some of the samples and the relatively low depth of the rest of the samples.

### Calculating filter based on Read Depth

Mean read depth in consecutive 10kb windows along the genome was calculated using the software [samtools](http://www.htslib.org/doc/samtools.html) ([Li et al. 2009](https://academic.oup.com/bioinformatics/article/25/16/2078/204688)) in a custom [script](./pop_depth_10kwin.sh).

```{bash}
pop_list=($(cat /mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression/lp_ll_introgression_populations.txt | cut -f2 | sort -u))
for pop in ${pop_list[@]}
 do
  sbatch -t 02-00:00 -n 1 -c 1 --mem=90GB pop_depth_10kwin.sh ${pop}
done
```

The script outputs a bed file for each chromosome of each population, reporting the mean read depth for each 10kb window. These bed files can be analysed using an R [script](./filter_rd.R), that calculates and plots the distribution of mean read depth values for all windows. Based on these distributions, the same [script](./filter_rd.R) will also calculate a maximum depth value and output a bed file for each population with the windows to be excluded based on this maximum depth.

The maximum depth is calculated as the overall mean read depth + 0.5 times the standard deviation of mean read depth values. This ends up excluding ~1100 windows in each population, which correspond to around 0.5% of all windows.

### Applying the filters

Before we can divide the phased VCF into the three desired population-pair VCFs, we need to quickly adapt the phased VCF's header, that saw most of the important information stripped away during phasing (will give problems in GATK if we skip this). To "fix" the header we simply copy the pre-phased header, add the few new fields added during phasing, and finally add the phased part of the table:

```{bash}
# take pre-phased header
grep "##" lp_ll_introgression_LyCa_ref.sorted.filter5.vcf \
 > lp_ll_introgression_LyCa_ref.sorted.filter5.phased.fixed.vcf

# add info and format added in phasing
grep -E "##INFO|##FORMAT" lp_ll_introgression_LyCa_ref.sorted.filter5.phased.vcf \
 >> lp_ll_introgression_LyCa_ref.sorted.filter5.phased.fixed.vcf

# add phased vcf table
grep -v "##" lp_ll_introgression_LyCa_ref.sorted.filter5.phased.vcf \
 >> lp_ll_introgression_LyCa_ref.sorted.filter5.phased.fixed.vcf
```

After this we can use a custom [script](./split_miss_rd_filter.sh) to split the VCF into the three population-pair VCFs and apply the specific missing data and read depth filters.

## Removing Genes for Demographic Inference

Demographic inferences are based on neutral genomic history and it's therefore safer to remove all of the genes from the data when running, to avoid any biases that selection might be introducing in allele frequency trajectories across time.

We have a list of genomic coordinates for gene and their surroundings (5kb down- and up-stream) that was generated during our previous work ([Bazzicalupo et al., 2022](https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13439), [Lucena et al., 2020](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15366)). We filter these from the VCF as follows:

```{bash}
for pop in wel eel sel
 do
  prefix=/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression_LyCa_ref.sorted.filter5.phased.fixed.lpa-${pop}.miss.rd_fil
  echo "removing genes from lpa-${pop} population pair VCF"
  bedtools subtract -header -a ${prefix}.vcf \
   -b /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.genes.plus5000.bed \
   > ${prefix}.intergenic.vcf
done
```

## Deciding window size for detecting Introgression

In order to make some type of inference regarding a window being introgressed from one lineage into another, we need sufficient information to be contained in that window. We aim to have an average of around 200 SNPs per window in order to make our inferences solid.

We explore the number of SNPs present in 10kb and 100kb windows. Bed files of consecutive windows of these sizes along the genome are generated for the autosomic scaffolds using the following bedtools command:

```{bash}
# 10kb - 228047
grep -v "Super_Scaffold_10" /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.bed |
 bedtools makewindows -b stdin -w 10000 > \
 /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.10kb_wins.noX.bed

# 100kb - there are 22813
grep -v "Super_Scaffold_10" /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.bed |
 bedtools makewindows -b stdin -w 100000 > \
 /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.100kb_wins.noX.bed
```

Before proceeding with counting the number of SNPs in each window for the two window sizes, we need to check the amount of 10kb windows that are to be excluded for excessive read depth (see Additional Filtering above):

```{bash}
# check number of 100kb windows that overlap with at least one 10kb window 
# that has been filtered out (0 SNPs) because of read depth
bedtools intersect \
 -a /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.100kb_wins.noX.bed \
 -b /GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression/filter_beds/allpops_rd_filter.bed \
 -c | awk '{FS="\t"; OFS="\t"; if ($4>0) print $0}' | wc -l
```

There are a total of 910 of the 22813 100kb windows that overlap the regions filtered because of read depth. If we exclude those 100kb windows we lose around 4% of the genome.

```{bash}
# extract non-overlapping windows
bedtools intersect \
 -a /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.100kb_wins.noX.bed \
 -b /GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression/filter_beds/allpops_rd_filter.bed \
 -c | awk '{FS="\t"; OFS="\t"; if ($4==0) print $0}' | cut -f1-3 \
 > /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.100kb_wins.noX.no_rd.bed
```

Once filtered, we calculate mean number of SNPs in the 10kb and 100kb windows

```{bash}
for pop in wel eel sel
 do
  vcf=/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression_LyCa_ref.sorted.filter5.phased.fixed.lpa-${pop}.miss.rd_fil.vcf
  echo "lpa-${pop} pair 10kb:" 
  # 10kb 
  bedtools intersect \
   -a /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.10kb_wins.noX.bed \
   -b $vcf \
   -c | awk -v N=4 '{ sum += $N } END { print sum / NR }'

  echo "lpa-${pop} pair 100kb:" 
  # 100kb 
  bedtools intersect \
 -a /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.100kb_wins.noX.no_rd.bed \
 -b $vcf \
 -c | awk -v N=4 '{ sum += $N } END { print sum / NR }'
done
```

The results are:

lpa-wel pair 10kb: 25.9727
lpa-wel pair 100kb: 261.661
lpa-eel pair 10kb: 25.9331
lpa-eel pair 100kb: 261.357
lpa-sel pair 10kb: 26.0593
lpa-sel pair 100kb: 262.759

The window size of 100kb, although offering less resolution on a genome-wide scale and loosing additional genome due to read depth, seems better suited for our needs, as the average number of reliable SNPs in 100kb windows is >250.