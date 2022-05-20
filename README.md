# Lynx-Introgression-Scans

In this repository you can find all of the workflows and scripts I wrote and ran in order to detect introgressed windows in the genomes of *Lynx pardinus* and *Lynx lynx* populations.

The repository is divided into subfolders, where all of the stages of the analysis are described and the respective scripts are stored.

The overall workflow was performed as follows:

1.  [Preparing the DataSet](./1-prepare_dataset/1-prepare_dataset.md)

    1.1.[Variant Calling using GATK](./1-prepare_dataset/1-prepare_dataset.md#./1-prepare_dataset/1-prepare_dataset.md#Variant%20Calling%20using%20GATK)

    1.2.[Standard Variant Filtering](./1-prepare_dataset/1-prepare_dataset.md#./1-prepare_dataset/1-prepare_dataset.md#Standard%20Variant%20Filtering)

    1.3.[Phasing variants](./1-prepare_dataset/1-prepare_dataset.md#./1-prepare_dataset/1-prepare_dataset.md#Phasing%20variants)

2.  [Inferring Demographic History](./2-infer_demographic_history/2-infer_demographic_history.md)

    2.1.[dadi](./2-infer_demographic_history/dadi/dadi.md)

3.  [Running Simulations](./3-simulations/3-simulations.md)
