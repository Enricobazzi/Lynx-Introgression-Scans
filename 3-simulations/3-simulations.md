---
title: "3-simulations"
author: "Enrico"
date: "5/20/2022"
output: html_document
---

In order to generate a training set for the deep learning pipeline of [FILET](https://github.com/kr-colab/FILET) ([Schrider et al. 2018](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007341)) we will simulate haplotype data under the demographic model we inferred using dadi (see [dadi.md](../2-infer_demographic_history/dadi/dadi.md) for details).

The software chosen for the coalescent simulations is [msmove](https://github.com/genevalab/msmove) ([Garridan & Geneva, 2014](http://dx.doi.org/10.6084/m9.figshare.1060474)), a modified version of Hudson's coalescent simulator [ms](http://home.uchicago.edu/~rhudson1/source/mksamples.html) ([Hudson 2002](https://academic.oup.com/bioinformatics/article/18/2/337/225783?login=true)).

Msmove offers the advantage of allowing for finer control and tracking of migrant genealogies, useful in our study.

example translation https://dadi.readthedocs.io/en/latest/examples/YRI_CEU/YRI_CEU/
build ms command https://dadi.readthedocs.io/en/latest/api/dadi/Misc.html#dadi.Misc.ms_command
