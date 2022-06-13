import sys
import os
import numpy
import dadi
import pylab
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("pops", type=str,
                    help="population abbreviations of pair chosen")
args = parser.parse_args()

pops = args.pops

# vcf file
vcf_prefix = '/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression_LyCa_ref.sorted.filter5.phased.fixed.'
vcf_suffix = '.miss.rd_fil.vcf'
datafilename = [vcf_prefix, pops, vcf_suffix]

datafile = ''.join(datafilename)

# pop file
pop_prefix = '/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression/'
pop_suffix = '_dadi_popfile.txt'
popfilename = [pop_prefix, pops, pop_suffix]

popfile = ''.join(popfilename)

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict_vcf(datafile, popfile)

#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids = pops.split("-")

# projections to test in a loop

# range of LPA:
range1 = [30, 32, 34, 36, 38]
# range of WEL:
if pop_ids[1] == "wel":
    range2 = [32, 34, 36, 38, 40]
# range of EEL:
if pop_ids[1] == "eel":
    range2 = [30, 32, 34, 36, 38]
# range of SEL:
if pop_ids[1] == "sel":
    range2 = [16, 18, 20, 22, 24]

proj_list = []
fs_list = []

for i in range(0,5):
    print("testing projection in iberian: ", range1[i])
    proj_pop1 = range1[i]
    
    for n in range(0,5):
        print("with projection in eurasian: ", range2[n])
        proj_pop2 = range2[n]
        
        proj = [proj_pop1, proj_pop2]
        fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)
        fs_s = fs.S()
        proj_list.append(proj)
        fs_list.append(fs_s)
        print("number of SNPs: ", fs_s)

max_fs_s = max(fs_list)
max_fs_index = fs_list.index(max_fs_s)
max_fs_proj = proj_list[max_fs_index]

print("best fs value is : ", max_fs_s, " given by projecting : ", max_fs_proj)
