# Master program that process genotypic data from file: ????

import argparse
import sys
import os
import pandas as pd
import numpy as np
import subprocess

PARSER = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter)

PARSER.add_argument("-v", "--vcf", required=True, help = "provide vcf file to transform SNP scaffold and position into new cord scheme")
PARSER.add_argument("-o", "--outdir", required=True, help = "provide directory for writing the vcf file")

PARSER.add_argument("-i", "--input", required = False, help = "input file provides cord transform order.  It has 5 columns with tabs: out_scaff in_scaff direction beginning end")
PARSER.add_argument("-m", "--matrix", required= False, action = "store_true", help = "USE THIS FLAG TO MAKE MATRIX")
PARSER.add_argument("-k", "--skip", required= False, action = "store_true", help ="USE THIS FLAG TO SKIP MAKING COVERAGE FILE")

PARSER.add_argument("-over", "--coverage_over", required=False, default=1.50,
                    help="Maximum coverage depth as a multiple of genome-wide average")
PARSER.add_argument("-under", "--coverage_under", required=False, default=0.25,
                    help="Minimum coverage depth as a multiple of genome-wide average")
PARSER.add_argument("-qds", "--qds", required=False, default=2,
                    help="Mininum variant score")
PARSER.add_argument("-sor", "--sor", required=False, default=3,
                    help="Maximum strand bias score")
PARSER.add_argument("-mq", "--mps", required=False, default=50,
                    help="Minimum mean mapping quality score")
PARSER.add_argument("-mqrs", "--mqrs", required=False, default=-8,
                    help="Minimum mean quality rank sum score")
PARSER.add_argument("-rprs", "--rprs", required=False, default=-8,
                    help="Minimum read pos rank sum")

PARSER.add_argument("-n", "--minora", required= False, default = 1, help ="minumim number of individuals that need to have the allele")
PARSER.add_argument("-t", "--thinner", required= False, default = 1, help ="not count alleles that are a certain number of bp away from the previous one")
PARSER.add_argument("-c", "--chunk", required= False, default = 0, help ="minimum number of alleles per window. automatically scales to nstrains*7 in a 100k window")

PARSER.add_argument("-e", "--equality", required = False, action = "store_true", help ="USE THIS FLAG TO MAKE PAIRWISE COMPARISON CSV FILE")
PARSER.add_argument("-w", "--window", required= False, default = 100000, help ="what is the window size?")
PARSER.add_argument("-s", "--slide", required= False, default = 10000, help ="what is the window size?")
PARSER.add_argument("-p", "--percentpass", required= False, default = 20, help ="what is %age of alleles that have to pass QC when making pairwise comparisons?")

PARSER.add_argument("-y", "--haplotypes", required = False, action = "store_true", help ="USE THIS FLAG TO MAKE HAPLOTYPE FILE")
PARSER.add_argument("-l", "--similarity", required= False, default = 99, help ="what is the minimum similarity to be considered part of the same haplotype?")

argies = PARSER.parse_args()

input_file = argies.input                                                                                                                                
vcf = argies.vcf
outdir = argies.outdir
outdir1 = outdir+"/info_files"
outdir2 = outdir+"/matrices"

matrix = argies.matrix
skip = argies.skip

coverage_over = float(argies.coverage_over)
coverage_under = float(argies.coverage_under)
qds = float(argies.qds)
sor = float(argies.sor)
mps = float(argies.mps)
mqrs = float(argies.mqrs)
rprs = float(argies.rprs)

chunk = int(argies.chunk)
minora = float(argies.minora)
thinner = float(argies.thinner)

equality = argies.equality
window = int(argies.window)
slide = int(argies.slide)
percent_pass_min = float(argies.percentpass)

haplotypes = argies.haplotypes
similarity = int(argies.similarity)                                

# This part actually gets the directory with all the files:

indir = os.path.dirname(os.path.realpath(__file__))+"/"


########################################################################### COODRDINATE PORTION ####################################################################################################
#Code will coordinate transform a VCF file based on an input file with cord transform information. 

coordinate_code = indir+"coordinate_functions.py"

scaffold_output_code = indir+"coordinate_output.py"

if not skip and not input_file:
    execfile(scaffold_output_code)
    print("*************************")
    print("running coordinate code !")
    coordinates(vcf, outdir)
    print("*************************")

if not skip and input_file:
    execfile(coordinate_code)
    # this runs the code
    print("*************************")
    print("running coordinate code !")
    coordinates(input_file, vcf, outdir)
    print("*************************")
    vcf = '%s/%s.%s.vcf' % (outdir, vcf.split(".")[0].split("/")[-1], "transformed")

############################################### now, let's get the actual transformed VCF file and the chromosome length information file! ##########################################################

sin = outdir+"/chrom_file.txt"
print(sin)

########################################################################### MATRIX PORTION ##########################################################################################################

matrix_code = indir+"matrix_code.py"
execfile(matrix_code)

if matrix:
    if not skip:
        # creates the coverage info file
        print("creating coverage info file")
        coverage(vcf, outdir1)
    coveragefile = outdir1+"/coverageinfo.txt"
    #print(coveragefile, "coveragefile")
    print("creating matrix")
    try:
        shader = scale_dict(sin)
    except IOError:
        print("YOU SHOULD NOT USE -k UNTIL YOU RUN COORDINATE TRANSFORM")
        sys.exit()
    try:
        get_matrix(vcf, shader, coveragefile, coverage_under, coverage_over, 
               qds, sor, mps, mqrs, rprs, minora, outdir2)
    except IOError:
        print("YOU SHOULD NOT USE -k WHEN THERE'S NO COVERAGE FILE")
        sys.exit()
    matrixfile = outdir1+"/"+"matrix_for_clustrering.txt"
else:
    if not skip:
        print("creating coverage info file")
        coverage(vcf, outdir1)
        coveragefile = outdir1+"/coverageinfo.txt"
    else:
        if not equality and not haplotypes:
            print("FIGURE OUT EXECUTION")
            sys.exit()

################################## again, let's generate files that come out fo this ##########################################################################################################################

infile = outdir2+"/matrix_for_clustering.txt"
failfile = outdir2+"/fail_data.txt"

############################################################################# CREATING SNP SIMILARITY CSV FILE  ################################################################################################

if equality:
    compare_code = indir+"compare_strains.py"
    execfile(compare_code)
    print("making pairwise comparison csv file")
    print("*************************")
    diff_calculator(infile, outdir, window, slide, sin, chunk, percent_pass_min, failfile)

comps = outdir + "/pairwise.csv"

############################################################################### RUNNING HAPLOTYPE PREDICTIONS ##################################################################################################

if haplotypes:
    haplo_code = indir+"shared_haplos.py"
    execfile(haplo_code)
    print("making haplotype file")
    print("*************************")
    find_haplos(comps, window, similarity, slide, sin, outdir)

################################################################################################################################################################################################################





