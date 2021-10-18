# B49_B75_B96_haplotypes

This repository contains Python programs used for comparative analyses of haplotypes between the maize inbred lines B49, B75, and B96 as reported in a study with the following title:

“Maize Inbred Line B96 Is the Source of Large-Effect Loci for Resistance to Generalist but Not Specialist Spider Mites.”

The authors of this study are Huyen Bui, Robert Greenhalgh, Gunbharpur S. Gill, Meiyuan Ji, Andre H. Kurlovs, Christian Ronnow, Sarah Lee, Ricardo A. Ramirez, and Richard M. Clark. The manuscript has been published in the journal "Frontiers in Plant Science," is open access, and can be freely downloaded here:

https://doi.org/10.3389/fpls.2021.693088

Briefly, this repository has the code, written in Python 2.7, needed for the genome-wide haplotype analyses, and their visual display, as presented in Figure 6 and referenced in the manuscript's Results section.

Instructions for downloading the requisite dataset of variant calls for B49, B75, and B96, as assessed from alignments of high-throughput sequencing reads from each line to the maize B73 reference genome, as well as documentation/instructions for running the Python code, can be found in file “Pipeline_Commands_Documentation.txt”. For data downlaod, go to this figshare repository: https://doi.org/10.6084/m9.figshare.13708375.v1.

Additionally, python code (written for Python 3) to generate the input file used to plot Supplemental Figure S5 in the manuscript is provided (code file "Allele_freq_analysis_chr6_shared_haplotype.py"). The two input files for this program can be downloaded from the respective figshare data repository (see above; execute this program in a directory containing the two input files).

The development and release of code hosted on this Github repository, as well as the datasets that have been made available via the figshare repository, were funded in whole or in part by a US National Science Foundation award (1444449, Plant Genome Research Program) to Richard M. Clark (University of Utah) and Ricardo Ramirez (Utah State University).
