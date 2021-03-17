# B49_B75_B96_haplotypes

This repository contains Python programs used for comparative analyses of haplotypes between the maize inbred lines B49, B75, and B96 as reported in a study with the following title:

“Maize inbred line B96 is the source of large-effect loci for resistance to generalist but not specialist spider mites.”

The authors of this study are Huyen Bui, Robert Greenhalgh, Gunbharpur S. Gill, Meiyuan Ji, Andre H. Kurlovs, Christian Ronnow, Sarah Lee, Ricardo A. Ramirez, and Richard M. Clark. The manuscript is being submitted for journal review. A preprint is available on the bioRxiv preprint server at: 

https://doi.org/10.1101/2021.02.04.429847

Documentation on this GitHub repository will be updated with journal information once the manuscript is accepted and published.

Briefly, this repository has the code, written in Python 2.7, needed for the genome-wide haplotype analyses, and their visual display, as presented in Figure 6 and referenced in the manuscript's Results section.

Instructions for downloading the requisite dataset of variant calls for B49, B75, and B96, as assessed from alignments of high-throughput sequencing reads from each line to the maize B73 reference genome, as well as documentation/instructions for running the Python code, can be found in file “Pipeline_Commands_Documentation.txt”. For data downlaod, go to this figshare repository: https://doi.org/10.6084/m9.figshare.13708375.v1.

Additionally, python code (written for Python 3) to generate the input file used to plot Supplemental Figure S5 in a revision of the manuscript is provided (code file "Allele_freq_analysis_chr6_shared_haplotype.py"). The two input files for this program can be downloaded from the respective figshare data repository (see above; execute this program in a directory containing the two input files).
