
def coordinates(vcf, outdir):
    chrom_file = open(outdir + "/chrom_file.txt", "w")
    with open(vcf, "r") as openvcf: 
         for line in openvcf:
             if line.split("=")[0] == "##contig":
                 vcfcord = line.split("=")
                 vcf_scaff = vcfcord[2].split(",")[0]
                 vcf_scaff_len = int(vcfcord[3].split(">")[0])
                 if vcf_scaff_len > 100000:
                     chrom_file.write("%s\t%s\n"%(vcf_scaff, str(vcf_scaff_len)))
             if line[0] != "#":
                 break
    chrom_file.close()    
        
#################################################################################################