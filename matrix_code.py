# COVERAGE
def coverage(vcf, outdir1):
    cov = {}
    subprocess.call("mkdir %s"%(outdir1), shell=True)
    outcov = open(outdir1+"/coverageinfo.txt", "w")
    with open(vcf, "r") as openvcf:
        for line in openvcf:
            if line[0:2] == "##": #this contains the lengths of the scaffold for the original pass I don't care about it. 
                pass
            elif line[0:6] == "#CHROM": #need to store all the variable names for the samples and which columns correspond to it. 
                header = line.strip().split("\t")
                header_strains = header[9:]
                for strain in header_strains:
                    cov[strain] = []
            else:
                if line[0] != "#": #need average coverage for parent1, parent2, sel, and unsel
                    SNP_line = line.strip().split("\t")    
                    SNP_score = float(SNP_line[5])
                    SNP_scaff = SNP_line[0]
                    SNP_pos = SNP_line[1]
                    ref_allele = SNP_line[3]
                    alt_allele = SNP_line[4]
                    if len(ref_allele) == 1 and len(alt_allele) == 1:
                        haplotypes = SNP_line[9:]
                        index = 0
                        for strain in header_strains:
                            choice = haplotypes[index]
                            if ("./" in choice)==True:
                                index = index +1
                                continue
                            reads = choice.split(":")[1].split(",")
                            parent_cov = sum([int(i) for i in reads])
                            cov[strain].append(parent_cov)
                            index = index +1
    for strain in header_strains:
        outcov.write(strain + "\t" + str(int(sum(cov[strain])/len(cov[strain]))) + "\n")
    outcov.close()

def line_parser(line):
    """Parses VCF string to get mapping quality info"""
    linedict = {}
    line = line.split(";")
    for keyval in line:
        keyval = keyval.split("=")
        if len(keyval[1].split(",")) == 1:
            linedict[keyval[0]] = float(keyval[1])
    return(linedict)


# MATRIX
def get_matrix(vcf, shader, coveragefile, coverage_under, coverage_over, 
               qds, sor, mps, mqrs, rprs, minora, outdir2):
    cov = {}
    subprocess.call("mkdir %s"%(outdir2), shell=True)
    with open(coveragefile, "r") as opencov:
        for line in opencov:
            linestrip = (line.rstrip()).split("\t")
            st = linestrip[0]
            co = int(linestrip[1])
            cov[st] = co
    outfile = open(outdir2+"/matrix_for_clustering"+".txt", "w")
    outfile_fail = open(outdir2+"/fail_data"+".txt", "w")
    scaffs = []
    with open(vcf, "r") as openvcf:
        for line in openvcf:
            if line[0:2] == "##": #this contains the lengths of the scaffold for the original pass I don't care about it. 
                pass
            elif line[0:6] == "#CHROM": #need to store all the variable names for the samples and which columns correspond to it. 
                header = line.strip().split("\t")
                header_strains = header[9:]
                outfile.write("\t".join(["chrom", "pos", "medpos"]+header_strains)+"\n")
            else:
                if line[0] != "#":
                    SNP_line = line.strip().split("\t")    
                    SNP_score = float(SNP_line[5])
                    SNP_scaff = SNP_line[0]
                    SNP_pos = SNP_line[1]
                    int_pos = int(SNP_pos)
                    if (SNP_scaff in shader)== False:
                        continue
                    if (SNP_scaff in scaffs)==False:
                        scaffs.append(SNP_scaff)
                        old_pos = None
                    SNP_medpos = str(int(SNP_pos) + shader[SNP_scaff])
                    ref_allele = SNP_line[3]
                    alt_allele = SNP_line[4]
                    if len(ref_allele) == 1 and len(alt_allele) == 1:
                        SNP_metrics = line_parser(SNP_line[7])
                        if ("QD" in SNP_metrics 
                                and "MQ" in SNP_metrics 
                                and "SOR" in SNP_metrics
                                and "MQRankSum" in SNP_metrics 
                                and "ReadPosRankSum" in SNP_metrics
                                and SNP_metrics["QD"] >= qds
                                and SNP_metrics["MQ"] >= mps
                                and SNP_metrics["SOR"] < sor
                                and SNP_metrics["MQRankSum"] >= mqrs
                                and SNP_metrics["ReadPosRankSum"] >= rprs):
                            haplotypes = SNP_line[9:]
                            prescores = [SNP_scaff, SNP_pos, SNP_medpos]
                            scores = []
                            index = 0
                            for choice in haplotypes:
                                # is the variant there. 2 = not there
                                if ("./" in choice)==True:
                                    scores.append("F")
                                    index = index + 1
                                    continue
                                # is the coverage okay?
                                reads = choice.split(":")[1].split(",")
                                covey = sum([int(i) for i in reads])
                                strain = header_strains[index]
                                avecov = cov[strain]
                                if covey < avecov*coverage_under or covey > avecov*coverage_over:
                                    scores.append("F")
                                    index = index +1
                                    continue        
                                # now the genotype pass
                                genotype = choice.split(":")[0]
                                alleles = genotype.split("/")
                                allela = alleles[0]
                                all_all = len(alleles)
                                all_uniques = len(list(set(alleles)))
                                if all_uniques == 1:
                                     scores.append(allela)
                                     index = index +1
                                else:
                                     scores.append("F")
                                     index = index + 1
                             # accept alleles that pass qc in every strain
                            if ("F" in scores)==False: 
                                originals = list(set(scores))
                                firsties = originals[0]
                                if len(originals) <= 2:
                                    firsties = originals[0]
                                    totalsies = scores.count(firsties)
                                    if old_pos == None:
                                        if len(scores)-minora >=totalsies >= minora:
                                            outfile.write("\t".join(prescores+scores)+"\n")
                                            outfile_fail.write("\t".join(prescores+["PASS"])+"\n")
                                            old_pos = int_pos
                                    else:
                                        if len(scores)-minora >=totalsies >= minora and int_pos-old_pos > thinner:
                                            outfile.write("\t".join(prescores+scores)+"\n")
                                            outfile_fail.write("\t".join(prescores+["PASS"])+"\n")
                                            old_pos = int_pos
                            else:
                                outfile_fail.write("\t".join(prescores+["FAIL"])+"\n")
    outfile.close()
    outfile_fail.close()
                                
# this is to get the correct overall genomic position
def scale_dict(info_file):
    info_read = open(info_file, "r")
    info_dict = {}
    final_lenny = 0
    lenny = 0
    for line in info_read:
        linetab = (line.rstrip()).split("\t")
        scaffy = linetab[0]
        final_lenny = final_lenny + lenny
        lenny = int(linetab[1])
        info_dict[scaffy] = final_lenny
    info_read.close()
    return(info_dict)
            
#####################################################################################################################