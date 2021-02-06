#!/usr/bin/env python
        
def error(message, line):                                                             
    '''error function will occur at anytime the input file is in the incorrect format'''                                                                                 
    sys.exit("%s\n%s" % (message, line))
    
def warning(message, line):                                                             
    '''error funtion for the scaffolds that are split up with parts left out; however it will not end the code it'll just give you a warning and continue'''
    sys.stderr.write("%s\n%s" % (message, line))

def coordinates(input, vcf, outdir):

    order_dict = {}                                                                         
    out_scaff_seq_dict = {}                                                                 
    out_scaff_seq_write_dict = {}
    out_scaff_order_list = []                                                                         
    in_scaff_list = []
    in_scaff_cord_dict = {}                                                                     
    vcf_input_info = {}
    vcf_dict = {}
    vcf_gap_dict = {}
    header_length_dict = {}
    header_gap_dict = {}
    header_gap_length_dict = {}
    header_uncalled_dict = {}

    fasta_dict = {}
    fasta_order_list= []
    run_over_input_file = False


    outfile_vcf = open('%s/%s.%s.vcf' % (outdir, vcf.split(".")[0].split("/")[-1], "transformed"), 'w')
    with open(vcf, "r") as openvcf: 
         for line in openvcf:
             if line[0:6]=="#CHROM":
                 header = line
                 outfile_vcf.write("%s" % (header))
             if line.split("=")[0] == "##contig":
                 vcfcord = line.split("=")
                 vcf_scaff = vcfcord[2].split(",")[0]
                 vcf_scaff_len = int(vcfcord[3].split(">")[0])
                 fasta_dict[vcf_scaff] = vcf_scaff_len
                 fasta_order_list.append(vcf_scaff)
             if line[0] != "#":
                 if run_over_input_file == False:
                    with open(input, "r") as input_file:                                                                                        
                        for line in input_file:                                                                 
                            if line[0] != "#":                                                                     
                                line = line.strip().split("\t")                                                 
                                if len(line) != 5:                                                                 
                                    error("Make 5 columns", "\t".join(line))
                                else:                                                                             
                                    out_scaff = line[0]                                                         
                                    in_scaff = line[1]                                                             
                                    direction = line[2]                                                         
                                    start = line[3]                                                             
                                    end = line[4]                                                                 
                                    if direction not in ["f", "r"]:                                             
                                        error("Letter in direction column", "\t".join(line))
                                    if start == "b":                                                             
                                        start = 1
                                    elif start.isdigit():                                                         
                                        start = int(start)
                                    else:                                                                         
                                        error("Problem with start", "\t".join(line))
                                    if end == "e":                                                                  
                                        end = fasta_dict[in_scaff]
                                    elif end.isdigit():
                                        end = int(end)
                                    else:                                                                         
                                        error("Problem with end", "\t".join(line))
                                if out_scaff not in order_dict:                                                 
                                    order_dict[out_scaff] = []                                                     
                                    out_scaff_seq_dict[out_scaff] = []                                             
                                    out_scaff_seq_write_dict[out_scaff] = []
                                    out_scaff_order_list.append(out_scaff)                                              
                                    vcf_dict[out_scaff] = {}
                                    header_length_dict[out_scaff] = {}
                                order_dict[out_scaff].append([in_scaff, direction, start, end])                    
                                if in_scaff not in in_scaff_cord_dict:                                             
                                    in_scaff_cord_dict[in_scaff] = []
                                    in_scaff_list.append(in_scaff)
                                    vcf_gap_dict[in_scaff] = {}
                                    header_gap_dict[in_scaff] = {}
                                    header_gap_length_dict[in_scaff] = 0
                                in_scaff_cord_dict[in_scaff].append([start, end])                                 
                    for in_scaff in in_scaff_cord_dict:
                        in_scaff_cord_dict[in_scaff].sort(key=lambda x: x[0])                                    
                        prevend = 0
                        for cord in in_scaff_cord_dict[in_scaff]:                                                
                            start = cord[0]                                                                        
                            if start <= prevend:                                                                
                                error("Raise for troubleshooting...", "%s %s %s %s %s" % (in_scaff, "start =", start, "previous end =", prevend))
                            if start > prevend + 1:
                                warning("parts of a scaffold/s missing", 
                                "%s%s\n" % ("broken scaffold=", in_scaff))
                            prevend = cord[1]
                        length = fasta_dict[in_scaff]
                        if prevend != length:
                            warning("parts of a scaffold/s missing", 
                                "%s%s\n" % ("broken scaffold=", in_scaff))
                    run_over_input_file = True
                    chrom_file = open(outdir + "/chrom_file.txt", "w")
                    for out_scaff in out_scaff_order_list:
                        length = 0
                        for in_scaff_info in order_dict[out_scaff]:
                            length_beg = in_scaff_info[2]
                            length_end = in_scaff_info[3]
                            length += ((length_end+1)-length_beg)
                        header_length = "".join(["##contig=<ID=", out_scaff, ",length=", str(length), ">"])
                        chrom_file.write("%s\t%s\n" %(out_scaff, str(length)))
                        outfile_vcf.write("%s\n" % (header_length))
                    chrom_file.close()
                    for in_scaff in in_scaff_cord_dict:    
                        prevend = 0
                        for cord in in_scaff_cord_dict[in_scaff]:
                            start = cord[0]    
                            if start > prevend + 1:
                                mis_start = prevend +1
                                mis_end = start 
                                header_gap_length_dict[in_scaff] += (mis_end - mis_start)
                            prevend = cord[1]
                        length = fasta_dict[in_scaff]
                        if prevend != length:
                            header_gap_length_dict[in_scaff] += ((length+1) - prevend)
                        if header_gap_length_dict[in_scaff] != 0 :
                            header_gap = "".join(["##contig=<ID=", in_scaff, ",length=", str(header_gap_length_dict[in_scaff]), ">"])
                            outfile_vcf.write("%s\n" % (header_gap))
                 
                 else:
                    SNP_line = line.strip().split("\t")
                    SNP_scaff = SNP_line[0]
                    SNP_pos = int(SNP_line[1])
                    for out_scaff in out_scaff_order_list:
                        for in_scaff_info in order_dict[out_scaff]:
                            input_scaff = in_scaff_info[0]
                            input_direction = in_scaff_info[1]
                            input_beg = in_scaff_info[2]
                            input_end = in_scaff_info[3]
                            if SNP_scaff == input_scaff:
                                if input_beg <= SNP_pos <= input_end:
                                    length = 0
                                    for in_scaff_info in order_dict[out_scaff]:
                                        length_scaff = in_scaff_info[0]
                                        length_direction = in_scaff_info[1]
                                        length_beg = in_scaff_info[2]
                                        length_end = in_scaff_info[3]
                                        if length_scaff == SNP_scaff:
                                            if length_beg <= SNP_pos <= length_end:
                                                if length_direction == "f":
                                                    SNP_transform = length + ((SNP_pos+1) - length_beg)
                                                else:
                                                    SNP_transform = length + ((length_end+1) - SNP_pos)
                                        length += ((length_end+1)-length_beg)
                                    SNP_line[0] = out_scaff
                                    SNP_line[1] = str(SNP_transform)
                                    transform_line = "\t".join(SNP_line)
                                    vcf_dict[out_scaff][SNP_transform] = transform_line
                    for in_scaff in in_scaff_cord_dict:
                        if in_scaff == SNP_scaff:
                            prevend = 0
                            for cord in in_scaff_cord_dict[in_scaff]:
                                start = cord[0]    
                                if start > prevend + 1:
                                    mis_start = prevend +1
                                    mis_end = start
                                    if mis_start < SNP_pos < mis_end:
                                        transform_line = "\t".join(SNP_line)
                                        vcf_gap_dict[SNP_scaff][SNP_pos] = transform_line
                                prevend = cord[1]
                            length = fasta_dict[in_scaff]
                            if prevend != length:
                                if prevend < SNP_pos < length:
                                    transform_line = "\t".join(SNP_line)    
                                    vcf_gap_dict[SNP_scaff][SNP_pos] = transform_line            
    print "printing file"
    for out_scaff in out_scaff_order_list:
        for positions in sorted(vcf_dict[out_scaff]):
            outfile_vcf.write("%s\n" % vcf_dict[out_scaff][positions])
    for scaffold in fasta_order_list:
        if scaffold in in_scaff_list:
            if vcf_gap_dict[scaffold] != []:
                for positions in sorted(vcf_gap_dict[scaffold]):
                    outfile_vcf.write("%s\n" % vcf_gap_dict[scaffold][positions])
    outfile_vcf.close()
        
        
        
#################################################################################################