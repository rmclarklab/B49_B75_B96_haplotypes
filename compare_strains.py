
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

def diff_calculator(infile, outdir, win, slide, scale_file, chunk, percentpass, failfile):
    scale = scale_dict(scale_file)
    fails = pd.read_table(failfile, header = None, low_memory=False)
    snps = pd.read_table(infile, sep="\t", low_memory=False)
    filey = outdir+"/pairwise.csv"
    outfail = open(outdir+"/window_params.txt", "w")
    outfile = open(filey, "w")
    scaffolds = snps['chrom'].unique()
    nstrains = len(snps.columns[3:])
    if not chunk:
        minchunk = nstrains*7
    else:
        minchunk = chunk
    header = ["scaff","medpos", "total_med"]
    head = False
    for scaff in scaffolds:
        scaff_seg = snps[snps['chrom']==scaff]
        fail_seg = fails[fails[0]==scaff]
        beg = 0
        end = win
        dead = True
        while end < max(scaff_seg['pos']):
            chunk = scaff_seg[(scaff_seg['pos'] > beg) & (scaff_seg['pos'] < end)]
            failseg = fail_seg[(fail_seg[1] > beg) & (fail_seg[1] < end)]
            med = beg + win/2
            total_med = scale[str(scaff)]+med
            if len(failseg) > 0:
                fragfail = len(failseg[failseg[3] == "FAIL"])/float(len(failseg))
            else:
                fragfail = 1
            outfail.write("\t".join([str(total_med), str(len(chunk)), str(fragfail)])+"\n")
            if len(chunk) >= minchunk and 100-fragfail*100 >= percentpass:
                dead = False
                out = [str(scaff), str(med), str(total_med)]
                for col in chunk.columns[3:]:
                    outline = [col]
                    section = chunk[col]
                    for colcomp in chunk.columns[3:]:
                        # goes over all the columns and compared them
                        if col == colcomp:
                            continue
                        compy = sorted([col,colcomp])
                        compy_right = str(compy[0]+";"+compy[1])
                        if compy_right not in header:
                            header.append(compy_right)
                            compection = chunk[colcomp]
                            compare = (section == compection)
                            simindex = len(compare[compare==True])/float(len(compare))
                            simindex = round(simindex, 5)*100.00
                            out.append(str(simindex))
                        else:
                            continue
                if head == False:
                    outfile.write(",".join(header)+"\n")
                    compy_totals = len(header)-3
                    head = True
                outfile.write(",".join(out)+"\n")
                header = ["scaff","medpos", "total_med"]
            else:
                if dead == False and head == True:
                    outfile.write(",".join([str(scaff), str(med), str(total_med)] + ["-1"]*compy_totals)+"\n")
                    dead = True
            beg = beg + slide
            end = end + slide
    outfile.close()
    outfail.close()
    