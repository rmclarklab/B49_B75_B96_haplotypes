
def setter(s, adj): # basically tells you where you have shared stuff
	bella = set(s.split(";"))
	for stuff in adj:
		news = set(stuff.split(";"))
		if len(bella & news) > 0:
			bella = bella | news
	return(bella)

# pass 2
def setter2(s, adj): # basically tells you where you have shared stuff
	bella = s
	for stuff in adj:
		if len(bella & stuff) > 0:
			bella = bella | stuff
	return(bella)

# this is to get the correct overall genomic position
def scale_dict(info_file):
	info_read = open(info_file, "r")
	info_dict = {}
	for line in info_read:
		linetab = (line.rstrip()).split("\t")
		scaffy = linetab[0]
		lenny = int(linetab[1])
		info_dict[scaffy] = lenny
	info_read.close()
	return(info_dict)
	
def find_haplos(comps, window, minimum, slide, chromfile, outdir):
    outfile1 = open(outdir + "/haplotypes_chrom.txt", "w")
    outfile2 = open(outdir + "/haplotypes_absolute.txt", "w")
    outfile3 = open(outdir + "/pairwise_comparisons.txt", "w")

    compy = pd.read_csv(comps)
    all_strains = set()
    strainies = list(compy.columns[3:])
    for compppp in strainies:
        both = set(compppp.split(";"))
        all_strains = all_strains | both
    total_pairwise_possible = float(sum(range(len(all_strains))))
    
    
    winny = window/2

    dicta = scale_dict(chromfile)

    chroms = list(compy['scaff'].unique())

    outfile1.write("\t".join(["chrom","beg", "end", "nhaplos", "nstrains", "haplos"])+"\n")
    outfile2.write("\t".join(["chrom","beg", "end", "nhaplos", "nstrains", "haplos"])+"\n") 

    for chrom in chroms: # process each chrom separately
        compy_adj = compy[compy['scaff'] == chrom]
        totals = len(compy_adj)
        nkey = 1
        for i in range(totals):
            meddy = compy_adj.iloc[i,1]
            totalmed = compy_adj.iloc[i,2]
            beg = meddy - winny
            end = meddy + winny
            if beg == 0:
                beg = 1
            if end > dicta[str(chrom)]:
                end = dicta[str(chrom)]
            overbeg = totalmed - winny
            overend = totalmed + winny
            stria = compy_adj.iloc[i,3:]
            higher = stria[stria > minimum] # only select the ones that are of high similarity
            higher_list = list(higher.index)
            higher_ser = pd.Series(higher_list)
            unique_sets = []
            # pass 1
            sets = list(higher_ser.apply(setter, args = (higher_list,)))
            sets_series = pd.Series(sets)
            # pass 2
            setlist = list(sets_series.apply(setter2, args = (sets,)))
            totalsets = 0
            totalstrains = 0
            total_pairwise = 0
            for setty in setlist:
                setlist = sorted(list(setty))
                if (setlist in unique_sets) == False:
                    totalsets = totalsets+1
                    nstr = len(setty)
                    totalstrains = totalstrains+nstr
                    total_pairwise = total_pairwise + sum(range(len(setty)))
                    unique_sets.append(setlist)
            outfile1.write("\t".join([str(chrom), str(beg), str(end), str(totalsets), str(totalstrains), str(unique_sets)[1:-1]])+"\n")
            outfile2.write("\t".join([str(overbeg), str(overend), str(totalsets), str(totalstrains), str(unique_sets)[1:-1]])+"\n")
            outfile3.write("\t".join([str((overbeg+overend)/2), str(total_pairwise/total_pairwise_possible)])+"\n")

    outfile1.close()
    outfile2.close()	
	

				
						
						
						
						
					
		
				
				
	
