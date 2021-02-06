
import sys

filey = sys.argv[1]
info_filey = sys.argv[2]
freq = int(sys.argv[3])
outfiley = sys.argv[4]

fileyh = open(filey)
outfileyh = open(outfiley, "w")


def scale_dict(info_file):
    info_read = open(info_file, "r")
    info_dict = {}
    final_lenny = 0
    lenny = 0
    chromlist = []
    for line in info_read:
        linetab = (line.rstrip()).split("\t")
        scaffy = linetab[0]
        final_lenny = final_lenny + lenny
        lenny = int(linetab[1])
        info_dict[scaffy] = final_lenny
        chromlist.append(scaffy)
    info_read.close()
    return(info_dict, chromlist)

scale_up = scale_dict(info_filey)[0]
chromlist = scale_dict(info_filey)[1]

current = 1
for line in fileyh:
    liney = line.split("\t")
    chroma = liney[0]
    locpos = liney[1]
    posse = liney[2]
    lastesies = len(liney[3:])
    if posse != "medpos":
        posse = int(posse)
        if posse >= current:
            for value in range(current, posse, freq):
                if scale_up[chroma] >= value:
                    print(value, scale_up[chroma])
                    where = chromlist.index(chroma)
                    current_chroma = chromlist[where-1]
                else:
                    current_chroma = chroma
                locposse = value - scale_up[current_chroma]
                outlist = [str(current_chroma), str(locposse), str(value)] + ["0"] * lastesies
                outstring = "\t".join(outlist) + "\n"
                outfileyh.write(outstring)
            outfileyh.write(line) 
            current = value+freq
        else:
            outfileyh.write(line)      
    else:
        outfileyh.write(line)

fileyh.close()
outfileyh.close()
            

            


