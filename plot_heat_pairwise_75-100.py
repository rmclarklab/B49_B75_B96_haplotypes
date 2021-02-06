from numpy import mean
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plot
import matplotlib.mlab as mlab
import matplotlib.pylab as lab
import matplotlib.patches as patches
import matplotlib.ticker as plticker
from matplotlib import rcParams
from matplotlib import gridspec
from matplotlib import cm
import sys
import time

rcParams['font.sans-serif'] = 'Arial'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams["axes.edgecolor"] = "0.15"
rcParams["axes.linewidth"]  = 1

infile = sys.argv[1] # which file to go over
window = int(sys.argv[2]) # what size window was used?
slide = int(sys.argv[3]) # what slide
scaffold_file = sys.argv[4] # file with the scaffold and length
color_scheme = sys.argv[5] # what color scheme to use?
outdir = sys.argv[6] # where to output the pdf

# this is to process the file and output the shading correctly
def scale_dict(info_file):
    info_read = open(info_file, "r")
    info_dict = {}
    final_lenny = 0
    lenny = 0
    lista = []
    for line in info_read:
        linetab = (line.rstrip()).split("\t")
        scaffy = linetab[0]
        final_lenny = final_lenny + lenny
        lenny = int(linetab[1])
        info_dict[scaffy] = final_lenny
        lista.append(scaffy)
    info_read.close()
    return(info_dict, lista)

# processing the scaffold information
oyoy = scale_dict(scaffold_file)
shader = oyoy[0]
lista = oyoy[1]

midwin = window/2

outfiley = open(outdir+"/connected_fragments.txt", "w")

b = pd.read_csv(infile)
maxx = max(b["total_med"])
comparisons = list(b.columns[3:])

strains = []

ncomps = len(comparisons)

for comp in comparisons:
    strainies = comp.split(";")
    if (strainies[0] in strains)==False:
        strains.append(strainies[0])
    if comp == comparisons[ncomps-1]:
        strains.append(strainies[1])


nstrains = len(strains)
print("nstrains", nstrains)

grids = 8
print("grids", grids)

# this needs to change
fig = plot.figure(figsize=(6.85, 4), dpi=1200)

gs = gridspec.GridSpec(grids, grids)

begger = 1
mover = nstrains-1
ender = begger + mover

lastfirst = None

boxcount = pd.Series()

## let's do colors
## have to normalize the data set first so that they're continuous

norm = matplotlib.colors.Normalize(vmin=60, vmax=100)
m = cm.ScalarMappable(norm=norm, cmap=color_scheme)


for comparison in comparisons:
    print("%s: %s"%("processing", comparison))
    compare = comparison.split(";")
    compcol = b[comparison]
    first = compare[0]
    if lastfirst == None:
        lastfirst = first
        ax = plot.subplot(gs[begger:ender, :])
        standing = strains.index(first)
        remains = strains[standing+1:][::-1]
        ncomp = len(remains)-1
        tickies = list(pd.Series(range(len(remains)))+0.5)
        # label offset needs to be .2 to worky well
        plotsizies = ender - begger
        tickplace = (plotsizies+0.25)/plotsizies
        ax.set_title(first, y=tickplace)
        ####
        ax.set_ylim(0,ender-begger)
        ax.set_xlim(0, max(b['total_med']))
        yloc = plticker.FixedLocator(tickies)
        ax.yaxis.set_ticklabels(remains)
        ax.yaxis.set_major_locator(yloc)
        xloc = []
        prevval = None
        for scaff in lista:
            if prevval == None:
                prevval = shader[scaff]
            else:
                midscaff = prevval + int((shader[scaff]-prevval)/2)
                xloc.append(midscaff)
                prevval = shader[scaff]
        midscaff = prevval + int((max(b['total_med'])-prevval)/2)
        xloc.append(midscaff)
        plot.xticks(xloc, lista)
        lockyx = list(ax.xaxis.get_major_locator().tick_values(0,max(b['total_med'])))
        # this are the variables for the shading below
        old = None
        shade = False
        # here comes the shading
        for contig in lista:
            val = shader[contig]
            if old != None and shade == True:
                plot.axvspan(old, val, color='0.80', alpha=0.5, linewidth=0)
                shade = False
            else:
                if old != None:
                    shade = True
            old = shader[contig]
        # the last one
        if shade == True:
            plot.axvspan(old, maxx, color='0.80', alpha=0.5, linewidth=0)
    else:
        if first == lastfirst:
            ncomp = ncomp- 1
            pass
        else:
            standing = strains.index(first)
            remains = strains[standing+1:][::-1]
            ncomp = len(remains)-1
            begger = ender + 4
            mover = mover - 1
            ender = begger + mover
            lastfirst = first
            ax = plot.subplot(gs[begger:ender, :])
            tickies = list(pd.Series(range(len(remains)))+0.5)
            # label offset needs to be .2 to work well
            plotsizies = ender - begger
            tickplace = (plotsizies+0.25)/plotsizies
            ax.set_title(first, y=tickplace)
            ####
            ax.set_ylim(0,ender-begger)
            ax.set_xlim(0, max(b['total_med']))
            yloc = plticker.FixedLocator(tickies)
            ax.yaxis.set_ticklabels(remains)
            ax.yaxis.set_major_locator(yloc)
            xloc = []
            prevval = None
            for scaff in lista:
                if prevval == None:
                    prevval = shader[scaff]
                else:
                    midscaff = prevval + int((shader[scaff]-prevval)/2)
                    xloc.append(midscaff)
                    prevval = shader[scaff]
            midscaff = prevval + int((max(b['total_med'])-prevval)/2)
            xloc.append(midscaff)
            plot.xticks(xloc, lista)
            # this are the variables for the shading below
            old = None
            shade = False
            # here comes the shading
            for contig in lista:
                val = shader[contig]
                if old != None and shade == True:
                    plot.axvspan(old, val, color='0.80', alpha=0.5, linewidth=0)
                    shade = False
                else:
                    if old != None:
                        shade = True
                old = shader[contig]
            # the last one    
            if shade == True:
                plot.axvspan(old, maxx, color='0.80', alpha=0.5, linewidth=0)
    second = compare[1]
    ############################## 75 ####################################
    asel = compcol[(compcol >= 75) & (compcol < 83)]
    samies = list(asel.index)
    singletons = set(samies)
    if len(samies) > 0:
        print("processing 75-83% similar")
        backsies = pd.Series(samies[1:]+[-20])
        fronties = pd.Series([-20]+samies[:-1])
        normies = pd.Series(samies)
        topgo = backsies-samies
        bottomgo = samies-fronties
        tops = pd.Series([-20]+list(backsies[topgo == 1]))
        bottoms = pd.Series(list(fronties[bottomgo == 1]))
        cancels = tops-bottoms
        endies = list(tops[cancels != 0])[1:] # no need to +1 because end included with .loc
        beggies = list(bottoms[cancels != 0])        
        for terrier in range(len(endies)):
            slicey = b.loc[beggies[terrier]:endies[terrier],'total_med']
            aver = int(mean(b.loc[beggies[terrier]:endies[terrier], comparison]))
            colori = m.to_rgba(aver)
            singletons = singletons - set(range(beggies[terrier],endies[terrier]+1))
            begpos =  min(slicey)-midwin
            endpos = max(slicey)+midwin
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(75)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
        # let's deal with singletons:
        for single in singletons:
            medwinpos = b.loc[single,'total_med']
            begpos =  medwinpos-midwin
            endpos = medwinpos+midwin
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(75)))
            patchy = patches.Rectangle((begpos, ncomp),endpos-begpos, 1)
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))    
    ############################## 83 ####################################
    asel = compcol[(compcol >= 83) & (compcol < 90)]
    samies = list(asel.index)
    singletons = set(samies)
    if len(samies) > 0:
        print("processing 83-90% similar")
        backsies = pd.Series(samies[1:]+[-20])
        fronties = pd.Series([-20]+samies[:-1])
        normies = pd.Series(samies)
        topgo = backsies-samies
        bottomgo = samies-fronties
        tops = pd.Series([-20]+list(backsies[topgo == 1]))
        bottoms = pd.Series(list(fronties[bottomgo == 1]))
        cancels = tops-bottoms
        endies = list(tops[cancels != 0])[1:] # no need to +1 because end included with .loc
        beggies = list(bottoms[cancels != 0])        
        for terrier in range(len(endies)):
            slicey = b.loc[beggies[terrier]:endies[terrier],'total_med']
            aver = int(mean(b.loc[beggies[terrier]:endies[terrier], comparison]))
            colori = m.to_rgba(aver)
            singletons = singletons - set(range(beggies[terrier],endies[terrier]+1))
            begpos =  min(slicey)-midwin
            endpos = max(slicey)+midwin
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(83)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
        # let's deal with singletons:
        for single in singletons:
            medwinpos = b.loc[single,'total_med']
            begpos =  medwinpos-midwin
            endpos = medwinpos+midwin
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(83)))
            patchy = patches.Rectangle((begpos, ncomp),endpos-begpos, 1)
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))    
    ############################## 90 #############################################
    asel = compcol[(compcol >= 90) & (compcol < 95)]
    samies = list(asel.index)
    singletons = set(samies)
    if len(samies) > 0:
        print("processing 90-95% similar")
        backsies = pd.Series(samies[1:]+[-20])
        fronties = pd.Series([-20]+samies[:-1])
        normies = pd.Series(samies)
        topgo = backsies-samies
        bottomgo = samies-fronties
        tops = pd.Series([-20]+list(backsies[topgo == 1]))
        bottoms = pd.Series(list(fronties[bottomgo == 1]))
        cancels = tops-bottoms
        endies = list(tops[cancels != 0])[1:] # no need to +1 because end included with .loc
        beggies = list(bottoms[cancels != 0])        
        for terrier in range(len(endies)):
            slicey = b.loc[beggies[terrier]:endies[terrier],'total_med']
            aver = int(mean(b.loc[beggies[terrier]:endies[terrier], comparison]))
            colori = m.to_rgba(aver)
            singletons = singletons - set(range(beggies[terrier],endies[terrier]+1))
            begpos =  min(slicey)-midwin
            endpos = max(slicey)+midwin
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(90)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
        # let's deal with singletons:
        for single in singletons:
            medwinpos = b.loc[single,'total_med']
            begpos =  medwinpos-midwin
            endpos = medwinpos+midwin
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(90)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
    ############################## 95 #############################################
    asel = compcol[(compcol >= 95) & (compcol < 99)]
    samies = list(asel.index)
    singletons = set(samies)
    if len(samies) > 0:
        print("processing 95-99% similar")
        backsies = pd.Series(samies[1:]+[-20])
        fronties = pd.Series([-20]+samies[:-1])
        normies = pd.Series(samies)
        topgo = backsies-samies
        bottomgo = samies-fronties
        tops = pd.Series([-20]+list(backsies[topgo == 1]))
        bottoms = pd.Series(list(fronties[bottomgo == 1]))
        cancels = tops-bottoms
        endies = list(tops[cancels != 0])[1:] # no need to +1 because end included with .loc
        beggies = list(bottoms[cancels != 0])        
        for terrier in range(len(endies)):
            slicey = b.loc[beggies[terrier]:endies[terrier],'total_med']
            aver = int(mean(b.loc[beggies[terrier]:endies[terrier], comparison]))
            colori = m.to_rgba(aver)
            singletons = singletons - set(range(beggies[terrier],endies[terrier]+1))
            begpos =  min(slicey)-midwin
            endpos = max(slicey)+midwin
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(95)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
        # let's deal with singletons:
        for single in singletons:
            medwinpos = b.loc[single,'total_med']
            begpos =  medwinpos-midwin
            endpos = medwinpos+midwin
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(95)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
    ############################# 99 ########################################
    # these first, then the ones that are the same
    asel = compcol[compcol >= 99]
    samies = list(asel.index)
    singletons = set(samies)
    # let's plot these high values first
    # start with a boolean approach to identify long fragments of similarity
    if len(samies) > 0:
        print("processing 99% similar and higher")
        backsies = pd.Series(samies[1:]+[-20])
        fronties = pd.Series([-20]+samies[:-1])
        normies = pd.Series(samies)
        topgo = backsies-samies
        bottomgo = samies-fronties
        tops = pd.Series([-20]+list(backsies[topgo == 1]))
        bottoms = pd.Series(list(fronties[bottomgo == 1]))
        cancels = tops-bottoms
        endies = list(tops[cancels != 0])[1:] # no need to +1 because end included with .loc
        beggies = list(bottoms[cancels != 0])
        keep_score = set([])        
        for terrier in range(len(endies)):
            slicey = b.loc[beggies[terrier]:endies[terrier],'total_med']
            aver = int(mean(b.loc[beggies[terrier]:endies[terrier], comparison]))
            colori = m.to_rgba(aver)
            singletons = singletons - set(range(beggies[terrier],endies[terrier]+1))
            begpos =  min(slicey)-midwin
            endpos = max(slicey)+midwin
            ############ now let's go onsies
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(99)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
        # let's deal with singletons:
        for single in singletons:
            medwinpos = b.loc[single,'total_med']
            begpos =  medwinpos-midwin
            endpos = medwinpos+midwin
            ############ now let's go onsies
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(99)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
outfiley.close()

boxcount.to_csv(outdir+"/skipes.txt", sep="\t")

fig.savefig(outdir+"/pairwise_heat.pdf", bbox_inches="tight")

