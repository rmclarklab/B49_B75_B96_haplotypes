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

infile = sys.argv[1] # which file to go over
window = int(sys.argv[2]) # what size window was used?
slide = int(sys.argv[3]) # what slide
bottompanel = int(sys.argv[4]) # this is how often to plot points on the bottom panel
scaffold_file = sys.argv[5] # file with the scaffold and length
color_scheme = sys.argv[6] # what color scheme to use?
outdir = sys.argv[7] # where to output the pdf

# this is to process the file and output the shading correctly
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

# processing the scaffold information
shader = scale_dict(scaffold_file)
maxx = sum(pd.read_csv(scaffold_file, sep="\t", header=None)[1])
print(maxx)

midwin = window/2

outfiley = open(outdir+"/connected_fragments.txt", "w")

b = pd.read_csv(infile)
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

grids = sum(range(nstrains+4))+25
print("grids", grids)

fig = plot.figure(figsize=(grids, grids), dpi=10000)

gs = gridspec.GridSpec(grids, grids)

begger = 1
mover = nstrains-1
ender = begger + mover

lastfirst = None

boxcount = pd.Series()

## let's do colors
## have to normalize the data set first so that they're continuous

norm = matplotlib.colors.Normalize(vmin=75, vmax=100)
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
        ax.title.set_fontsize(80) # this is very awkward, but it the only way to do this
        ax.set_ylim(0,ender-begger)
        ax.set_xlim(0, max(b['total_med']))
        yloc = plticker.FixedLocator(tickies)
        ax.yaxis.set_ticklabels(remains)
        ax.yaxis.set_tick_params(labelsize=60)
        ax.yaxis.set_major_locator(yloc)
        xloc = plticker.MultipleLocator(2000000)
        ax.xaxis.set_major_locator(xloc)
        lockyx = list(ax.xaxis.get_major_locator().tick_values(0,max(b['total_med'])))
        ax.xaxis.set_tick_params(labelsize=150, colors='black')
        # for better labeling:
        new_lockyx = [int(i) for i in lockyx] # this is to create labels with numbers
        xlabs = []
        for i in new_lockyx:
            j = str(i)
            if len(str(j)) <= 2:
                xlabs.append(i/1)
            elif 3 <= len(str(j)) <= 6:
                xlabs.append(i/1000)
            elif  3 <= len(str(j)) <= 9:
                xlabs.append(i/1000000)
            else:
                xlabs.append(round(i/float(1000000000), 1))
        ax.xaxis.set_ticklabels(xlabs)
        # this are the variables for the shading below
        old = None
        shade = True
        # here comes the shading
        for contig in sorted(shader):
            val = shader[contig]
            if old != None and shade == True:
                plot.axvspan(old, val, color='0.85', alpha=0.5)
                shade = False
            else:
                if old != None:
                    shade = True
            old = shader[contig]
        # the last one    
        if shade == True:
            plot.axvspan(old, maxx, color='0.85', alpha=0.5)
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
            # label offset needs to be .2 to worky well
            plotsizies = ender - begger
            tickplace = (plotsizies+0.25)/plotsizies
            ax.set_title(first, y=tickplace)
            ####
            ax.title.set_fontsize(80) # this is very awkward, but it the only way to do this
            ax.set_ylim(0,ender-begger)
            ax.set_xlim(0, max(b['total_med']))
            yloc = plticker.FixedLocator(tickies)
            ax.yaxis.set_ticklabels(remains)
            ax.yaxis.set_tick_params(labelsize=60)
            ax.yaxis.set_major_locator(yloc)
            xloc = plticker.MultipleLocator(2000000)
            ax.xaxis.set_major_locator(xloc)
            lockyx = list(ax.xaxis.get_major_locator().tick_values(0,max(b['total_med'])))
            ax.xaxis.set_tick_params(labelsize=150, colors='black')
            # for better labeling:
            new_lockyx = [int(i) for i in lockyx] # this is to create labels with numbers
            xlabs = []
            for i in new_lockyx:
                j = str(i)
                if len(str(j)) <= 2:
                    xlabs.append(i/1)
                elif 3 <= len(str(j)) <= 6:
                    xlabs.append(i/1000)
                elif  3 <= len(str(j)) <= 9:
                    xlabs.append(i/1000000)
                else:
                    xlabs.append(round(i/float(1000000000), 1))
            ax.xaxis.set_ticklabels(xlabs)
            # this are the variables for the shading below
            old = None
            shade = True
            # here comes the shading
            for contig in sorted(shader):
                val = shader[contig]
                if old != None and shade == True:
                    plot.axvspan(old, val, color='0.85', alpha=0.5)
                    shade = False
                else:
                    if old != None:
                        shade = True
                old = shader[contig]
            # the last one    
            if shade == True:
                plot.axvspan(old, maxx, color='0.85', alpha=0.5)
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
    # these first, then the ones that are the same
    asel = compcol[compcol >= 95]
    samies = list(asel.index)
    singletons = set(samies)
    # let's plot these high values first
    # start with a boolean approach to identify long fragments of similarity
    if len(samies) > 0:
        print("processing 95% similar and higher")
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

            ########## this is for frequency plotting later
            testbeg = begpos%bottompanel
            if testbeg == 0:
                shift = 0
            else:
                shift = bottompanel-testbeg
            actualbeg = begpos + shift
            #end
            testend = endpos%bottompanel
            shift = bottompanel-testend
            actualend = endpos + shift
            # now get the stuff
            stellar = [i for i in range(actualbeg, actualend, bottompanel) if i not in keep_score]
            newvals = pd.Series(stellar).value_counts()
            boxcount = boxcount.add(newvals, fill_value = 0)
            keep_score = set(stellar) | keep_score

            ############ now let's go onsies
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(95)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
        # let's deal with singletons:
        for single in singletons:
            medwinpos = b.loc[single,'total_med']
            begpos =  medwinpos-midwin
            endpos = medwinpos+midwin

            ########## this is for frequency plotting later
            testbeg = begpos%bottompanel
            if testbeg == 0:
                shift = 0
            else:
                shift = bottompanel-testbeg
            actualbeg = begpos + shift
            #end
            testend = endpos%bottompanel
            shift = bottompanel-testend
            actualend = endpos + shift
            # now get the stuff
            stellar = [i for i in range(actualbeg, actualend, bottompanel) if i not in keep_score]
            newvals = pd.Series(stellar).value_counts()
            boxcount = boxcount.add(newvals, fill_value = 0)
            keep_score = set(stellar) | keep_score
            ############ now let's go onsies
            outfiley.write("%s\t%s\t%s\t%s\t%s\n" %(first, second, str(begpos), str(endpos), str(99)))
            ax.add_patch(patches.Rectangle((begpos, ncomp),endpos-begpos, 1, fc = colori, ec=None))
outfiley.close()

### final box with the plot

futureindex = range(0,int(max(boxcount.index)),bottompanel)
repa = len(futureindex)
repa_stuff = [0]*repa

greatness = pd.Series(data= repa_stuff, index = futureindex)

boxcount = boxcount.add(greatness, fill_value = 0)

boxcount = boxcount/float(ncomps)

begger = ender + 3
ender = grids - 5
ax = plot.subplot(gs[begger:ender, :])
print("last plotting")
colori = m.to_rgba(100)
ax.plot(boxcount.index, boxcount, color=colori, lw=8)
ax.set_xlim(0, max(b['total_med']))
yloc = plticker.MultipleLocator(0.1)
ax.yaxis.set_tick_params(labelsize=150)
ax.yaxis.set_major_locator(yloc)
xloc = plticker.MultipleLocator(2000000)
ax.xaxis.set_major_locator(xloc)
lockyx = list(ax.xaxis.get_major_locator().tick_values(0,max(b['total_med'])))
ax.xaxis.set_tick_params(labelsize=150, colors='black')
ax.set_xlabel("Genomic Position (Mb)", size = 120)
ax.set_ylabel("Fraction 95%+ identical", size = 120)
# for better labeling:
new_lockyx = [int(i) for i in lockyx] # this is to create labels with numbers
xlabs = []
for i in new_lockyx:
    j = str(i)
    if len(str(j)) <= 2:
        xlabs.append(i/1)
    elif 3 <= len(str(j)) <= 6:
        xlabs.append(i/1000)
    elif  3 <= len(str(j)) <= 9:
        xlabs.append(i/1000000)
    else:
        xlabs.append(round(i/float(1000000000), 1))
ax.xaxis.set_ticklabels(xlabs)
# this are the variables for the shading below
old = None
shade = True
# here comes the shading
for contig in sorted(shader):
    val = shader[contig]
    if old != None and shade == True:
        plot.axvspan(old, val, color='0.85', alpha=0.5)
        shade = False
    else:
        if old != None:
            shade = True
    old = shader[contig]
# the last one    
if shade == True:
    plot.axvspan(old, maxx, color='0.85', alpha=0.5)


boxcount.to_csv(outdir+"/skipes.txt", sep="\t")

fig.savefig(outdir+"/multiplot_heat.pdf", bbox_inches="tight")


