import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import patches
from decimal import Decimal

rcParams['font.sans-serif'] = 'Arial'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

BEG = int(sys.argv[1])
END = int(sys.argv[2])
SCHEMA = sys.argv[3]
TICKGAP = int(sys.argv[4])
OUTDIR = sys.argv[5]

def arange(start, end, step):
    """A range function that works for floats"""
    split_start = str(start).split(".")
    if len(split_start) > 1:
        dec_start = len(split_start[1])
    else:
        dec_start = 0
    split_step = str(step).split(".")
    if len(split_step) > 1:
        dec_step = len(split_step[1])
    else:
        dec_step = 0
    if dec_step > dec_start:
        decider = str(step)
    else:
        decider = str(start)
    final_list = []
    current = float(start)
    while current < end:
        current = int(current) if float(current) == int(current) else current
        final_list.append(current)
        current = current + step
        current = float(Decimal(current).quantize(Decimal(decider)))
    return(final_list)


def afill(start, end, ntries):
    """A function that fill evenly spaced values between two numbers"""
    step = (end-start)/float(ntries+1) if ntries > 0 else 0
    final_list = [float(start) + (i+1)*step for i in range(ntries)]
    return(final_list)

def heatmap_panel(beg, end, schema, tickgap, outdir):
    baby = [beg] + afill(beg, end, 100) + [end]
    heighta = baby[1] - baby[0]
    widtha = len(baby)/10.0
    fig = plt.figure(figsize=(1, 5), dpi=10000)
    axes = plt.axes()
    norm = matplotlib.colors.Normalize(vmin=60, vmax=100)
    mappa = cm.ScalarMappable(norm=norm, cmap=schema)
    for scream in baby:
        new_patch = patches.Rectangle(xy = (0, scream), width = widtha, height = heighta, fc=mappa.to_rgba(scream), ec="none")
        axes.add_patch(new_patch)
    axes.set_xticks([])
    ytickies = arange(beg, end, tickgap) + [end]
    ylabbies = [str(ytick) for ytick in ytickies]
    ytickies = [ytick + heighta/2.0 if ytick not in (beg, end) else ytick for ytick in ytickies ]
    ytickies = [ytick + heighta if ytick == end else ytick for ytick in ytickies]
    print("where ticks are at: ", ytickies)
    axes.set_ylim(75,100+heighta)
    axes.set_yticks(ytickies)
    axes.set_yticklabels(ylabbies, size=15)
    outpdf = outdir+"/legend.pdf"
    fig.savefig(outpdf, bbox_inches="tight")


heatmap_panel(BEG, END, SCHEMA, TICKGAP, OUTDIR)

