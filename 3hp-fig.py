import json
from math import floor, ceil
from statistics import mean, median, stdev
from os import listdir

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

FIXED_KEY = None # (pathway, seed) or None
HIGH_RESOLUTION = True

fontLight = mpl.font_manager.FontProperties()
fontRegular = mpl.font_manager.FontProperties()

fontRegular.set_size(20)
fontLight.set_size(20)

def createPlotFromSources(leftTitle, title, sources):
    global fontRegular, fontLight

    coupling = {} # group: [values, ...]
    disruption = {} # group: [values, ...]
    for group, filename in sources.items():
        f = open(filename, 'r')
        data = json.load(f)
        f.close()

        coupling[group] = data['coupling']
        disruption[group] = data['disruption']

    groups = list(sources.keys()) # doesn't preserve order necessarily

    print(title)
    if True:
        for group in groups:
            samples = disruption[group]
            print('%s: %.2f%% mean, %.2f%% median, %.2f%% stdev' %
                  (group, 100 * mean(samples), 100 * median(samples), 100 * stdev(samples)))
    else:
        print(','.join('%.2f%%' % (100 * median(disruption[group])) for group in groups))
    print()

    fig, axs = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [5, 1]},
                            figsize=(12, 8), dpi=(100 if not HIGH_RESOLUTION else 300))
    scatterAx, distAx = axs

    prevFontSize = fontRegular.get_size()
    fontRegular.set_size(prevFontSize * 1.2)
    fig.suptitle(title.upper(), fontproperties=fontRegular)
    fig.text(0, 0.98, leftTitle, { # using suptitle defaults
        'fontsize': mpl.rcParams['figure.titlesize'],
        'fontweight': mpl.rcParams['figure.titleweight'],
        'verticalalignment': 'top',
        'horizontalalignment': 'left'
    }, fontproperties=fontRegular)
    fontRegular.set_size(prevFontSize)

    minY = min(min(disruption[group]) for group in groups)
    minY = (floor if minY < 0 else ceil)(minY / 0.1) * 0.1
    scatterAx.set_ylim(minY - 0.025, 0 + 0.025)

    scatterAx.set_xticklabels(scatterAx.get_xticks(), font_properties=fontLight)
    scatterAx.set_yticklabels(scatterAx.get_yticks(), font_properties=fontLight)
    distAx.set_xticklabels(distAx.get_xticks(), font_properties=fontLight)
    distAx.set_yticklabels(distAx.get_yticks(), font_properties=fontLight)

    labelpad = 15
    scatterAx.set_xlabel('Mean Coupling', fontproperties=fontRegular, labelpad=labelpad)
    scatterAx.set_ylabel('Change in 3-HP Yield', fontproperties=fontRegular, labelpad=labelpad)
    scatterAx.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))
    scatterAx.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
    scatterAx.grid(linestyle=':')

    for group in groups:
        scatterAx.scatter(x=coupling[group], y=disruption[group], s=4, alpha=0.5, edgecolors="none")

    distAx.set_xlabel('Density', fontproperties=fontRegular, labelpad=labelpad)
    distAx.yaxis.set_tick_params(labelbottom=True)
    #distAx.yaxis.set_label_position('right')
    distAx.yaxis.tick_right()
    distAx.yaxis.set_ticks_position('none')
    distAx.yaxis.grid(linestyle=':')

    for group in groups:
        sns.distplot(disruption[group], vertical=True, ax=distAx)

    scatterAx.legend(groups, markerscale=8, prop=fontRegular)

    plt.tight_layout(pad=3, w_pad=0.5, rect=(0, 0, 1, 0.95))

sns.set_style('white')

path = 'output'
keys = set()
if FIXED_KEY is None:
    for name in listdir(path):
        if not name.startswith('3hp-') or not name.endswith('.json'):
            print('Skipped', name)
            continue
        p = None
        seed = None
        for segment in name.split('-'):
            if segment.startswith('seed'):
                assert(seed is None)
                seed = int(segment[4:])
            if segment.startswith('p'):
                assert(p is None)
                p = int(segment[1:])
                assert(p in [1, 2])
        assert((p is None) == (seed is None))
        if p is not None:
            keys.add((p, seed))
else:
    keys = set([FIXED_KEY])

for key in sorted(keys):
    title = 'Pathway %d, Seed %d' % key
    print()
    print(title)
    print('=' * len(title))

    createPlotFromSources('A)' if key[0] == 1 else 'D)', 'Pathway %d: Scenario 1' % key[0], {
        '10% reactions active': path + '/3hp-p%d-seed%d-s1-10%%.json' % key,
        '25% reactions active': path + '/3hp-p%d-seed%d-s1-25%%.json' % key,
        '50% reactions active': path + '/3hp-p%d-seed%d-s1-50%%.json' % key
    })
    plt.savefig('output/fig-3hp-p%d-s1-seed%d.png' % key) # might use SVG for print copy

    createPlotFromSources('B)' if key[0] == 1 else 'E)', 'Pathway %d: Scenario 2' % key[0], {
        '10% reactions active': path + '/3hp-p%d-seed%d-s2-10%%.json' % key,
        '25% reactions active': path + '/3hp-p%d-seed%d-s2-25%%.json' % key,
        '50% reactions active': path + '/3hp-p%d-seed%d-s2-50%%.json' % key
    })
    plt.savefig('output/fig-3hp-p%d-s2-seed%d.png' % key)

    createPlotFromSources('C)' if key[0] == 1 else 'F)', 'Pathway %d: S1 and S2 Compared' % key[0], {
        'Only S1, 10 reactions active': path + '/3hp-p%d-seed%d-s1-10.json' % key,
        'Only S2, 10 reactions active': path + '/3hp-p%d-seed%d-s2-10.json' % key,
        'Both, 5 reactions active for each': path + '/3hp-p%d-seed%d-both-5.json' % key
    })
    plt.savefig('output/fig-3hp-p%d-both-seed%d.png' % key)
