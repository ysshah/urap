#!/usr/bin/env python
"""
Analysis of mystery on 8.2.0_44b.

@author Yash Shah
"""
import numpy as np, pickle, os, sys, json
from matplotlib.pyplot import *
import matplotlib.patches as mpatches
from matplotlib.widgets import CheckButtons
from matplotlib.collections import PatchCollection
from matplotlib.collections import RegularPolyCollection
import AnalysisBackend.whwp.libinspect as libinspect
from AnalysisBackend.whwp.lib_low_freq import psd_and_bin
from timestream import timestream

date = '20140921_010507'
hdf5 = '/scratch/ngoecknerwald/resonance_timestreams/%s.hdf5' %date
gain = '/scratch/ngoecknerwald/largepatch/gain/'
# cuts = '/data/pb1/neil/alpha_cut_all/%s.pkl' %date
cuts = '/scratch/ngoecknerwald/resonance_timestreams/alpha_cut_yash/%s.pkl' %date
hwMap = pickle.load(open('PB1_hardwaremap.pkl'))

tChan = np.array([hwMap['boloid'].index('8.2.0_%dt' %i) for i in range(1,92)])
bChan = np.array([hwMap['boloid'].index('8.2.0_%db' %i) for i in range(1,92)])
allChannels = np.concatenate((tChan, bChan))

cut = pickle.load(open(cuts))
flags = cut[0]

BPF_LO = 8.0
BPF_HI = 15.0

def plotWafer(ax):
    """Plot top and bottom pixels of wafer 8.2.0 on axis AX.

    Returns
    -------
    tpatches, bpatches : numpy.array
        Numpy arrays of the top and bottom patches.
    """
    def addSquare(ax, x0, y0, i, w):
        """Add square I to AX at position (X0, Y0) with wafer number W. Returns
        the matplotlib.patches patch plotted.
        """
        x = x0 - i*s/2 - i*xsep/2
        y = y0 - i*s
        ax.text(x + s/2, y + s/2, str(wafer[w]), ha='center', va='center')
        return mpatches.Rectangle((x, y), s, s, picker=True)

    def wafer8_2_0(ax, x0, y0):
        """Plot wafer 8.2.0 on AX at position (X0, Y0). Returns a numpy.array
        of all matplotlib.patches plotted.
        """
        patches = []
        w = 0

        for length in range(6,12):
            x0 += s + xsep
            for i in range(length):
                patches.append(addSquare(ax, x0, y0, i, w))
                y0 -= ysep
                w += 1
            y0 += length * ysep
        y0 -= ysep

        for length in range(10, 5, -1):
            x0 += (s + xsep) / 2
            y0 -= s
            for i in range(length):
                patches.append(addSquare(ax, x0, y0, i, w))
                w += 1
                y0 -= ysep
            y0 += (length-1) * ysep

        return np.array(patches)

    wafer = np.arange(1, 92)

    s = 1.0
    xsep = 0.1
    ysep = 0.1
    hexSep = 14*s

    x0 = -s - xsep
    y0 = 0

    ax.text(3*s + 2.5*xsep, 1.5*s, "8.2.0 top", ha='center')
    ax.text(hexSep + 3*s + 2.5*xsep, 1.5*s, "8.2.0 bottom", ha='center')
    tpatches = wafer8_2_0(ax, x0, y0)
    bpatches = wafer8_2_0(ax, x0 + hexSep, y0)

    return tpatches, bpatches


def eigenGrid(labels, amplitudes, component, full=False):
    """Plot and color an eigenvector amplitude grid.

    Returns
    -------
    fig, ax : matplotlib Figure and Axis of the eigengrid

    * tcolors, bcolors : numpy.array
        Numpy arrays of the eigenvector amplitudes of the top and bottom
        pixels, where 0.0 indicates a dark pixel (Only returned if full=True)
    """
    pixels = map(lambda s: s[6:], labels)
    tcolors = np.zeros(91)
    bcolors = np.zeros(91)

    for i, p in enumerate(pixels):
        num = int(p[:-1]) - 1
        if p[-1] == 't':
            tcolors[num] = abs(amplitudes[component][i])
        elif p[-1] == 'b':
            bcolors[num] = abs(amplitudes[component][i])

    fig, ax = subplots(figsize=(14, 7))
    tpatches, bpatches = plotWafer(ax)

    lightPatches = np.concatenate((tpatches[tcolors != 0.0],
                                   bpatches[bcolors != 0.0]))
    darkPatches = np.concatenate((tpatches[tcolors == 0.0],
                                  bpatches[bcolors == 0.0]))
    lightColors = np.concatenate((tcolors[tcolors != 0.0],
                                  bcolors[bcolors != 0.0]))

    lightCollection = PatchCollection(lightPatches, cmap=cm.jet, alpha=0.5,
        picker=True)
    lightCollection.set_array(lightColors)
    ax.add_collection(lightCollection)
    darkCollection = PatchCollection(darkPatches, facecolors='white')
    ax.add_collection(darkCollection)
    fig.colorbar(lightCollection)

    axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Eigenvector Amplitudes of Component %d' %component)

    if full:
        return fig, ax, tcolors, bcolors
    else:
        return fig, ax


def plotEigenGrid(component, amplitudes, labels, save=False, loc=''):
    """Plot the eigengrid of COMPONENT. Save the figure if SAVE and LOC."""
    fig, ax = eigenGrid(labels, amplitudes, component)
    saveOrShow(fig, save, loc, 'grid%d.png' % component)


def plotEigenGridAndComponent(component, data, components, amplitudes, labels):
    """Plot the eigengrid of COMPONENT, alongside the PSD of the component.
    Click individual pixels to overlay a plot of their PSD on the plot of the
    component's PSD; click the pixels again to remove their plot.
    """
    fig, ax, tcolors, bcolors = eigenGrid(labels, amplitudes, component, True)
    lightChannels = np.concatenate((tChan[tcolors != 0.0],
                                    bChan[bcolors != 0.0]))

    # Component plot
    fig2, ax2 = subplots()
    freq, psd, ones = psd_and_bin(components[component],
                                  data['sample_rate'], BPF_LO, BPF_HI)
    ax2.plot(freq, psd, label='Component %d' %component)
    ax2.set_title('Component %d' %component)
    ax2.set_xlabel('Frequency')
    ax2.set_ylabel('Power')
    ax2.set_yscale('log')
    ax2.margins(x=0.02, y=0.02)
    lines = {}

    def onpick(event):
        channel = lightChannels[event.ind[0]]
        boloid = hwMap['boloid'][channel]
        if channel in lines:
            lines[channel].remove()
            lines.pop(channel)
            print 'Removing', boloid
        else:
            freq, psd, ones = psd_and_bin(data['channels'][channel]['nohwps'],
                data['sample_rate'], BPF_LO, BPF_HI)
            lines[channel], = ax2.plot(freq, psd, label='Bolo ID %s' %boloid)
            print 'Plotting', boloid
        ax2.legend(loc='lower right')
        ax2.margins(x=0.02, y=0.02)
        fig2.canvas.draw()
    fig.canvas.mpl_connect('pick_event', onpick)

    show()


def flagGrid(numFlags):
    fig, ax = subplots(figsize=(14, 7))
    tpatches, bpatches = plotWafer(ax)
    patches = np.concatenate((tpatches, bpatches))

    cmax = np.max(numFlags)
    cmin = np.min(numFlags)

    cmap = get_cmap('jet', cmax-cmin+1)

    collection = PatchCollection(patches, cmap=cmap, alpha=0.5, picker=True)
    collection.set_clim(vmin=cmin-0.5, vmax=cmax+0.5)
    collection.set_array(numFlags)
    ax.add_collection(collection)
    colorbar(collection, ticks=np.arange(cmin,cmax+1))

    def onpick(event):
        channel = allChannels[event.ind[0]]
        print '\nFlags for Bolo ID ' + hwMap['boloid'][channel] + ':'
        for i in range(18, 39):
            flag = (flags[channel] & (1 << i)) >> i
            if flag[-1]:
                print '\x1b[31;1mBit %2d - %d - %s\x1b[0m' %(i, flag[-1], cut[1][i])
            else:
                print 'Bit %2d - %d - %s' %(i, flag[-1], cut[1][i])
    fig.canvas.mpl_connect('pick_event', onpick)

    axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])

    return fig, ax


def getOffsets():
    s = 1.0
    xsep = 0.1
    ysep = 0.1

    def patch(x0, y0):
        offsets = []
        for length in range(6,12):
            x0 += s + xsep
            for i in range(length):
                x = x0 - i*s/2 - i*xsep/2
                y = y0 - i*s
                offsets.append([x, y])
                y0 -= ysep
            y0 += length * ysep
        y0 -= ysep
        for length in range(10, 5, -1):
            x0 += (s + xsep) / 2
            y0 -= s
            for i in range(length):
                x = x0 - i*s/2 - i*xsep/2
                y = y0 - i*s
                offsets.append([x, y])
                y0 -= ysep
            y0 += (length-1) * ysep
        return offsets

    return patch(0, 0) + patch(14, 0)


def plotAllFlags():
    fig, ax = subplots(figsize=(16,8))

    offsets = getOffsets()

    start = 22
    numFlags = [(flags[allChannels][:,-1] & (1 << start)) >> start]

    cmax = np.max(numFlags[0])
    cmin = np.min(numFlags[0])
    cmap = get_cmap('jet', cmax-cmin+1)

    collection = RegularPolyCollection(
        numsides=4,
        rotation=np.pi/4,
        sizes=(1000,),
        linewidths = (1,),
        offsets = offsets,
        transOffset = ax.transData,
        alpha=0.5,
        cmap=cmap
        )

    collection.set_clim(vmin=cmin-0.5, vmax=cmax+0.5)

    for i, w in enumerate(range(91) + range(91)):
        ax.annotate(str(w + 1), offsets[i], ha='center', va='center')

    subplots_adjust(left=0.2)

    ax.add_collection(collection)
    axis('equal')

    collection.set_array(numFlags[0])
    cb = colorbar(collection, ticks=np.arange(cmin,cmax+1))

    rax = axes([0.05, 0.1, 0.1, 0.8])

    status = [True, False, False]
    check = CheckButtons(rax, ('22', '23', '24'), tuple(status))

    def func(label):
        bit = int(label)
        i = bit - start
        status[i] = not status[i]
        if status[i]:
            numFlags[0] += (flags[allChannels][:,-1] & (1 << bit)) >> bit
        else:
            numFlags[0] -= (flags[allChannels][:,-1] & (1 << bit)) >> bit

        cmax = np.max(numFlags[0])
        cmin = np.min(numFlags[0])
        cmap = get_cmap('jet', cmax-cmin+1)

        ax.collections[0].set_array(numFlags[0])
        ax.collections[0].set_clim(vmin=cmin-0.5, vmax=cmax+0.5)
        ax.collections[0].set_cmap(cmap)
        cb.set_ticks(np.arange(cmin, cmax+1))
        fig.canvas.draw()

    check.on_clicked(func)
    show()

    # fig.subplots_adjust(left=0.2)
    # rax = axes([0.05, 0.1, 0.1, 0.8])
    # status = [True] + [False]*20
    # check = CheckButtons(rax, [str(i) for i in range(18, 39)], status)


def plotCutGrid(save=False, loc=''):
    numFlags = np.zeros(91 + 91)
    for i in range(18, 39):
        numFlags += (flags[allChannels][:,-1] & (1 << i)) >> i
    fig, ax = flagGrid(numFlags)
    ax.set_title('All flags for %s' %date)
    saveOrShow(fig, save, loc, 'cutgrid.png')


def plotFlagGrid(flagBit, save=False, loc=''):
    numFlags = (flags[allChannels][:,-1] & (1 << flagBit)) >> flagBit
    if np.any(numFlags):
        fig, ax = flagGrid(numFlags)
        ax.set_title('Flag %d - %s' %(flagBit, cut[1][flagBit]), fontsize=12)
        saveOrShow(fig, save, loc, 'bit%d.png' %flagBit)
    else:
        print 'No flags set for bit %d - %s' %(flagBit, cut[1][flagBit])


def mainPSD():
    """Do a PCA using the PSDs of the data."""
    if os.path.isfile('data.pkl'):
        data = pickle.load(open('data.pkl'))
    else:
        d = libinspect.make_args('--input %s --gain %s --detrend 2 --nofft'
            +' --wafer %s' %(hdf5, gain, '8.2.0'))
        data = libinspect.load_TOD(d, do_PSD=True)

    ch = data['channels'].keys()
    labels = [hwMap['boloid'][c] for c in ch]

    # Channel number for Bolo ID 8.2.0_44b
    channel = 652
    # Get frequency (x axis)
    freq = data['channels'][channel]['freq_nohwps']
    # Make a 2D array of Channels x Data
    m = np.array([data['channels'][c]['psd_nohwps'] for c in ch])
    # Covariance matrix
    cov = np.cov(m)
    # Compute eigenvalues and eigenvectors
    evals, evecs = np.linalg.eig(cov)
    # 2D array of Data x Channels
    mT = m.transpose()
    # Principle components of data
    components = [np.dot(mT, evecs[:,i]) for i in range(len(ch))]
    # Eigenvector amplitudes
    # amplitudes = [[np.dot(comp, data['channels'][c]['psd_nohwps']
        # ) for c in ch] for comp in components]

    n = np.linalg.norm(m, axis=1)
    mNorm = np.array([m[i] / n[i] for i in range(len(m))])
    amplitudes = [[np.dot(comp, psd) for psd in mNorm] for comp in components]

    # hold('on')

    anomalyComponents = [0, 1, 3, 5, 7, 15, 16, 17, 19, 20, 21, 25, 26, 28, 31,
        32, 35, 37, 39, 41, 43, 44, 48, 49, 50, 51, 54, 55, 59, 61, 65, 68, 69,
        70, 71, 74, 78]


def plotAllComponents(compPSD, initialComp = 0):
    """Plot the PSD of each principle component.
    Use up and down arrow keys to switch the plotted component.
    """
    fig, ax = subplots()
    global num
    num = initialComp

    def press(event):
        if event.key == 'up' or event.key == 'down':
            global num
            num = num - 1 if event.key == 'up' else num + 1
            if num == -1:
                num = len(compPSD) - 1
            elif num == len(compPSD):
                num = 0
            ax.get_lines()[0].remove()
            ax.plot(compPSD[num][0], compPSD[num][1])
            ax.set_title('PSD of component %d' %num)
            fig.canvas.draw()

    fig.canvas.mpl_connect('key_press_event', press)
    line, = ax.plot(compPSD[num][0], compPSD[num][1])
    ax.set_yscale('log')
    ax.set_title('PSD of component %d' %num)
    xlabel('Frequency')
    ylabel('Power')
    show()


def saveOrShow(fig, save, loc, filename):
    """Depending on SAVE boolean, save a FIG to LOC+FILENAME, or show it."""
    if save:
        if loc:
            fig.savefig(loc + filename, dpi=200)
        else:
            print 'Please supply save location'
        close(fig)
    else:
        fig.show()

def saveAllCompPSD(compPSD, components, loc):
    """Save the PSDs of all principle COMPONENTS."""
    for c in components:
        figure()
        plot(compPSD[c][0], compPSD[c][1])
        yscale('log')
        title('PSD of component %d' %c)
        xlabel('Frequency')
        ylabel('Power')
        savefig('%scomp%d.png' %(loc, c), dpi=200)
        close()
        print 'Saved component', c

def saveAllEigenGrids(components, amplitudes, labels, loc):
    """Save all eigengrids of COMPONENTS."""
    for c in components:
        plotEigenGrid(c, amplitudes, labels, save=True, loc=loc)
        print 'Saved eigengrid of component', c

def rms(x):
    """Return the root-mean-square of X."""
    return np.sqrt(np.mean(x**2))

def out(s):
    """Write string S to stdout."""
    sys.stdout.write(s)
    sys.stdout.flush()


def dataToJSON(compPSD):
    with open('data.csv', 'wb') as f:
        f.write(','.join(['freq'] + ['comp%d' %c for c in range(len(compPSD))]) + '\n')
        for i,freq in enumerate(compPSD[0][0]):
            f.write(','.join([str(freq)] + [str(compPSD[c][1][i]) for c in range(len(compPSD))]) + '\n')
    # d = {}
    # d["freq"] = list(compPSD[0][0])
    # for i in range(len(compPSD)):
    #     d[i] = list(compPSD[i][1])

def findOutliers(amplitudes, labels):
    weight = 30
    for i in range(len(amplitudes)):
        outliers = amplitudes[i] > weight*np.median(amplitudes[i])
        numOutliers = np.sum(outliers)
        if numOutliers:
            print 'Component %2d - %2d - %s' %(i, numOutliers, labels[outliers])

if __name__ == '__main__':

    if os.path.isfile('data_ts.pkl'):
        out('\nLoading data from data_ts.pkl...')
        data = pickle.load(open('data_ts.pkl'))
        out(' Done.\n')
    else:
        # d = libinspect.make_args('--input %s --gain %s --detrend 2 --nofft'
        #     +' --wafer %s' %(hdf5, gain, '8.2.0'))
        # data = libinspect.load_TOD(d, no_demod=True)
        ch8_2_0 = [
            513, 514, 515, 516, 517, 518, 519, 520, 556, 522, 523, 652, 526,
            527, 528, 529, 531, 532, 533, 534, 537, 538, 667, 668, 669, 542,
            543, 672, 673, 549, 666, 509, 555, 684, 685, 558, 559, 688, 689,
            690, 563, 564, 565, 567, 568, 569, 570, 571, 572, 573, 574, 706,
            709, 585, 586, 678, 594, 595, 596, 603, 604, 605, 606, 677, 609,
            550, 691, 619, 613, 562, 630, 631, 504, 506, 507, 508, 618, 510,
            511
        ]
        d = libinspect.make_args('--input %s --gain %s --detrend 2 --nofft' %(
            hdf5, gain))
        data = libinspect.load_TOD(d, channels=ch8_2_0)

    channels = np.array(data['channels'].keys())
    time = data['time']

    labels = np.array([hwMap['boloid'][channel] for channel in channels])

    bpf = timestream.make_bpf(BPF_LO, BPF_HI, n=1024, nyq=data['sample_rate']/2.)
    dataMatrix = np.array([timestream.apply_filter(
        bpf, data['channels'][channel]['nohwps']) for channel in channels])

    out('Computing principle components...')
    covarianceMatrix = np.cov(dataMatrix)
    evals, evecs = np.linalg.eig(covarianceMatrix)
    components = [np.dot(dataMatrix.T, evecs[:,i]
        ) for i in range(len(channels))]
    components.sort(key=lambda x: rms(x), reverse=True)
    out(' Done.\n')

    out('Calculating eigenvector amplitudes...')
    amplitudes = [[rms((np.dot(comp, ts) / np.dot(comp, comp)) * comp
        ) for ts in dataMatrix] for comp in components]
    out(' Done.\n')

    out('Creating PSDs of components...')
    compPSD = [psd_and_bin(ts, data['sample_rate'], BPF_LO, BPF_HI
        )[:2] for ts in components]
    out(' Done.\n')

    out('Creating PSDs of data...')
    dataPSD = [psd_and_bin(data['channels'][c]['nohwps'], data['sample_rate'],
        BPF_LO, BPF_HI)[:2] for c in channels]
    out(' Done.\n')

    bandFeatureComponents = [0]
    spike12HzComponents = [1, 2, 3, 4, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19,
        21, 22, 27, 28]

    if False:
        # Count number that principle components reduced to depending on weight
        # of median comparison
        for w in range(10, 101, 10):
            count = 0
            for i in range(len(amplitudes)):
                if np.any(amplitudes[i] > w*np.median(amplitudes[i])):
                    count += 1
            print '%3d   %2d' %(w, count)
