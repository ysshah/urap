#!/usr/bin/env python
"""
Analysis of mystery on 8.2.0_44b.

@author Yash Shah
"""
import numpy as np, pickle, os
from matplotlib.pyplot import *
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import AnalysisBackend.whwp.libinspect as libinspect
from AnalysisBackend.whwp.lib_low_freq import psd_and_bin
from timestream import timestream

date = '20140921_010507'
hdf5 = '/scratch/ngoecknerwald/resonance_timestreams/%s.hdf5' %date
gain = '/scratch/ngoecknerwald/largepatch/gain/'
cuts = '/data/pb1/neil/alpha_cut_all/%s.pkl' %date
hwMap = pickle.load(open('PB1_hardwaremap.pkl'))
bIDs = np.array(hwMap['boloid'])

def singleCase():
    # Channel number for Bolo ID 8.2.0_44b
    channel = 652

    d = libinspect.make_args('--input %s --gain %s --detrend 2 --nofft' %(
        hdf5, gain))
    data = libinspect.load_TOD(d, channels=[0, channel], do_PSD=True)
    time = data['time']

    plot(data['channels'][channel]['raw_freq'],
       data['channels'][channel]['raw_psd'])
    yscale('log')
    show()

    # plot(data['channels'][channel]['freq_nohwps'],
         # data['channels'][channel]['psd_nohwps'])
    # yscale('log')
    # show()

    

    # Demodulated data
    # plot(data['channels'][channel]['freqs_demod'][0],
    #    data['channels'][channel]['psds_demod'][0])
    # plot(data['channels'][channel]['freqs_demod'][1],
    #    data['channels'][channel]['psds_demod'][1])
    # plot(data['channels'][channel]['freqs_demod'][2],
    #    data['channels'][channel]['psds_demod'][2])
    # yscale('log')
    # show()


def plotComponents(comps, freq, components):
    figure()
    for i in comps:
        plot(freq, components[i], label='Component %d' %i)
    legend(loc='upper left')
    yscale('log')
    show()


def plotEigenAmps(comp, amplitudes, labels):
    figure()
    title('Eigenvector Amplitude of Component %d' %comp)
    xlabel('Bolo ID')
    ylabel('Eigenvector Amplitude')
    plot(range(len(amplitudes[comp])), amplitudes[comp],
        label='Component %d' %(comp + 1))
    scatter(range(len(amplitudes[comp])), amplitudes[comp],
        label='Component %d' %(comp + 1))
    xticks(range(len(amplitudes[comp])), labels, rotation='vertical')
    tight_layout()
    axis('tight')
    show()


def plotWafer(ax):
    wafer = np.arange(1, 92)

    s = 1.0
    xsep = 0.1
    ysep = 0.1
    hexSep = 14*s

    def addSquare(ax, x0, y0, i, w):
        x = x0 - i*s/2 - i*xsep/2
        y = y0 - i*s
        ax.text(x + s/2, y + s/2, str(wafer[w]), ha='center', va='center')
        return mpatches.Rectangle((x, y), s, s, picker=True)

    def wafer8_2_0(ax, x0, y0):
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

    x0 = -s - xsep
    y0 = 0

    ax.text(3*s + 2.5*xsep, 1.5*s, "8.2.0 top", ha='center')
    ax.text(hexSep + 3*s + 2.5*xsep, 1.5*s, "8.2.0 bottom", ha='center')
    tpatches = wafer8_2_0(ax, x0, y0)
    bpatches = wafer8_2_0(ax, x0 + hexSep, y0)

    return tpatches, bpatches


def plotEigenGrid(comp, data, components, amplitudes, labels):
    pixels = map(lambda s: s[6:], labels)

    tcolors = np.zeros(91)
    bcolors = np.zeros(91)

    for i, p in enumerate(pixels):
        num = int(p[:-1]) - 1
        if p[-1] == 't':
            tcolors[num] = abs(amplitudes[comp][i])
        elif p[-1] == 'b':
            bcolors[num] = abs(amplitudes[comp][i])

    fig, ax = subplots(figsize=(14, 7))
    tpatches, bpatches = plotWafer(ax)

    lightPatches = np.concatenate((tpatches[tcolors != 0.0],
                                   bpatches[bcolors != 0.0]))
    darkPatches = np.concatenate((tpatches[tcolors == 0.0],
                                  bpatches[bcolors == 0.0]))
    lightColors = np.concatenate((tcolors[tcolors != 0.0],
                                  bcolors[bcolors != 0.0]))

    lightCollection = PatchCollection(lightPatches, cmap=cm.jet, alpha=0.5)
    lightCollection.set_array(lightColors)
    ax.add_collection(lightCollection)

    darkCollection = PatchCollection(darkPatches, facecolors='white')
    ax.add_collection(darkCollection)

    colorbar(lightCollection)
    axis('equal')
    setp(ax.get_xticklabels(), visible=False)
    setp(ax.get_yticklabels(), visible=False)

    title('Eigenvector Amplitudes of Component %d' %comp)
    savefig('bpfEigenGrids/grid%d.png' %comp, dpi=200)

    close(fig)

    # figure()
    # plot(freq, components[comp], label='Component %d' %comp)
    # title('Component %d' %comp)
    # legend()
    # yscale('log')
    # savefig('eigenGrids/comp%d.png' %comp, dpi=200)
    # close()

    # figure()
    # freq, psd, ones = psd_and_bin(components[comp], data['sample_rate'], 11, 13)
    # plot(freq, psd)
    # title('Component %d' %comp)
    # xlabel('Frequency')
    # ylabel('Power')
    # yscale('log')
    # margins(x=0.02, y=0.02)
    # show()

    # print 'Saved plots for component', comp
    # show()


def plotCutGrid():
    cut = pickle.load(open(cuts))
    flags = cut[0]

    tChan = np.array([hwMap['boloid'].index('8.2.0_%dt' %i) for i in range(1,92)])
    bChan = np.array([hwMap['boloid'].index('8.2.0_%db' %i) for i in range(1,92)])
    channels = np.concatenate((tChan, bChan))

    numFlags = np.zeros(91 + 91)
    for i in range(18, 39):
        numFlags += (flags[channels][:,-1] & (1 << i)) >> i

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
        channel = channels[event.ind[0]]
        print '\nFlags for Bolo ID ' + hwMap['boloid'][channel] + ':'
        for i in range(18, 39):
            flag = (flags[channel] & (1 << i)) >> i
            if flag[-1]:
                print '\x1b[31;1m%s - %s\x1b[0m' %(flag[-1], cut[1][i])
            else:
                print '%s - %s' %(flag[-1], cut[1][i])
                

    fig.canvas.mpl_connect('pick_event', onpick)

    axis('equal')
    setp(ax.get_xticklabels(), visible=False)
    setp(ax.get_yticklabels(), visible=False)
    # savefig('flags/cutgrid.png', dpi=200)
    # close(fig)
    show()


def plotFlagGrid(flagBit):
    cut = pickle.load(open(cuts))
    flags = cut[0]

    tChan = np.array([hwMap['boloid'].index('8.2.0_%dt' %i) for i in range(1,92)])
    bChan = np.array([hwMap['boloid'].index('8.2.0_%db' %i) for i in range(1,92)])
    channels = np.concatenate((tChan, bChan))

    chanFlags = (flags[channels][:,-1] & (1 << flagBit)) >> flagBit

    if np.any(chanFlags):
        fig, ax = subplots(figsize=(14, 7))
        tpatches, bpatches = plotWafer(ax)
        patches = np.concatenate((tpatches, bpatches))

        cmax = np.max(chanFlags)
        cmin = np.min(chanFlags)

        cmap = get_cmap('jet', cmax-cmin+1)

        collection = PatchCollection(patches, cmap=cmap, alpha=0.5, picker=True)
        collection.set_clim(vmin=cmin-0.5, vmax=cmax+0.5)
        collection.set_array(chanFlags)
        ax.add_collection(collection)
        colorbar(collection, ticks=np.arange(cmin,cmax+1))

        def onpick(event):
            channel = channels[event.ind[0]]
            print '\nFlags for Bolo ID ' + hwMap['boloid'][channel] + ':'
            for i in range(18, 39):
                flag = (flags[channel] & (1 << i)) >> i
                if flag[-1]:
                    print '\x1b[31;1mBit %2d - %d - %s\x1b[0m' %(i, flag[-1], cut[1][i])
                else:
                    print 'Bit %2d - %d - %s' %(i, flag[-1], cut[1][i])

        fig.canvas.mpl_connect('pick_event', onpick)

        axis('equal')
        setp(ax.get_xticklabels(), visible=False)
        setp(ax.get_yticklabels(), visible=False)
        title('Flag %d - %s' %(flagBit, cut[1][flagBit]), fontsize=12)
        savefig('flags/bit%d.png' %flagBit, dpi=200)

        close(fig)
        print 'Saved plots for bit %d' %flagBit
        # show()
    else:
        print 'No flags set for bit %d - %s' %(flagBit, cut[1][flagBit])


# def plotPSD(channels, data):
#     figure()
#     for c in channels:
#         plot(data['channels'][c]['freq_nohwps'],
#              data['channels'][c]['psd_nohwps'], label=str(c))
#     yscale('log')
#     legend()
#     show()


def main():
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


def plotPSD(ts, data):
    freq, psd, ones = psd_and_bin(ts, data['sample_rate'], 11., 13.)
    figure()
    plot(freq, psd)
    yscale('log')
    margins(x=0.02, y=0.02)
    show()


def plotAllComponents(compPSD, initialComp = 0):
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


def saveAllCompPSD(compPSD, anomalyComponents):
    for i in anomalyComponents:
        figure()
        plot(compPSD[i][0], compPSD[i][1])
        yscale('log')
        title('PSD of component %d' %i)
        xlabel('Frequency')
        ylabel('Power')
        savefig('bpfEigenGrids/comp%d.png' %i, dpi=200)
        print 'Saved component %d' %i
        close()


if __name__ == '__main__':

    if os.path.isfile('data_ts.pkl'):
        print 'Loading data from data_ts.pkl...'
        data = pickle.load(open('data_ts.pkl'))
        print 'Done.'
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
            550, 691, 619, 613, 562, 630, 631, 504, 506, 507, 508, 618, 510, 511
        ]
        d = libinspect.make_args('--input %s --gain %s --detrend 2 --nofft' %(hdf5, gain))
        data = libinspect.load_TOD(d, channels=ch8_2_0)

    ch = data['channels'].keys()
    time = data['time']

    labels = [hwMap['boloid'][c] for c in ch]

    bpf = timestream.make_bpf(11., 13., n=1024, nyq=data['sample_rate']/2.)
    # m = np.array([timestream.apply_filter(bpf, data['channels'][c]['hwpss']) for c in ch])
    m = np.array([timestream.apply_filter(bpf, data['channels'][c]['nohwps']) for c in ch])
    cov = np.cov(m)
    evals, evecs = np.linalg.eig(cov)
    mT = m.transpose()
    components = [np.dot(mT, evecs[:,i]) for i in range(len(ch))]

    n = np.linalg.norm(m, axis=1)
    mNorm = np.array([m[i] / n[i] for i in range(len(m))])
    amplitudes = [[np.dot(comp, ts) for ts in mNorm] for comp in components]

    compPSD = [psd_and_bin(ts, data['sample_rate'], 11., 13.) for ts in components]

    anomalyComponents = [1, 2, 3, 4, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19,
        22, 23, 33, 34]
