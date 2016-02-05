#!/usr/bin/env python
"""
Analysis of mystery on 8.2.0_44b.

@author Yash Shah
"""
import numpy as np, pickle, os
import AnalysisBackend.whwp.libinspect as libinspect
from matplotlib.pyplot import *
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

hdf5 = '/scratch/ngoecknerwald/resonance_timestreams/20140921_010507.hdf5'
gain = '/scratch/ngoecknerwald/largepatch/gain/'
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
    legend()
    yscale('log')
    show()


def plotEigenAmps(comp, amplitudes, labels):
    figure()
    title('Eigenvector Amplitude of Component %d' %i)
    xlabel('Bolo ID')
    ylabel('Eigenvector Amplitude')
    plot(range(len(amplitudes[i])), amplitudes[i],
        label='Component %d' %(i + 1))
    scatter(range(len(amplitudes[i])), amplitudes[i],
        label='Component %d' %(i + 1))
    xticks(range(len(amplitudes[i])), labels, rotation='vertical')
    tight_layout()
    axis('tight')
    show()


def plotWafer(tcolors, bcolors):
    fig, ax = subplots(figsize=(20, 10))

    wafer = np.arange(1, 92)

    s = 1.0
    xsep = 0.1
    ysep = 0.1
    hexSep = 14*s

    def addSquare(x0, y0, i, w):
        x = x0 - i*s/2 - i*xsep/2
        y = y0 - i*s
        text(x + s/2, y + s/2, str(wafer[w]), ha='center', va='center')
        return mpatches.Rectangle((x, y), s, s)

    def wafer8_2_0(x0, y0):
        patches = []
        w = 0

        for length in range(6,12):
            x0 += s + xsep
            for i in range(length):
                patches.append(addSquare(x0, y0, i, w))
                y0 -= ysep
                w += 1
            y0 += length * ysep
        y0 -= ysep

        for length in range(10, 5, -1):
            x0 += (s + xsep) / 2
            y0 -= s
            for i in range(length):
                patches.append(addSquare(x0, y0, i, w))
                w += 1
                y0 -= ysep
            y0 += (length-1) * ysep

        return np.array(patches)

    x0 = -s - xsep
    y0 = 0

    text(3*s + 2.5*xsep, 1.5*s, "8.2.0 top", ha='center')
    text(hexSep + 3*s + 2.5*xsep, 1.5*s, "8.2.0 bottom", ha='center')
    tpatches = wafer8_2_0(x0, y0)
    bpatches = wafer8_2_0(x0 + hexSep, 0)

    # import pdb
    # pdb.set_trace()

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

    # show()


def plotGrid(comp, freq, components, amplitudes, labels):
    pixels = map(lambda s: s[6:], labels)

    tcolors = np.zeros(91)
    bcolors = np.zeros(91)

    for i, p in enumerate(pixels):
        num = int(p[:-1]) - 1
        if p[-1] == 't':
            tcolors[num] = abs(amplitudes[comp][i])
        elif p[-1] == 'b':
            bcolors[num] = abs(amplitudes[comp][i])

    plotWafer(tcolors, bcolors)

    # figure()
    # plot(freq, components[comp], label='Component %d' %comp)
    # legend()
    # yscale('log')
    show()

def plotPSD(channels, data):
    figure()
    for c in channels:
        plot(data['channels'][c]['freq_nohwps'],
             data['channels'][c]['psd_nohwps'], label=str(c))
    yscale('log')
    legend()
    show()

if __name__ == '__main__':

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
    amplitudes = [[np.dot(comp, data['channels'][c]['psd_nohwps']
        ) for c in ch] for comp in components]

    # hold('on')

    anomalyComponents = [0, 3, 4, 6, 11, 12, 16, 18, 19]

    # plotComponents(range(11, 21), freq, components)

    # Eigenvector amplitudes of component 4
    # plotEigenAmps(0, amplitudes, labels)

    # plotGrid(0, freq, components, amplitudes, labels)
