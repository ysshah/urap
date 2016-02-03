#!/usr/bin/env python
"""
Analysis of mystery on 8.2.0_44b.

@author Yash Shah
"""
from matplotlib.pyplot import *
import numpy as np, pickle
import AnalysisBackend.whwp.libinspect as libinspect

hdf5 = '/scratch/ngoecknerwald/resonance_timestreams/20140921_010507.hdf5'
gain = '/scratch/ngoecknerwald/largepatch/gain/'
hwMap = pickle.load(open('PB1_hardwaremap.pkl'))
bIDs = np.array(hwMap['boloid'])

def singleCase():
    # Channel number for Bolo ID 8.2.0_44b
    channel = 652

    d = libinspect.make_args('--input %s --gain %s --detrend 2 --nofft' %(hdf5, gain))
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


if __name__ == '__main__':

    d = libinspect.make_args('--input %s --gain %s --detrend 2 --nofft --wafer %s' %(
        hdf5, gain, '8.2.0'))
    data = libinspect.load_TOD(d, do_PSD=True)

    ch = data['channels'].keys()

    freq = data[channel]['freq_nohwps']

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

    hold('on')

    # First 10 principle components
    figure()
    for i in range(10):
        plot(freq, components[i], label='Component %d' %(i+1))
    legend()
    yscale('log')
    show()

    # Eigenvector amplitudes of component 4
    i = 3
    figure()
    title('Eigenvector Amplitude of Component %d' %(i+1))
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
