#!/usr/bin/env python
"""
Find outliers using the eigenvector amplitudes of channels.

@author Yash Shah

Questions
---------
- Input file? Is that necessary?
- By default, which channels should we load? All channels?
- Option to load specific channels?
- Option for date range?
- Option for HDF5 file (required)?
- Do I need the timestream filtering arguments? Do those have to do with
  bandpass filter?
- Bandpass filter arguments (low freq, high freq, n)?
- W option for max vs. W * median comparison?

"""
import numpy as np, pickle, os, sys, argparse
from matplotlib.pyplot import *
import matplotlib.patches as mpatches
from matplotlib.widgets import CheckButtons
from matplotlib.collections import PatchCollection
from matplotlib.collections import RegularPolyCollection
import AnalysisBackend.whwp.libinspect as libinspect
from AnalysisBackend.whwp.lib_low_freq import psd_and_bin
from timestream import timestream

def plotPSD(freq, power, plotTitle, filename):
    figure()
    plot(freq, power)
    yscale('log')
    title(plotTitle)
    xlabel('Frequency (Hz)')
    ylabel('Power')
    savefig(filename, dpi=200)
    close()

def plotComponent(component, location, compPSDs):
    freq, power = compPSDs[component]
    plotPSD(freq, power, 'PSD of Component %d' %component,
        '%scomp%d' %(location, component))

def plotChannel(channel, location, dataPSDs, labels):
    freq, power = dataPSDs[channel]
    plotPSD(freq, power, 'PSD of Channel %d' %labels[channel],
        '%schan%s' %(location, labels[channel]))

def ces_main(args):
    lib_args = '--input %s --gain %s --detrend 2 --nofft' %(args.hdf5, args.gain)
    if args.wafer:
        lib_args += ' --wafer %s' %args.wafer
        d = libinspect.make_args(lib_args)
        data = libinspect.load_TOD(d)
    else:
        d = libinspect.make_args(lib_args)
        data = libinspect.load_TOD(d, channels=[0])

    channels = np.array(data['channels'].keys())
    time = data['time']

    bpf = timestream.make_bpf(BPF_LO, BPF_HI, n=1024, nyq=data['sample_rate']/2.)
    dataMatrix = np.array([timestream.apply_filter(
        bpf, data['channels'][channel]['nohwps']) for channel in channels])

    print 'Computing principle components...'
    covarianceMatrix = np.cov(dataMatrix)
    eigenvalues, eigenvectors = np.linalg.eig(covarianceMatrix)
    components = [np.dot(dataMatrix.T, eigenvectors[:,i]
        ) for i in range(len(channels))]
    # Sort principle components by highest RMS first
    components.sort(key=lambda x: rms(x), reverse=True)

    print 'Calculating eigenvector amplitudes...'
    amplitudes = [[rms((np.dot(comp, ts) / np.dot(comp, comp)) * comp
        ) for ts in dataMatrix] for comp in components]

    print 'Creating PSDs of components...'
    compPSDs = [psd_and_bin(ts, data['sample_rate'], BPF_LO, BPF_HI
        )[:2] for ts in components]

    print 'Creating PSDs of data...'
    dataPSDs = [psd_and_bin(data['channels'][c]['nohwps'], data['sample_rate'],
        BPF_LO, BPF_HI)[:2] for c in channels]

    weight = 30
    for i in range(len(amplitudes)):
        outliers = amplitudes[i] > weight*np.median(amplitudes[i])
        numOutliers = np.sum(outliers)
        if numOutliers:
            plotComponent(i, '', compPSDs)
            for c in channels[outliers]:
                plotChannel(c, '', dataPSDs, labels)

def addargs(parser):
    # Default arguments
    parser.add_argument('-o','--output',dest='output_folder',help='Output folder')
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f',dest='in_path',help='Input file')
    parser.add_argument('-v','--verbose',dest='verbose',action='store_true',help='Verbose output')
    parser.add_argument('--nocache',dest='cache',action='store_false',help='Only run if output file does not exist')

    # Include timestream filtering options (not sure if needed)
    libmap_whwp.add_filter_timestream_args_whwp(parser)

    # My arguments
    # parser.add_argument('--hdf5',help='HDF5 file',required=True)
    parser.add_argument('--gain',help='Gain folder',required=True)
    parser.add_argument('--wafer',help='Use this entire wafer')

def grabargs(args_param=None):
    parser = argparse.ArgumentParser(description='Find outliers using the eigenvector amplitudes of channels.')
    addargs(parser)
    args = parser.parse_args(args_param)
    return args

def main(args_param=None):
    args = grabargs(args_param)
    ces_main(args)

if __name__=='__main__':
    main()
