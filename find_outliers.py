#!/usr/bin/env python
"""
Find outliers using the eigenvector amplitudes of channels.

@author Yash Shah

Todo
----
- Add more arguments
- Save output summary file
- Save plots

Questions
---------
- Is it ok to give channels argument if wafer is supplied?
- Where is PB1_hardwaremap.pkl?

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

def plot_PSD(freq, power, plot_title, filename):
    figure()
    plot(freq, power)
    yscale('log')
    title(plot_title)
    xlabel('Frequency (Hz)')
    ylabel('Power')
    savefig(filename, dpi=200)
    close()

def plot_component(component, location, comp_PSDs):
    freq, power = comp_PSDs[component]
    plotPSD(freq, power, 'PSD of Component %d' %component,
        '%scomp%d' %(location, component))

def plot_channel(channel, location, data_PSDs, labels):
    freq, power = data_PSDs[channel]
    plotPSD(freq, power, 'PSD of Channel %d' %labels[channel],
        '%schan%s' %(location, labels[channel]))

def ces_main(args):
    greg = os.path.splitext(os.path.split(args.in_path)[1])[0]
    stats_file = os.path.join(args.psd_folder, '%s.pkl'%greg)
    print('Loading the stats file %s'%stats_file)
    assert os.path.exists(stats_file)
    lowfreq_data = pickle.load(open(stats_file))

    data = libinspect.load_TOD(args, channels=lowfreq_data['stacked']['chanlist'])

    channels = np.array(data['channels'].keys())
    time = data['time']

    bpf = timestream.make_bpf(BPF_LO, BPF_HI, n=1024, nyq=data['sample_rate']/2.)
    data_matrix = np.array([timestream.apply_filter(
        bpf, data['channels'][channel]['nohwps']) for channel in channels])

    print 'Computing principle components...'
    covariance_matrix = np.cov(data_matrix)
    eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
    components = [np.dot(data_matrix.T, eigenvectors[:,i]
        ) for i in range(len(channels))]
    # Sort principle components by highest RMS first
    components.sort(key=lambda x: rms(x), reverse=True)

    print 'Calculating eigenvector amplitudes...'
    amplitudes = [[rms((np.dot(comp, ts) / np.dot(comp, comp)) * comp
        ) for ts in data_matrix] for comp in components]

    print 'Creating PSDs of components...'
    comp_PSDs = [psd_and_bin(ts, data['sample_rate'], BPF_LO, BPF_HI
        )[:2] for ts in components]

    print 'Creating PSDs of data...'
    data_PSDs = [psd_and_bin(data['channels'][c]['nohwps'], data['sample_rate'],
        BPF_LO, BPF_HI)[:2] for c in channels]

    weight = 30
    for i in range(len(amplitudes)):
        outliers = amplitudes[i] > weight*np.median(amplitudes[i])
        num_outliers = np.sum(outliers)
        if num_outliers:
            plotComponent(i, '', comp_PSDs)
            for c in channels[outliers]:
                plotChannel(c, '', data_PSDs, labels)

def addargs(parser):
    # Default arguments
    parser.add_argument('-o','--output',dest='output_folder',help='Output folder')
    parser.add_argument('-f',dest='in_path',help='Input file',required=True)
    parser.add_argument('-v','--verbose',dest='verbose',action='store_true',help='Verbose output')
    parser.add_argument('--nocache',dest='cache',action='store_false',help='Only run if output file does not exist')

    # My arguments
    parser.add_argument('--gain_folder',help='Gain folder',required=True)
    parser.add_argument('--psd_folder',help='Statistics input folder',required=True) # Loads the PSD as part of the cw
    parser.add_argument('--wafer',help='Only map this wafer')

    # Time stream filter
    libmap_whwp.add_filter_timestream_args_whwp(parser)

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
