#!/usr/bin/env python
"""
Timestream analysis of POLARBEAR data.

@author Yash Shah
"""
from matplotlib.pyplot import *
import numpy as np, pickle, libinspect


def PCAplots(time, components):
	"""Plot the first few principle components."""
	figure()
	title('First 3 Principle Components')
	xlabel('Time')
	for i in range(3):
		plot(time, components[i], label='Component %d' %(i+1))
	legend()

	figure()
	title('Second 3 Principle Components')
	xlabel('Time')
	for i in range(3, 6):
		plot(time, components[i], label='Component %d' %(i+1))
	legend()

	show()


def PSDplots(components, data):
	figure()
	for i in range(len(amplitudes)):
		print "Plotting PSDs of component %d" %(i+1)
		title('PSD of Component %d' %(i+1))
		xlabel('Frequency')
		freq, psd, ones = psd_and_bin(components[i], data['sample_rate'], 0, 4)
		plot(freq, psd)
		yscale('log')
		margins(x=0.02, y=0.02)
		savefig('psd%d.png' %(i+1), dpi=200)
		clf()


def ampPlots(hwMap, ch, amplitudes):
	labels = [hwMap['boloid'][c] for c in ch]
	figure()
	for i in range(len(amplitudes)):
		print "Plotting eigenvector amplitudes of component %d" %(i+1)
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
		savefig('amp%d.png' %(i+1), dpi=200)
		clf()


def sanityCheck(components):
	"""
	Plot covariance matrix of principle components to confirm orthogonality.
	"""
	figure()
	imshow(abs(np.cov(components)), interpolation='none',
		norm=matplotlib.colors.LogNorm())
	colorbar()
	show()


if __name__ == '__main__':

	# Get light channels
	with open('lights.txt') as f:
		lights = eval(f.readline())

	# Get channels at center
	centerNums = range(25, 29) + range(34, 39) + range(44, 49
			 ) + range(54, 59) + range(64, 68)
	centerIDs = ['10.2_%dt' %(v) for v in centerNums]
	hwMap = pickle.load(open('PB1_hardwaremap.pkl'))
	bIDs = np.array(hwMap['boloid'])
	ch = np.argwhere(np.array([bIDs == v for v in centerIDs]))[:,1]
	ch = [c for c in ch if c in lights]

	# Define data directories
	# hdf5 = '/scratch/sample_largepatch/20140701_031318.hdf5'
	# gain = '/scratch/sample_gain/'
	largepatch = '/scratch/ngoecknerwald/largepatch/'
	hdf5 = largepatch+'hdf5_ces_downsampled/20140803_232307/20140805_010231.hdf5'
	gain = largepatch+'gain/'

	# Load data
	d = libinspect.make_args('--input %s --gain %s --detrend 2 --nofft' %(hdf5, gain))
	## Check if channels are in "cw"
	# od, df, cw = load_TOD(d, full=True, channels=[])
	# ch = [c for c in ch if c in cw.channels]
	data = libinspect.load_TOD(d, channels=ch, do_PSD=True)

	# Do whatever with data
	time = data['time']
	# Make a 2D array of Channels x Data
	# m = np.array([data['channels'][c]['raw'] for c in ch])
	m = np.array([data['channels'][c]['nohwps'] for c in ch])
	# Covariance matrix
	cov = np.cov(m)
	# Compute eigenvalues and eigenvectors
	evals, evecs = np.linalg.eig(cov)

	# 2D array of Data x Channels
	mT = m.transpose()
	# Principle components of data
	components = [np.dot(mT, evecs[:,i]) for i in range(10)]

	# amplitudes = [[np.dot(comp, data['channels'][c]['raw']
		# ) for c in ch] for comp in components]
	amplitudes = [[np.dot(comp, data['channels'][c]['nohwps']
		) for c in ch] for comp in components]

	# PCAplots(time, components)

	# PSDplots(components, data)

	# ampPlots(hwMap, ch, amplitudes)

	


	figure()
	# chOI = 987
	# plot(time, data['channels'][chOI]['nohwps'], label='NoHWPS channel '+bIDs[chOI])
	# plot(time, data['channels'][chOI]['raw'], label='Raw channel '+bIDs[chOI])

	for c in ch:
		plot(time, data['channels'][c]['nohwps'], label='Channel '+bIDs[c])
	plot(time, components[0], label='First component')
	title('All Channels with no HWPS')
	legend()
	show()
