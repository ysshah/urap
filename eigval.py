from libinspect import *
from matplotlib.pyplot import *
import numpy as np, pickle, socket

if __name__ == '__main__':

	# Get light channels
	with open('lights.txt') as f:
		lights = eval(f.readline())

	# Get channels at center
	centerNums = range(25, 29) + range(34, 39) + range(44, 49) + range(54, 59) + range(64, 68)
	centerIDs = ['10.2_%dt' %(v) for v in centerNums]
	hwMap = pickle.load(open('PB1_hardwaremap.pkl'))
	bIDs = np.array(hwMap['boloid'])
	ch = np.argwhere(np.array([bIDs == v for v in centerIDs]))[:,1]
	ch = [c for c in ch if c in lights]

	# Define data directories
	if socket.gethostname().startswith('gordita'):
		hdf5 = '/scratch/sample_largepatch/20140701_031318.hdf5'
		gain = '/scratch/sample_gain/'
	else:
		dataDir = '/scratch1/scratchdirs/takakura/largepatch/'
		hdf5 = dataDir + 'hdf5_ces_downsampled/20140707_025357.hdf5'
		gain = dataDir + 'gain/'

	# Load data
	d = make_args('--input %s --detrend 2 --nofft --gain %s' %(hdf5, gain))
	## Check if channels are in "cw"
	# od, df, cw = load_TOD(d, full=True, channels=[])
	# ch = [c for c in ch if c in cw.channels]
	data = load_TOD(d, channels=ch, do_PSD=True)

	# Do whatever with data
	time = data['time']
	# Make a 2D array of Channels x Data
	m = np.array([data['channels'][c]['raw'] for c in ch])
	# Covariance matrix
	cov = np.cov(m)
	# Compute eigenvalues and eigenvectors
	evals, evecs = np.linalg.eig(cov)

	# scatter(ch, evals)
	# plot(ch, evals)
	# yscale('log')
	# show()

	# 2D array of Data x Channels
	mT = m.transpose()

	components = [np.dot(mT, evecs[:,i]) for i in range(10)]

	freq, psd, ones = psd_and_bin(components[0], data['sample_rate'], 0, 4)

	figure()
	plot(freq, psd)
	yscale('log')
	show()
