#!/usr/bin/env python

'''

Code to examine WHWP timestreams, loads TOD into ipython

@Author Neil Goeckner-Wald, ngoecknerwald@berkeley.edu

'''

import numpy as np
import argparse
from AnalysisBackend.whwp import libmap_whwp
from AnalysisBackend.misc import numerical_analysis
from AnalysisBackend.deglitching import flaghdf5
from AnalysisBackend.pointing import pointingwrap
from AnalysisBackend.mapping import flatmap
from AnalysisBackend.deglitching import packetdrop
from AnalysisBackend.hdf5 import libautohdf5

##############################################################################
#
# Arguments: xs is frequency, ys is PSD values, k is initial bin size, r is growth rate
#
# Returns a log binned PSD, nmodes is number of frequency samples in a given bin
#
##############################################################################

def bin_PSD(xs,ys,k,r):
        '''
        k:      Initial bin size
        r:      Bin growth rate
        '''

        xb = []
        yb = []
        nmodes = []

        b0 = 0
        b1 = k
        binwidth = float(k)
        while True:
                if b0 >= len(xs): break

                x = np.mean(xs[b0:b1])
                y = np.mean(ys[b0:b1])

                xb.append(x)
                yb.append(y)
                nmodes.append(len(xs[b0:b1]))

                binwidth = binwidth*r
                b0 = b1
                b1 = b1 + int(binwidth)

        xb = np.array(xb)
        yb = np.array(yb)
        nmodes = np.array(nmodes)

        assert xb.shape == yb.shape and nmodes.shape == yb.shape

        return (xb, yb, nmodes)
################################################################
#
# Arguments: ts is timestream values as an array, sample rate is sampling rate of the timestream
#
# Band min and band max are the frequency limits to 
#
# Returns a PSD in K^2 s as a single sided PSD
#
################################################################

def psd_and_bin(ts, sample_rate, band_min, band_max):

	freq, psd = numerical_analysis.psd_nochunk(ts, Fs=sample_rate, zeropad=True)				
	ok = np.logical_and(freq < band_max, freq > band_min)

	#Check our sanity that we are actually returning something
	assert sum(ok) > 0

	freq = freq[ok]
	psd = psd[ok]

	return freq, psd, np.ones(psd.shape)

################################################################
#
# Main method, returns a python dictionary from loading data
#
# Filtering, gain parameters are specified by the argparse argument args
#
# List channels, if none specified it returns every active channel on the focal plane
#
# full=True specifies that the Calwrap and DataFolderhdf5 objects are also returned along with the python dictionary
#
# no_demod specifies that no demodulated timestreams will be returned, only the raw unfiltered timestreams. This makes the code dramatically faster
#
# do_PSD returns a PSD along with every timestream. This is slow!
#
# gainmode sets the calibration to use, 'RJ' is Rayleigh-Jeans temperature, 'CMB' is in CMB temperature units, None returns raw ADC counts
#
################################################################

def load_TOD(args, channels=None, full=False, no_demod=False, do_PSD=False, gainmode='RJ'):

	#Load the data
	df,cw = libmap_whwp.load_data_demod(args,open_mode='r',gain_mode=gainmode,require_hwp=False,nopsd=True)
	
	print 'Data folder keys: '	
	print df.__dict__.keys()
	print 'Calwrap keys: '
	print cw.__dict__.keys()

	if not cw.gains:
		print('No gains found! Be careful here...')
		return None

	print('Successfully loaded data folder and calwrap')

	#This section of code computes boresight pointing for the purpose of masking off point sources
	detpnt = pointingwrap.BoresightPointing(df.bolo_time,df.az,df.el,cw)
	detpnt.args = args
	detpnt.args.nofieldrot=False
	detpnt.args.azel=False
	detpnt.ra,detpnt.dec,detpnt.pa = pointingwrap.azel2radec_slalib.azel2radecpa(detpnt.t,detpnt.cw.ut1utc,detpnt.az,detpnt.el)
	detpnt.ra_src = np.median(np.unwrap(detpnt.ra))
	detpnt.dec_src = np.median(np.unwrap(detpnt.dec))
	detpnt.q = pointingwrap.quat_pointing.offset_radecpa_makequat(detpnt.ra,detpnt.dec,detpnt.pa,detpnt.ra_src,detpnt.dec_src)
	ra,dec,pa = pointingwrap.quat_pointing.offset_radecpa_applyquat(detpnt.q,0.,0.)
	xmax= max(np.abs(ra).max(),np.abs(dec).max())+0.025#enough radius to cover focal plane
	pixelsize = 2.*np.pi/180./60.
	mapinfo = flatmap.FlatMapInfo(-xmax,xmax,-xmax,xmax,pixelsize)
	args.mask_center = False
	sourcemaskmap,sourcemaskmap_nocoadd = libmap_whwp.mask_sources(mapinfo,detpnt.ra_src,detpnt.dec_src,args,df.bolo_time)

	#Make a channel list if none is provided
	if channels is None:
		channels = cw.channels

	# Find first and last sample that are inside the science scan
	i0 = df.scanlist[0][0]
	i1 = df.scanlist[-1][0] + df.scanlist[-1][1]
	nt = df.bolo_time.size
	obsmask = np.zeros(nt,dtype=bool)
	obsmask[i0:i1]=True
	az = df.az
	el = df.el
	whwpangle=df.whwpangle
	scanlist=df.scanlist
	
	chdict = {}
	
	for ch in channels:
		
		boloname = cw.xmlmap.boloid[ch]
		print('Loading channel: %d, boloname: %s' % (ch, boloname))

		#Do things with dark channels
		if (not ch in cw.channels) or no_demod:
			d={}
			print('Unrecognized channel, probably a dark bolometer!')
			ts = libmap_whwp.get_tsblock(df,cw,[ch],verbose=True,keepDark=True)
			if do_PSD:	
				psdinfo=psd_and_bin(ts, df.sample_rate, 0, 4.)
				d['raw_freq'] = psdinfo[0]
                        	d['raw_psd'] = psdinfo[1]
                        	d['raw_nmodes'] = psdinfo[2]
			d['raw'] = ts
			chdict[ch] = d
			continue

		#do things with live channels
		(waferi1d,),waferpa,waferbeamorient = libmap_whwp.getwaferi1dpa_channel([mapinfo],detpnt,cw,[ch])
		wafermask_chan,wafermask_chan_filt = libmap_whwp.build_wafermask_channel(waferi1d,obsmask,sourcemaskmap, \
			sourcemaskmap_nocoadd,[ch],None,cw)
		drop = packetdrop.merge_packetdrops(packetdrop.construct_packetdrop_list(wafermask_chan_filt[0]))

		#load in the raw timestream
		ts = libmap_whwp.get_tsblock(df,cw,[ch],verbose=False)
		
		#demodulate the timestream
		ts_demod = libmap_whwp.demodulate_ts(args,np.copy(ts),[ch],df)
		libmap_whwp.calibrate_angle(ts_demod, cw, [ch])

		if not df.demoded:
                       	df.dm.chmask = wafermask_chan_filt[0]	
		for i in range(3):
			packetdrop.fill_packetdrop(ts_demod[0][i],drop)

		#Filter the timstreams
		filterinfo = libmap_whwp.filter_timestreams_demod(ts_demod,wafermask_chan,wafermask_chan_filt,[ch],args,df,cw)

		ts = ts[0]			
		ts_demod=ts_demod.reshape(3,nt)

		d={}

		if do_PSD:      
			#Make a PSD for the raw non-demodulated data
			psdinfo=psd_and_bin(ts, df.sample_rate, 0, 16.)
			d['raw_freq'] = psdinfo[0]
			d['raw_psd'] = psdinfo[1]
			d['raw_nmodes'] = psdinfo[2]

			#Compute power spectral density
        	        freqs_demod = []
        	        psds_demod = [] # in K^2 per Hz
        	        nmodes_demod = []

			for i in range(0, 3):
				psdinfo=psd_and_bin(ts_demod[i], df.sample_rate, 0, 4.)			
				freqs_demod.append(psdinfo[0])
				psds_demod.append(psdinfo[1])
				nmodes_demod.append(psdinfo[2])

			d['freqs_demod'] = freqs_demod
                	d['psds_demod'] = psds_demod
                	d['nmodes_demod'] = nmodes_demod

		#Save the output where we can find it
		d['filterinfo'] = filterinfo
		d['raw'] = ts
		d['filter_demod'] = ts_demod
		chdict[ch] = d

	outdict={}
	outdict['channels'] = chdict
	outdict['az'] = az
        outdict['el'] = el
        outdict['whwpangle'] = whwpangle
        outdict['scanlist'] = scanlist
	outdict['sample_rate'] = df.sample_rate
	outdict['time'] = np.asarray(range(az.shape[0]))/df.sample_rate
        outdict['obsmask'] = obsmask

	if full:
		return outdict, df, cw

	return outdict


###########################################################
#
# This returns the argparse object needed for the load_TOD method
#
# To see a list of available options run this method with the default '-h'
#
# Sample args_param may be '--input /scratch1/scratchdirs/ngoeckne/whwp/hdf5_ces_downsampled_10days/\
# 20140725_005302/20140731_142250.hdf5 --gain /scratch1/scratchdirs/ngoeckne/whwp/gain/'
#
###########################################################

def make_args(args_param='-h'):
	
	parser = argparse.ArgumentParser()
        parser.add_argument('--nocache',dest='cache',action='store_false',help='Only run if output file does not exist')
        parser.add_argument('--turnarounds',dest='turnarounds',action='store',type=float,help='Include turnarounds in TS masking. Put 1 for all',default=None)
        parser.add_argument('--require_pair',action='store_true',help='Require both bolos work')
        parser.add_argument('--input',dest='in_path',action='store',help='Input path',required=True)
        parser.add_argument('--gain_folder',action='store',help='Gain folder')
        parser.add_argument('--mask_radius', dest='mask_radius', action='store', type=float, help='Mask radius for point sources (arcmin) for filtering', default=5.0)
        parser.add_argument('--maskephem',dest='maskephem',action='append',default=[],help='Mask and do not map source by using ephemeris file.')
        parser.add_argument('--maskephemradius',dest='maskephemradius',type=float,default=15.0,action='store',help='Mask radius to use for maskephem sources')
        libmap_whwp.add_filter_timestream_args_whwp(parser)

	try:
		return parser.parse_args(args_param.split())
	except SystemExit:
		return None



# End of library
