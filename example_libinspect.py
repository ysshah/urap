#!/usr/bin/python

from libinspect import *

#Prints out a list of all possible arguments
#make_args('-h')

# Example usage to load a data file with some filtering options
args = make_args('--input /scratch1/scratchdirs/ngoeckne/whwp/path_CES/20140725_005302/20140726_045109.hdf5 \
	--poly 15 --polpoly 5 --detrend 2 --nofft --deconvolve --gain_folder /scratch1/scratchdirs/ngoeckne/whwp/gain/')

#This loads the timestream for channel 102, make sure that what you pass here is a list
d = load_TOD(args, [102,])

#Prints all of the keys in the dictionary
#print d.keys()

#Find the first subscan
#begin,end=d['scanlist'][0]

#Take a subset of the TOD
begin,end=(0,679)

#Example plotting code
figure()
plot(d['time'][begin:end], d['whwpangle'][begin:end], label='WHWP Angle', color='b')
legend()
xlabel('Time, s')
ylabel('Angle')
title('WHWP Angle')
savefig('whwp_glitch.png')

figure()
plot(d['time'][begin:end], d['channels'][102]['filter_demod'][0][begin:end], label='I', color='b')
plot(d['time'][begin:end], d['channels'][102]['filter_demod'][1][begin:end], label='Q', color='r')
plot(d['time'][begin:end], d['channels'][102]['filter_demod'][2][begin:end], label='U', color='g')
legend()
ylabel('T RJ, K')
xlabel('Time, s')
title('Demodulated Timestream')
savefig('whwp_demodulated.png')

figure()
plot(d['time'][begin:end], d['channels'][102]['raw'][begin:end], label='Raw TOD', color='b')
legend()
ylabel('T RJ, K')
xlabel('Time, s')
title('Raw Timstream')
savefig('raw_timestream.png')

show()
