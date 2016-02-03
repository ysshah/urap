import libinspect
from matplotlib.pyplot import *

largepatch = '/scratch/ngoecknerwald/largepatch/'
hdf5_file = 'hdf5_ces_downsampled/20140803_232307/20140805_010231.hdf5'
gain_dir = 'gain/'

args=libinspect.make_args('--input %s --gain %s --detrend 2' %(largepatch+hdf5_file, largepatch+gain_dir))

tod=libinspect.load_TOD(args, channels=[0,1], do_PSD=True)

zoom_4f=np.logical_and(tod['channels'][0]['raw_freq'] > 7.8, tod['channels'][0]['raw_freq'] < 8.2)


figure()
plot(tod['channels'][0]['raw_freq'], tod['channels'][0]['raw_psd'], label='raw')
plot(tod['channels'][0]['freq_nohwps'], tod['channels'][0]['psd_nohwps'], label='hwpss')
yscale('log')
xlabel('f, Hz')
ylabel('PSD, $K^2 / Hz$')
title('Channel 0 HWPSS test')
legend()

figure()
plot(tod['channels'][0]['raw_freq'][zoom_4f], tod['channels'][0]['raw_psd'][zoom_4f], label='raw')
plot(tod['channels'][0]['freq_nohwps'][zoom_4f], tod['channels'][0]['psd_nohwps'][zoom_4f], label='hwpss')
yscale('log')
xlabel('f, Hz')
ylabel('PSD, $K^2 / Hz$')
title('Channel 0 HWPSS test, zoom on 4f line')
legend()

figure()
plot(tod['channels'][1]['raw_freq'], tod['channels'][1]['raw_psd'], label='raw')
plot(tod['channels'][1]['freq_nohwps'], tod['channels'][1]['psd_nohwps'], label='hwpss')
yscale('log')
xlabel('f, Hz')
ylabel('PSD, $K^2 / Hz$')
title('Channel 1 HWPSS test')
legend()

show()
