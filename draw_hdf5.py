from AnalysisBackend.mapping import maprep
from matplotlib.pyplot import *
import matplotlib
import sys
import numpy as np

#define constants
radToDeg = 57.2958
MAP_PIX = 271

print('Usage: python draw_hdf5.py <map to be plotted>.hdf5')

#read in the file
filename = sys.argv[1]
m = maprep.MapMakingVectors.load(filename)

#read the patch description from the mapmaker
source_name = m.mapinfo.source
pixel_size = m.mapinfo.pixel_size*radToDeg

if(source_name == 'PB1RA23HAB'):
	ra_center = 23.03
	dec_center = -32.8
else:
	ra_center = 0.
	dec_center = 0.
ra_center *= 15.


print('The pixel size is:' + str(pixel_size))

#read the information from the hdf5 file
nHits_raw = m.mapinfo.view2d(m.nhit).T
I_raw = m.mapinfo.view2d(m.I).T
Q_raw = m.mapinfo.view2d(m.Q).T
U_raw = m.mapinfo.view2d(m.U).T
I_weight = m.mapinfo.view2d(m.w).T

def downsample_map(data, N=8):
	from numpy import average, split
	width = data.shape[0]
	height= data.shape[1]
	width_used = width - width % N	
	height_used = height - height % N
	data_used = data[0:width_used, 0:height_used]

	return average(split(average(split(data_used, width_used // N, axis=1), axis=-1), height_used // N, axis=1), axis=-1) 

#downsampled maps to see if there is anything there
nHits = downsample_map(nHits_raw, N=1)
I = downsample_map(I_raw, N=1)
Q = downsample_map(Q_raw, N=1)
U = downsample_map(U_raw, N=1)

#decide the translatiuon between pixels and the RA, DEC coordinates of the patches
ra_min = ra_center - (pixel_size*(float(MAP_PIX)/2.)/np.cos(dec_center / radToDeg))
ra_max = ra_center + (pixel_size*(float(MAP_PIX)/2.)/np.cos(dec_center / radToDeg))
dec_min = dec_center - pixel_size*(float(MAP_PIX)/2.)
dec_max = dec_center + pixel_size*(float(MAP_PIX)/2.)
map_extent=[ra_min, ra_max, dec_min, dec_max]

#Apply the mask where we actually have data
nHits_mask = np.ma.masked_where(nHits < 1, nHits)
weight_mask = np.ma.masked_where(nHits < 1, I_weight)
I_mask = np.ma.masked_where(np.abs(I) > 0.0005, I)
Q_mask = np.ma.masked_where(np.abs(Q) > 0.0005, Q)
U_mask = np.ma.masked_where(np.abs(U) > 0.0005, U)

#do all of the plotting
figure()
imshow(nHits_mask, cmap=cm.Blues, extent=map_extent, origin='lower')
colorbar()
title('nhits')

figure()
imshow(weight_mask, cmap=cm.Blues, extent=map_extent, origin='lower')
colorbar()
title('w')

figure()
imshow(I_mask, cmap=cm.jet, extent=map_extent, origin='lower')
colorbar()
title('I')

figure()
imshow(Q_mask, cmap=cm.jet, extent=map_extent, origin='lower')
colorbar()
title('Q')

figure()
imshow(U_mask, cmap=cm.jet, extent=map_extent, origin='lower')
colorbar()
title('U')

show()
