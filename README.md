#### Scripts

`case.py` - Main PCA script for investigating the anomaly

`eigval.py` - Initial PCA script



#### Data

`data.pkl` - PSD data

`data_ts.pkl` - Timestream data



#### Plot directories

`eigenGrids3/` - Same as `eigenGrids2/` but with BPF from 8.0 to 15.0 Hz

`eigenGrids2/` - Same as `eigenGrids1/` but with RMS calculation of eigenvector amplitudes instead of timestream normalization to magnitude 1

`eigenGrids1/` - PCA using band pass filtered timestream data (11.0 to 13.0 Hz)

`eigenGrids0/` - PCA of PSD data

`flags1/` - Flags of 8.2.0 wafer with PSD info (bit 30 off now)

`flags0/` - Flags of 8.2.0 wafer

`plots_hwp/` and `plots_nohwps/` - Initial testing of PCA code, with and without half wave plate signal



#### Miscellaneous

`PB1_hardwaremap.pkl` - Contains mapping of BoloIDs to channel numbers

`TarFileBuffer.py` - A utility to handle tar file buffer instead of making many small files

`test_hwpss.py` - Create PSD plots of the raw and HWPSS data

`draw_hdf5.py` - Draw an HDF5 image using `imshow()`

`example_libinspect.py` - Create some plots using libinspect.py
