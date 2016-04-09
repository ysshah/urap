#!/usr/bin/env python
"""
Flag information to JSON.

@author Yash Shah
"""
import numpy as np, pickle, os, sys, json

DATE = '20140921_010507'
CUTS = '/scratch/ngoecknerwald/resonance_timestreams/alpha_cut_yash/%s.pkl' %DATE
cut = pickle.load(open(CUTS))
flags = cut[0]
hwMap = pickle.load(open('PB1_hardwaremap.pkl'))
DEST = 'flags.json'

tChan = np.array([hwMap['boloid'].index('8.2.0_%dt' %i) for i in range(1,92)])
bChan = np.array([hwMap['boloid'].index('8.2.0_%db' %i) for i in range(1,92)])
allChannels = np.concatenate((tChan, bChan))


def flagsToJSON():
    d = {"date": DATE, "file": CUTS}
    for i in range(len(cut[1])):
        d["bit%d"%i] = [a.item() for a in (
            (flags[allChannels][:,1] & (1 << i)) >> i).astype(int)]
    for i,c in enumerate(cut[1]):
        d["desc%d"%i] = c
    json.dump(d, open(DEST, 'w'))


if __name__ == '__main__':
    main()
