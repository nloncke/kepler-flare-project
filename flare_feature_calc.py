import numpy as np
from scipy import stats
from lightcurves import *
import itertools

#def shift(l, n):
#    return l[n:] + l[:n]

def flareFeatures(files, flarfiles)
    """ Given a list of lightcurve files and flag files, 
    compute an array of features for each flare event"""
    # read in list of flare flag files
    with open(flarfiles, 'r') as f:
        myFlags = [line.strip() for line in f]
    f.close()

    # read in list of lightcurve files
    with open(files, 'r') as f:
        myFiles = [line.strip() for line in f]
    f.close()

    flags = getflags(flarfiles)
    
    for i in range(len(files)):

        # read in lightcurve & normalize to 1
        lightcurve = np.genfromtext(files[i])
        time = lightcurve[:, 0]
        flux = lightcurve[:, 1]
        normflux = flux / np.median(flux)
        
        # grab array of flagged indices
        lcflags = flags[i]
        floatflags = [float(j) for j in lcflags]

        # find boundaries bt events
        bnds = floatflags - np.roll(floatflags, 1)
        # indices --> having some confusion about tuples vs indices here
        tst = np.where(bnds != 1.)

        # compute descriptive stats (skew, kurtosis) for each flare 
        flareSkew = stats.skew(normflux)
        flareKurt = stats.kurtosis(normflux)

        # compute 2nd derivative of points flagged as flare points



        # test whether first and last points are lower than the midpoint

        # compute the abs value of the slope around the flare

        # compute the ratio of the abs value of the slopes on either side of the flare

        # compute the amplitude of the full lightcurve (w stellar variabiity)

        # compute the stddev of the flattened lightcurve 
        # (w stellar variability subtracted)

    return flareFeatureArray
