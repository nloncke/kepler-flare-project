import numpy as np
from scipy import stats
from lightcurves import *
import itertools

#def shift(l, n):
#    return l[n:] + l[:n]

def flareFeatures(files, flarfiles):
    """ Given files containing a list of lightcurve files and flag
    files, respectively, compute an array of features for each flare
    event
    """
    # generate list of flare flag files
    with open(flarfiles, 'r') as f:
        myFlags = [line.strip() for line in f]

    # generate list of lightcurve files
    with open(files, 'r') as f:
        myFiles = [line.strip() for line in f]

    flags = getflags(myFlags)
    
    for i in xrange(len(myFiles)):
        # read in lightcurve & normalize to 1
        lightcurve = np.genfromtxt(myFiles[i])
        time = lightcurve[:, 0]
        flux = lightcurve[:, 1]
        normflux = flux / np.median(flux)  # why median?
        
        # grab array of flagged indices
        lcflags = flags[i]
        floatflags = [float(j) for j in lcflags]

        # find boundaries bt events
        bnds = floatflags - np.roll(floatflags, 1)
        bnds = np.where(bnds != 1.)[0] # isolate individual flare event indices

        # compute relevant statistics for each event
        for j in xrange(len(bnds)):
            # get event indices into normflux array, NOT INCLUSIVE
            beg = lcflags[bnds[j]]
            if j == len(bnds) - 1:
                end = lcflags[-1] + 1
            else:
                end = lcflags[bnds[j+1] - 1] + 1

            flareSkew = stats.skew(normflux[beg:end])
            flareKurt = stats.kurtosis(normflux[beg:end])

            # compute 2nd derivative of points flagged as flare points
            firstDeriv = np.gradient(normflux[beg:end])
            secondDeriv = np.gradient(firstDeriv)

            # test whether first and last points are lower than the midpoint
            mid = int((end-beg) / 2)
            passedMidptCheck = True
            diff1 = normflux[mid] - normflux[beg]
            diff2 = normflux[mid] - normflux[end - 1]
            if diff1 <= 0 or diff2 <= 0:
                passedMidptCheck = False

            # compute the abs value of the slope around the flare

            # compute the ratio of the abs value of the slopes on either side of the flare

            # compute the amplitude of the full lightcurve (w stellar variabiity)


            # compute the stddev of the flattened lightcurve
            # (w stellar variability subtracted)

    return flareFeatureArray
