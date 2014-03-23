import numpy as np
from scipy import stats
from lightcurves import *
import itertools

#def shift(l, n):
#    return l[n:] + l[:n]

def window(beg, end, width, maxindex):
    """Generates a window around beg:end with at most +/-width points on
    either side.  Takes the maximum index of the array to be indexed into
    for bounds checking.  Returns tuple (newbeg, newend) or -1 if invalid
    arguments.
    """
    if (beg < 0) or (end >= length):
        return -1

    maxwidth = min(beg, maxindex-end)
    real_width = min(width, maxwidth)
    return (beg-real_width, end+real_width)


def flareFeatures(files, flarfiles):
    """ Given files containing a list of lightcurve files and flag
    files, respectively, compute an array of features for each flare
    event.  Output is a list of dictionaries, structured as follows:
    [...,
        {'id':, 'num_events':, 'amplitude':, 'stddev':, 'flare_features':},
    ...]

    The value of the flare_features key is a list of dictionaries (one for each event):
    [...,
        {'skew':, 'kurtosis':, 'second_deriv':, 'passed_midpt_check':, 'slope':, 'slope_ratio':},
    ...]
    """

    # generate list of flare flag files
    with open(flarfiles, 'r') as f:
        myFlags = [line.strip() for line in f]

    # generate list of lightcurve files
    with open(files, 'r') as f:
        myFiles = [line.strip() for line in f]

    flags = getflags(myFlags)
    flareFeatureArray = list()
    
    # gather data on events and create dictionary for each lightcurve
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

        ltcurve_dict = dict()
        # ltcurve_dict['id'] = 
        ltcurve_dict['num_events'] = len(bnds)
        ltcurve_dict['flare_features'] = list()

        # collect flare-specific data and create dictionary for each event 
        for j in xrange(len(bnds)):
            # get event indices into normflux array, ready for slicing
            beg = lcflags[bnds[j]]
            if j == len(bnds) - 1:
                end = lcflags[-1] + 1
            else:
                end = lcflags[bnds[j+1] - 1] + 1

            flareSkew = stats.skew(normflux[beg:end])
            flareKurt = stats.kurtosis(normflux[beg:end])

            # compute 2nd derivative of points flagged as flare points
            first_deriv = np.gradient(normflux[beg:end])
            second_deriv = np.gradient(first_deriv)

            # test whether first and last points are lower than the midpoint
            mid = int((end-beg) / 2)
            passed_midpt_check = True
            diff1 = normflux[mid] - normflux[beg]
            diff2 = normflux[mid] - normflux[end - 1]
            if diff1 <= 0 or diff2 <= 0:
                passed_midpt_check = False

            # compute the abs value of the slope around the flare


            # compute the ratio of the abs value of the slopes on either side of the flare
            pre_slope = (normflux[beg-1]-normflux[beg-5]) / (time[beg-1]-time[beg-5])
            post_slope = (normflux[end+4]-normflux[end]) / (time[end+4]-time[end])
            slope_ratio = np.abs(post_slope/pre_slope)

            # add entry to flare_features list
            event_dict = dict()
            event_dict['skew'] = flareSkew
            event_dict['kurtosis'] = flareKurt
            event_dict['second_deriv'] = second_deriv
            event_dict['passed_midpt_check'] = passed_midpt_check
            # event_dict['slope'] = 
            event_dict['slope_ratio'] = slope_ratio
            ltcurve_dict['flare_features'].append(event_dict)


        # compute the amplitude of the full lightcurve (w stellar variabiity)
        amp = np.max(normflux) - np.min(normflux)

        # compute the stddev of the flattened lightcurve
        # (w stellar variability subtracted)

        # update values in flare dictionary and add to flareFeatureArray
        ltcurve_dict['amplitude'] = amp
        # ltcurve_dict['stddev'] = stddev
        flareFeatureArray.append(ltcurve_dict)

    return flareFeatureArray
