import numpy as np
from scipy import stats
from lightcurves import *
import itertools
import string
#def shift(l, n):
#    return l[n:] + l[:n]

def window(beg, end, width, maxindex):
    """Generates a window around beg:end with at most +/-width points
    on either side.  Takes the maximum index of the array to be
    indexed into for bounds checking.  Returns tuple (newbeg, newend)
    or -1 if invalid arguments.
    """
    if (beg < 0) or (end >= maxindex):
        return -1

    maxwidth = min(beg, maxindex-end)
    real_width = min(width, maxwidth)
    return (beg-real_width, end+real_width)


def flareFeatures(files, flarfiles):
    """ Given files containing a list of lightcurve files and flag
    files, respectively, compute an array of features for each flare
    event.  Output is a list of dictionaries, structured as follows:

    [..., {'id':, 'num_events':, 'amplitude':, 'stddev':,
        'flare_features':}, ...]

    The value of the flare_features key is a list of dictionaries (one
    for each event):

    [..., {'skew':, 'kurtosis':, 'second_deriv':, 'slope':,
    'slope_ratio':, 'passed_midpt_check':, 'has_consec_points':}, ...]
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
        # strip non-digit chars from filename
        alltab=string.maketrans('','')
        nodigits=alltab.translate(alltab, string.digits)
        kid = myFiles[i].translate(alltab, nodigits)

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
        num_events = len(bnds)

        # compute the amplitude of the full lightcurve (w stellar variabiity)
        amp = np.max(normflux) - np.min(normflux)

        # compute the stddev of the flattened lightcurve
        # (w stellar variability subtracted)

        # update values in flare dictionary
        ltcurve_dict = dict()
        ltcurve_dict['id'] = kid
        ltcurve_dict['num_events'] = num_events
        ltcurve_dict['amplitude'] = amp
        # ltcurve_dict['stddev'] = stddev

        # collect flare-specific data and create dictionary for each event 
        ltcurve_dict['flare_features'] = list()
        for j in xrange(len(bnds)):
            # get event indices into normflux array, ready for slicing
            beg = lcflags[bnds[j]]
            if j == len(bnds) - 1:
                end = lcflags[-1] + 1
            else:
                end = lcflags[bnds[j+1] - 1] + 1

            # compute (somewhat arbitrary) window around flare
            wind_width = 4 
            [wind_beg, wind_end] = window(beg, end, width, len(flux))

            # skew, kurtosis, mean of 2nd deriv around window
            flareSkew = stats.skew(normflux[wind_beg:wind_end])
            flareKurt = stats.kurtosis(normflux[wind_beg:wind_end])
            first_deriv = np.gradient(normflux[wind_beg:wind_end])
            second_deriv = np.mean(np.gradient(first_deriv))

            # compute the abs value of the slope around the flare
            slope = np.abs((normflux[wind_end]-normflux[wind_beg]) /
                           (time[wind_end]-time[wind_beg]))

            # compute the ratio of the abs value of the slopes on either side of the flare
            pre_slope = (normflux[beg]-normflux[wind_beg]) / (time[beg]-time[wind_beg])
            post_slope = (normflux[wind_end]-normflux[end]) / (time[wind_end]-time[end])
            slope_ratio = np.abs(post_slope / pre_slope)

            # test whether first and last points are lower than the midpoint
            mid = int((end-beg) / 2)
            passed_midpt_check = True
            diff1 = normflux[mid] - normflux[beg]
            diff2 = normflux[mid] - normflux[end - 1]
            if diff1 <= 0 or diff2 <= 0:
                passed_midpt_check = False

            # test whether flagged points are evenly spaced
            has_consec_points = True
            diffs = (time[beg:end] - np.roll(time[beg:end],1))[1:]
            if np.where(diffs > 0.5)[0]:
                has_consec_points = False

            # add entry to flare_features list
            event_dict = dict()
            event_dict['skew'] = flareSkew
            event_dict['kurtosis'] = flareKurt
            event_dict['second_deriv'] = second_deriv
            event_dict['slope'] = slope
            event_dict['slope_ratio'] = slope_ratio
            event_dict['passed_midpt_check'] = passed_midpt_check
            event_dict['has_consec_points'] = has_consec_points
            ltcurve_dict['flare_features'].append(event_dict)

        # add data for current lightcurve to array
        flareFeatureArray.append(ltcurve_dict)
    return flareFeatureArray
