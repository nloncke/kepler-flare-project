from scipy import stats
from lightcurves import *
import numpy as np
import string
import scipy.interpolate as sp

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
    """ Given a single file containing a list of lightcurve files (arg
    files) and flag files (arg flarefiles), respectively, compute an
    array of features for each flare event.  Output is a list of
    dictionaries, structured as follows:

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

        # update values in flare dictionary
        ltcurve_dict = dict()
        ltcurve_dict['id'] = kid
        ltcurve_dict['num_events'] = num_events
        ltcurve_dict['amplitude'] = amp

        # collect flare-specific data and create dictionary for each event 
        ltcurve_dict['flare_features'] = list()
        ignore_in_smoothed = set() # collects pts to be ignored when smoothing
        for j in xrange(len(bnds)):
            # get event indices into normflux array, ready for slicing
            beg = lcflags[bnds[j]]
            if j == len(bnds) - 1:
                end = lcflags[-1] + 1
            else:
                end = lcflags[bnds[j+1] - 1] + 1

            # compute window around flare for future operations
            wind_width = 4 
            [wind_beg, wind_end] = window(beg, end, wind_width, len(flux))

            # compute window around flare to be deleted when smoothing lc
            [ignore_beg, ignore_end] = window(beg, end, 2, len(flux))
            ignore_in_smoothed = ignore_in_smoothed.union(set(range(ignore_beg,ignore_end)))
            
            # skew, kurtosis, mean of 2nd deriv around window
            flareSkew = stats.skew(normflux[wind_beg:wind_end])
            flareKurt = stats.kurtosis(normflux[wind_beg:wind_end])
            first_deriv = np.gradient(normflux[wind_beg:wind_end])

            # around the flare (ie, before and after like slope)
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

        ######### EXIT THE FLARE LOOP #############
        ######### ENTER LIGHTCURVE SCOPE ##########

        # compute the stddev of the flattened lightcurve (w stellar variability subtracted)
        # take out a few points on either side for better interpolation
        cutT = np.delete(time, list(ignore_in_smoothed))
        cutFlux = np.delete(normflux, list(ignore_in_smoothed))
        fflux = sp.interp1d(cutT, cutFlux, kind='nearest') # flux function
        fsmoothed = smooth(fflux(time), window='flat')
        stddev = stats.tstd(fsmoothed)
        
        # add data for current lightcurve to array
        ltcurve_dict['stddev'] = stddev
        flareFeatureArray.append(ltcurve_dict)

    return flareFeatureArray


def dict_to_arr(flareFeatureArray):
    """ Given a list of dictionaries with the following structure,
    output an array of flare features ready for scikitlearn.  Each row
    is a separate flare, each column is a feature.

    Input
    -----
    [..., {'id':, 'num_events':, 'amplitude':, 'stddev':,
        'flare_features':}, ...]

    The value of the flare_features key is a list of dictionaries (one
    for each event):

    [..., {'skew':, 'kurtosis':, 'second_deriv':, 'slope':,
    'slope_ratio':, 'passed_midpt_check':, 'has_consec_points':}, ...]

    Output
    ------
    [stuff]
    """

    data = list()

    # add number of flares per lightcurve to the features
    # collect flare-wide features first
    for curve in flareFeatureArray:
        curvespecs = [curve["amplitude"], curve["stddev"], curve["num_events"]]

        # now visit each flare within each lightcurve
        for flare in curve["flare_features"]:
            flarespecs = list(curvespecs)   # copy and add to template
            flarespecs.append(int(flare["has_consec_points"])) # bool
            flarespecs.append(flare["kurtosis"])
            flarespecs.append(int(flare["passed_midpt_check"])) # bool
            flarespecs.append(flare["second_deriv"])
            flarespecs.append(flare["skew"])
            flarespecs.append(flare["slope"])
            flarespecs.append(flare["slope_ratio"])

            data.append(flarespecs) # add flare data to array_midpt

    return np.array(data)
