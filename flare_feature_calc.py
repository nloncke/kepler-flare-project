from sklearn import svm
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


def feat_dict_to_bunch(flareFeatureArray, vetfile=None):
    """ Given a list of dictionaries with the following structure,
    output a bunch of flare features ready for scikitlearn.  Output
    has dimensions num_samples x num_features.  Each sample is an
    individual flare, each column is a feature (statistic).  If
    vetfile is provided, assumes that all the flares in the feature
    array have been vetted with responses recorded in vetfile.

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

    bunch = dict()  # the entire bundle 
    data = list()   # data array of features
    if vetfile: # we got labelled data
        all_vets = vetfile_to_dict(vetfile)
        target = list()

    # collect flare-wide features first
    for curve in flareFeatureArray:
        curvespecs = [curve["amplitude"], curve["num_events"], curve["stddev"]]

        # now visit each flare within each lightcurve
        for flare in curve["flare_features"]:
            flarespecs = list(curvespecs)   # copy and add to template
            feats = [int(flare["has_consec_points"]), flare["kurtosis"],
                     int(flare["passed_midpt_check"]), flare["second_deriv"],
                     flare["skew"], flare["slope"], flare["slope_ratio"]]
            flarespecs.extend(feats) # join flare with curve template
            data.append(flarespecs)  # add flare data to array

        if vetfile: # if labelled data
            default = ['m'] * curve["num_events"]
            vet = all_vets.get(curve["id"], default) # array of ynm for the lightcurve
            target.extend(vet)

    # package things nicely
    bunch["data"] = np.array(data)
    bunch["feature_names"] = ["amplitude","num_events","stddev",
                              "has_consec_points","kurtosis",
                              "passed_midpt_check", "second_deriv",
                              "skew", "slope", "slope_ratio"]
    if vetfile:
        bunch["target"] = vets_to_ints(target)
        bunch["target_names"] = np.array(['n','y','m'])

    return bunch


def vetfile_to_dict(vetfile):
    """Given vetfile in the format specified in lightcurves.py, this
    function converts the responses to a dictionary with Kepler IDs as
    keys and lists of y/n/m as values.
    """
    responses = dict()
    with open(vetfile, 'r') as f:
        for line in f:
            tokens = line.split()
            if tokens[0].isdigit():
                responses[tokens[0]] = tokens[1:]
            else:
                # is there a way to issue a warning and continue with
                # processing the rest of the lines in the file?
                raise ValueError('{} is not a valid Kepler id.'.format(tokens[0]))
    return responses


def vets_to_ints(target):
    """ Prepares target array of ynm to integers:
    n -> 0
    y -> 1
    m -> 2
    """
    vetmap = {'n': 0, 'y': 1, 'm': 2}
    mapped = list(target)
    for i in xrange(len(mapped)):
        mapped[i] = vetmap.get(mapped[i])
    return np.array(mapped)


def overlap(dictlist, vetfile):
    """ Collects the dictionaries in dictlist containing key-value
    pairs id: kid, where 'kid' is the first token per line in vetfile.
    """
    with open(vetfile, 'r') as f:
        kids = [line.split()[0] for line in f]
    return [item for item in dictlist if item["id"] in kids]


def run_classifier(files, flarefiles, vetfile):
    """ Script for training classifier and testing on training inputs
    """
    allFeats = flareFeatures(files, flarefiles)
    data = overlap(allFeats, vetfile)
    bunch = feat_dict_to_bunch(data, vetfile)

    # train classifier, test on labelled data
    clf = svm.SVC()
    clf.fit(bunch["data"], bunch["target"])
    predictions = clf.predict(bunch["data"])
    hits = predictions == bunch["target"]
    return sum(hits) / float(len(hits))
