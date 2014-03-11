import numpy as np
from scipy import stats

# read in light of flare flag files

# read in list of lightcurves that have flagged flares

for i in range(len(filelist)):

    # read in lightcurve & normalize to 1
    lightcurve = np.genfromtext(lightcurvefilename)
    time = lightcurve[:, 0]
    flux = lightcurve[:, 1]
    normflux = flux / np.median(flux)

    # read in associated flare flags - use getFlags?


    # compute descriptive stats (skew, kurtosis) for each flare 
    flareSkew = stats.skew
    flareKurt = stats.kurtosis

# compute 2nd derivative of points flagged as flare points

# test whether first and last points are lower than the midpoint

# compute the abs value of the slope around the flare

# compute the ratio of the abs value of the slopes on either side of the flare

# compute the amplitude of the full lightcurve (w stellar variabiity)

# compute the stddev of the flattened lightcurve 
# (w stellar variability subtracted)

