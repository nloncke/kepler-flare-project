"""Collection of functions for plotting light curves from and
searching for stellar flares in kepler data.

Some conventions:

flarfiles refers to a list of strings with entries 'kid#_flares.txt'.
These files contain the indices of potential flare events.

files refers to a list of strings with entries 'kid#.txt'.  These
files hold the brightness data.

flags refers to a nested list of ints.  Axis 0: each light curve; Axis
1: the indices of flare events.
"""

import utils
import commands
import numpy as np
import scipy.interpolate as sp
import scipy.integrate as si
import matplotlib.pyplot as plt

def getflags(flarfiles):
    """Given a list of the names of the files holding the flare
    flags, outputs a nested list of the indices of potential flaring
    events.
    """
    for i in xrange(len(flarfiles)):
        # read indices from a kid*_flares.txt file
        indcs = np.genfromtxt(flarfiles[i], dtype=int).tolist()

        # insert indices into 2d array
        if i == 0:
            flags = indcs
        else:                
            flags[i] = indcs

        # create next slot
        if i == 0:
            flags = [flags, []]
        elif i == len(flarfiles) - 1:
            pass
        else:
            flags.append([])
    return flags

def ltcurve(filename, display=True, style='-'):
    """Plots light curves from a single kepler data file (col1 = time;
    col2 = brightness).  Returns the time and brightness data arrays.
    """
    data = np.genfromtxt(filename)
    days = data[:, 0]
    flux = data[:, 1]
    normflux = flux / np.mean(flux)

    if display:
        title = 'Light curve for {0}'.format(filename)
        plt.plot(days, normflux, style)
        plt.title(title)
        plt.xlabel(r'$t (days)$')
        plt.ylabel(r'$brightness$')
        plt.show()

    return days, normflux


def ltcurves(files, flags=None):
    """Plots the light curves from the input list of kepler data
    files.  If passed 2D nested list of indices, overplots the
    indicated data points.  Still need to figure out how to go back a
    plot, though.
    """
    for i in range(len(files)):
        data = np.genfromtxt(files[i])

        days = data[:, 0]
        flux = data[:, 1]
        normflux = flux/np.mean(flux)

        title = 'Light curve for {0}'.format(files[i])
        
        plt.plot(days, normflux, '-')
        plt.title(title)
        plt.xlabel(r'$t (days)$')
        plt.ylabel(r'$brightness$')

        # just look at the flags, no vetting
        if flags is not None:
            plt.plot(days[flags[i]], normflux[flags[i]], 'y*',
        markersize=10)

        plt.show()

        # wait for user input to see next light curve
        print('Hit return for the next light curve. ')

        ans = str(raw_input())
        if ans != None:
            plt.close()
    print("Wow, you're done.")
    return


def kid_idx(vetfile, files):
    """Reads last kid from vet output file and returns index from
    input list of files.
    """
    f = open(vetfile, 'r')
    content = f.readlines()
    kid = ''

    # check for kid starting from last line
    for i in reversed(range(len(content))):
        if kid != '': break                 # kid already found
        elif content[i][0].isdigit():       # get first char in line
            line = content[i]
            for j in range(len(line)):
                if not line[:j+1].isdigit(): # substr contains
                    # non-digit
                    kid = line[:j]
                    break
    # obtained kid, now where in files is it?
    if kid == '':
        print "File empty or not in standard format"
        return -1
    else:
        for i in range(len(files)):
            if files[i].find(kid) != -1:
                # print "Last vetted {0}, index {1} of
                # {2}".format(files[i], i, len(files) - 1)
                return i
    return -1                   # last vetted file isn't in list


def flareshow(files, start=0, flags=None, vet=False, output='output.txt'):
    """Overplots the potential flares on their light curves, given a
    list of the names of the files containing brightness data and a
    nested list of the event indices.  Begins displaying from
    files[start].  If vet=True, writes user input y/n/m to output file
    and returns a nested list of responses.
    """
    assess = list()

    # open file for recording vets
    if vet:
        f = open(output, 'a')  # file for y, m, n
        filename = output.partition('.')
        vidx = open(filename[0] + '_indices' + filename[1] +
    filename[2], 'a') # for vetted event indices

    for i in range(start, len(files)):
        # plot the lightcurve
        if vet:
            plt.subplot(211)    # leave room to zoom
        days, normflux = ltcurve(files[i], display=True, style='-')
        if flags is not None:
            # overplot the potential flares
            plt.plot(days[flags[i]], normflux[flags[i]], 'y*',
    markersize=10)
            plt.show()

            # nonzero() returns a tuple, hence the slice
            # the +1 accounts for the first contiguous interval
            events = np.nonzero(np.diff(flags[i]) - 1)[0]
            bins = len(events)
            num_events = bins + 1

        if vet:
            # print 'bins is {0}'.format(bins)
            plt.subplot(212)    # switch to zoom plot
            ltcurve(files[i], display=True, style='-')
            ltcurve(files[i], display=True, style='.')
            plt.plot(days[flags[i]], normflux[flags[i]], 'y*',
        markersize=10)

            # start new row in file
            kid = files[i].lstrip('kid').rstrip('.txt')
            if i == start:
                f.write('{0:<9}'.format(kid))
                vidx.write('{0:<9}'.format(kid))
            else:
                f.write('\n{0:<9}'.format(kid))
                vidx.write('\n{0:<9}'.format(kid))

            # what if there's only one event?!
            if bins == 0:
                print '{0}, light curve {1} of {2}: {3} event found'.format(files[i], i+1, len(files), num_events)
                end = -1        # last flag index
                xmargin = 6 * (days[flags[i][end]] - days[flags[i][0]])
                ymargin = np.max(normflux[flags[i]]) - np.min(normflux[flags[i]]) # diff between hi and lo
                plt.xlim(days[flags[i][0]] - xmargin, days[flags[i][end]] + xmargin)
                # avg = 0.5 * (np.max(normflux[flags[i][:end]])
                #              + np.min(normflux[flags[i][:end]])) #
                #              hi/lo avg = center
                plt.ylim(np.min(normflux[flags[i]]) - 2*ymargin, np.max(normflux[flags[i]]) + ymargin)

                # prompt and add response to pseudo-array
                # print "i is {0}, start is {1}".format(i, start)
                ans = str(raw_input('Does this look like a flare? Say y, n, m.  '))

                if ans == 'q':
                    f.close()    
                    vidx.close()
                    print 'You stopped at index {0}'.format(i)
                    return assess

                if i == start:
                    assess.append(ans)
                else:                
                    assess[i - start].append(ans)

                # immediately record vet and flare indices
                f.write('{:^3}'.format(ans))
                # ynm?
                vidx.write(' {0:>4} {1:>02} '.format(flags[i][0],
                flags[i][end] - flags[i][0] + 1))

            # for multiple events
            else:
                beg = 0
                for j in range(num_events):
                    if j == 0:
                        print '{0}, light curve {1} of {2}: {3} events found'.format(files[i], i+1, len(files), num_events)
                        
                    # isolate individual events and prompt user
                    # response
                    if j == bins:
                        end = -1 # last flag
                    else:                    
                        end = events[j]

                    xmargin = 6 * (days[flags[i][end]] -
                                   days[flags[i][beg]])
                    ymargin = (np.max(normflux[flags[i][beg:end]]) -
                               np.min(normflux[flags[i][beg:end]])) # diff between hi and lo

                    plt.title('Event {0}'.format(j + 1))
                    plt.xlim(days[flags[i][beg]] - xmargin,
                             days[flags[i][end]] + xmargin)
                    # avg = 0.5 * (np.max(normflux[flags[i][beg:end]])
                    #              +
                    #              np.min(normflux[flags[i][beg:end]])) # hi/lo avg = center
                    plt.ylim(np.min(normflux[flags[i][beg:end]]) -
                             2*ymargin, np.max(normflux[flags[i][beg:end]]) +
                             ymargin)

                    # limit the keystrokes to y/n/m and take out the
                    # return here ideally
                    ans = str(raw_input('Does this look like a flare? Say y, n, m.  '))

                    # add response to pseudo-array
                    if ans == 'q':
                        f.close()    
                        vidx.close()
                        print 'You stopped at index {0}'.format(i)
                        return assess

                    if i == start:
                        assess.append(ans)
                    else:                
                        assess[start - i].append(ans)

                    # immediately write response and flares to output
                    f.write('{:^3}'.format(ans))
                    # add one to include endpoint (we want numpoints)
                    vidx.write(' {0:>4} {1:>02}'.format(flags[i][beg],
                                                        flags[i][end] - flags[i][beg] + 1))

                    # set up the beginning of the next event
                    if j == bins:
                        pass
                    else:
                        beg = events[j] + 1

            # create sublist for next light curve except last
            if i == start:
                assess = [assess, []]
            elif i == len(files) - 1:
                pass
            else:
                assess.append([])

        if flags is None:
            print 'Light curve {0} of {1}'.format(i+1, len(files))
        if flags is not None and vet is False:
            print 'Light curve {0} of {1}: {2} events found.'.format(i+1, len(files), num_events)
        x = str(raw_input('Hit return for the next light curve. \n'))
        if x == "q":
            f.close()    
            vidx.close()
            print 'You stopped at index {0}'.format(i)
            return assess
        if x == "":
            plt.close()

    if vet:
        f.close()
        vidx.close()
        return assess
    else:
        return

def smooth(x, window_len=11, window='flat'):
    """Smooth the data x using a window with requested size.
    
    This method is based on the convolution of a scaled window with
    the signal.  The signal is prepared by introducing reflected
    copies of the signal (with the window size) in both ends so that
    transient parts are minimized in the begining and end part of the
    output signal.
    
    input:
        x: the input signal 

        window_len: the dimension of the smoothing window; should be
        an odd integer

        window: the type of window from 'flat', 'hanning', 'hamming',
        'bartlett', 'blackman'. Flat window will produce a moving
        average smoothing.

    output:
        y: the smoothed signal
        
    example:
        t = linspace(-2, 2, 0.1)
        x = sin(t) + randn(len(t)) * 0.1
        y = smooth(x)
    
    See also: 

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array
    instead of a string
    NOTE: length(output) != length(input), to correct this: return
    y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth() only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is not of 'flat', 'hanning', "
                         + "'hamming', 'bartlett', 'blackman'")
    # not really much to do here...
    if window_len < 3:
        return x

    # concatenate beginning and end stubs for wrap-around 
    beg = x[(window_len - 1) : 0 : -1]
    end = x[-1 : -(window_len) : -1]
    concat = np.r_[beg, x, end]

    #print(len(concat))
    if window == 'flat':   #moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w/w.sum(), concat, mode='valid')

    # match up length of output with input
    if window_len % 2 != 0:
        return y[(window_len/2) : -(window_len/2)]
    else:
        return y[(window_len/2 - 1) : -(window_len/2)]


def intFlare(filenames, flaglist, stride=20, window='flat',
    plot=False, save=False, output="intstats.out"):
    """Isolates flare events using smoothing and integrates under
    flare. Returns list of stats from each event.

    Inputs:
    ------
    filenames: file containing a list of files that contain the brightness data

    flaglist: list of dictionaries with keys 'id' and 'flags.' The
    value of 'flags' is a flattened list of indices

    Output:
    ------
    allstats: list of dictionaries with keys 'id' and 'stats.'  The
    value of 'stats' is a num_events x 3 array.
        Column 1: integrated flux under the event
        Column 2: maximum flux value of event
        Column 3: duration of the event (in hours)
    """

    with open(filenames, 'r') as f:
        files = [line.strip() for line in f]
    kids = getNumeric(files)  # list of kids
    allstats = list()             # output list of dictionaries
    for d in flaglist:
        if d['id'] in kids:
            filename = files[kids.index(d['id'])]
            t, normflux = ltcurve(filename, display=False)
            # omit the flagged points and interpolate before smoothing
            if d["flags"]:
                # remove points on either side for better interpolation
                flags = d["flags"]
                events = np.flatnonzero(np.diff(flags) - 1)
                extFlags = list(flags)  # defensive copy with extra flagged points
                offset = 0              # how many flags have been inserted?
                for i in events:
                    extFlags.insert(i+1 + offset, flags[i+1] - 1)  # shave before next
                    extFlags.insert(i+1 + offset, flags[i] + 1)    # shave after current
                    offset += 2
                # corner cases
                if flags[0] != 0:
                    extFlags.insert(0, flags[0] - 1)
                if flags[-1] < len(t) - 1:
                    extFlags.insert(len(extFlags), flags[-1] + 1)
                cutT = np.delete(t, extFlags)
                cutFlux = np.delete(normflux, extFlags)            
                fflux = sp.interp1d(cutT, cutFlux, kind='nearest') # flux function
                fsmoothed = smooth(fflux(t), stride, window)
            else: # no flags argument
                fsmoothed = smooth(normflux, stride, window)

            if fsmoothed.size != normflux.size:
                raise ValueError("Smoothed and original data not of same length.")

            resids = normflux - fsmoothed # noise and flare events
            if plot:
                fig = plt.figure()
                fig.subplots_adjust(left=.06, bottom=.07, right=.98,
                                    top=.92, wspace=.1)

                plt.subplot(121)
                plt.plot(t, normflux, label='normalized')
                plt.plot(t[flags], normflux[flags], 'y*')
                plt.plot(t, fsmoothed, label='{}, stride {}'.format(window,stride))
                plt.title(filename.partition('.')[0])    # strip extension
                plt.legend()
            
                # plot residuals
                plt.subplot(122)
                plt.plot(t, resids, label='{}, stride {}'.format(window, stride))
                plt.plot(t[flags], resids[flags], 'y*')
                plt.title(filename.partition('.')[0] + ' events')
                plt.legend()

            if d["flags"]:
                # add 4 or 5 points on either side of event to be integrated over
                intBounds = np.flatnonzero(np.diff(extFlags) - 1) # bounds of the events

                stats = list()      # 2D array for each event
                for i in range(intBounds.size + 1):
                    if intBounds.size == 0:
                        beg = 0
                        end = -1
                    elif i == 0:
                        beg = 0
                        end = intBounds[i]
                    elif i == intBounds.size:
                        beg = intBounds[i-1] + 1
                        end = -1
                    else:
                        beg = intBounds[i-1] + 1
                        end = intBounds[i]
                        
                    [loBound, hiBound] = utils.window(beg, end, 4, len(t))
                    if ((loBound < 0) or (hiBound >= resids.size)):
                        print d['id']
                        raise ValueError("Invalid bounds of integration")

                    # duration (for now, just number of points * 0.5 hrs)
                    # account for +1 flag on either side of integration window
                    duration = 0.5 * (extFlags[end] - extFlags[beg] - 1)
                    sums = si.trapz(resids[loBound:hiBound], t[loBound:hiBound])
                    max_flux = np.max(resids[loBound:hiBound])
                    stats.append([sums, max_flux, duration])

                allstats.append({'id': d['id'], 'stats':
                                     np.array(stats)})
    if save:
        np.savez(output, int_stats=allstats)
    return allstats


def getEvents(file1, file2, option='y'):
    """Generates possible event flags given the names of vet output
    files file1 and file2, where
        file1: Kepler ids and y/n/m status
        file2: Kepler ids, first index of event, duration
        option: 'y' or 'm' depending on whether you want to integrate
                 the yesses or maybes
    """

    # open the vet files for reading
    f = open(file1, 'r')
    fidx = open(file2, 'r')

    # readline from file1, split by whitespace
    vetTokens = f.readline().split()
    lineLength = len(vetTokens)

    flags = list()
    # loop through ynm file
    while lineLength > 0:
        kid = vetTokens[0]
        # unfinished ynm file
        if lineLength == 1:
            print 'Reached line without vets\n'
            f.close()
            fidx.close()
            return flags
        # find corresponding kid in info file2
        infoTokens = None
        for line in fidx:
            if line.find(kid) != -1:
                infoTokens = line.split()
                fidx.seek(0)    # next search to start from beginning
                break
        # check to make sure valid line
        if infoTokens == None:
            raise ValueError('Could not find kid {} in {}'.format(kid,file2))
        if kid != infoTokens[0]:
            raise ValueError("Something went wrong with identifying the correct id")

        # traverse token arrays
        kidict = {'id': kid, 'flags': list()}
        for i in range(1, len(vetTokens)):
            # if y or m, generate flags:
            if vetTokens[i] == 'n':
                continue
            elif vetTokens[i] == option:
                # for h in range(len(infoTokens)):
                #     print infoTokens[h]
                idx = (i * 2) - 1            # corresponding index into info file2
                start = int(infoTokens[idx]) # beginning of run
                duration = int(infoTokens[idx + 1])
                for t in range(duration): 
                    kidict["flags"].append(start + t) # append to last nested list

        # update for next loop iteration
        flags.append(kidict)
        vetTokens = f.readline().split()
        lineLength = len(vetTokens)

    return flags


# include formatting for each file (the outputs of our functions)

# note to self: #13 is amazing, kid11546211
# super bonus: no return on ynm, or when switching to new curve
