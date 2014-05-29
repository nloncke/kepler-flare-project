"""Collection of commonly-used functions one might use for routine
tasks.
"""
import string

def trim(strlist, substr):
    """Given a list of strings, returns a copy with all occurences of
    substr removed from each list entry.
    """ 
    strings = list(strlist)
    for i in range(len(strings)):
        strings[i] = strings[i].replace(substr, '')
    return strings

def getNumeric(word):
    """Returns a string containing all the digits in word, strips away
    all non-numeric characters.
    """
    alltab=string.maketrans('','')
    nodigits=alltab.translate(alltab, string.digits)
    if type(word) == str:
        return word.translate(alltab, nodigits)
    if type(word) == list:
        return [w.translate(alltab, nodigits) for w in word]

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

def overlap(dictlist, vetfile):
    """ Collects the dictionaries in dictlist containing key-value
    pairs id: kid, where 'kid' is the first token per line in vetfile.
    Returns the intersection of the two sets and the relative
    complement of vetfile with respect to dictlist.
    """
    with open(vetfile, 'r') as f:
        kids = [line.split()[0] for line in f]
        
    joint = [item for item in dictlist if item["id"] in kids]
    comp  = [item for item in dictlist if item not in joint]
    return joint, comp
