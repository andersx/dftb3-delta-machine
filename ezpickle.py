#!/usr/bin/env python2

"""
Short wrapper functions to easily handle opening and 
saving Python objects to files in cPickle formats.

"""
import cPickle

def load(filename):
    """ Returns Python data from a cPickle formatted file.

    Arguments:
    filename -- name of the file to be parsed.

    """
    with open(filename, "rb") as f:
        p = cPickle.load(f)

    return p

def save(data, filename):
    """ Saves about any Python data into a cPickle formatted file.

    Arguments:
    data -- the data or object to be saved.
    filename -- name of the file to be parsed.

    NOTE: Will overwrites any pre-existing file with the
    name filename.

    """

    with open(filename, "wb") as f:
        cPickle.dump(data, f, protocol=2)

