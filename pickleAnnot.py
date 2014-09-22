#!/usr/bin/env python

# Read annotation file, write an AnnotationList object to a pickle
# file for use as input to matchAnnot.py.

import os
import sys
import optparse
import re        # regular expressions
import cPickle as pickle

from tt_log import logger
import Annotations as anno

def main ():

    opt, args = getParms()

    annotList = anno.AnnotationList (opt.gtf)

    handle = open (opt.output, 'w')
    pk = pickle.Pickler (handle, pickle.HIGHEST_PROTOCOL)
    pk.dump (annotList)
    handle.close()

    logger.debug('finished')

    return


def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <fasta_file> ... ')

    parser.add_option ('--gtf',    help='annotations in gtf format (required)')
    parser.add_option ('--output', help='output file name (required)')

    parser.set_defaults (gtf=None,
                         output=None,
                         )

    opt, args = parser.parse_args()

    return opt, args


if __name__ == "__main__":
    main()
