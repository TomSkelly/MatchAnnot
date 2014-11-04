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

VERSION = '20141104.01'

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    # Reading a pickled file and repickling it doesn't make sense. But one day it will...

    if opt.format == 'pickle':
        annotList   = anno.AnnotationList.fromPickle (opt.gtf)
    elif opt.format == 'alt':
        annotList   = anno.AnnotationList (opt.gtf, altFormat=True)
    else:
        annotList   = anno.AnnotationList (opt.gtf)

    annotList.toPickle(opt.output)

    logger.debug('finished')

    return


def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options]', version=VERSION)

    parser.add_option ('--gtf',       help='annotations file, in format specified by --format')
    parser.add_option ('--format',    help='annotations in alternate gtf format (def: %default)', \
                           type='choice', choices=['standard', 'alt', 'pickle'])
    parser.add_option ('--output', help='output file name (required)')

    parser.set_defaults (gtf=None,
                         format='standard',
                         output=None,
                         )

    opt, args = parser.parse_args()

    return opt, args


if __name__ == "__main__":
    main()
