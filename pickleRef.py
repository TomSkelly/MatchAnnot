#!/usr/bin/env python

# Read a reference fasta file, write a Reference object to a pickle
# file for use as input to findOverlaps.py, etc.

import os
import sys
import optparse
import re        # regular expressions
import cPickle as pickle

from tt_log import logger
import Reference as ref

VERSION = '20150105.01'

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    # Reading a pickled file and repickling it doesn't make sense. But one day it will...

    if opt.format == 'pickle':
        refObj = ref.Reference.fromPickle (opt.ref)
    else:
        refObj = ref.Reference (opt.ref)

    refObj.toPickle(opt.output)

    logger.debug('finished')

    return


def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options]', version=VERSION)

    parser.add_option ('--ref',       help='reference file, in format specified by --format')
    parser.add_option ('--format',    help='format of reference file (def: %default)', \
                           type='choice', choices=['standard', 'fasta', 'pickle'])      # 'standard' is the same as 'fasta'
    parser.add_option ('--output', help='output file name (required)')

    parser.set_defaults (format='fasta',
                         )

    opt, args = parser.parse_args()

    return opt, args


if __name__ == "__main__":
    main()
