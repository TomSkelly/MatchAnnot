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
import Reference   as ref

VERSION = '20150527.01'

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    # We may want to re-pickle a pickled annotation file. e.g., to add polyA annotations.

    if opt.format == 'pickle':
        annotList   = anno.AnnotationList.fromPickle (opt.gtf)
    elif opt.format == 'alt':
        annotList   = anno.AnnotationList (opt.gtf, altFormat=True)
    else:
        annotList   = anno.AnnotationList (opt.gtf)

    if opt.ref is not None:            # if a reference was specified, look for PolyA tracts in exons
        refObj = ref.Reference.fromPickle (opt.ref)
        annotList.annotatePolyA (refObj)

    annotList.toPickle(opt.output)

    logger.debug('finished')

    return


def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options]', version=VERSION)

    parser.add_option ('--gtf',    help='annotations file, in format specified by --format')
    parser.add_option ('--format', help='annotations in alternate gtf format (def: %default)', \
                           type='choice', choices=['standard', 'alt', 'pickle'])
    parser.add_option ('--ref',    help='reference file to search for genomic polyAs, in pickle format (def: none)')
    parser.add_option ('--output', help='output file name (required)')

    parser.set_defaults (format='standard',
                         )

    opt, args = parser.parse_args()

    return opt, args


if __name__ == "__main__":
    main()
