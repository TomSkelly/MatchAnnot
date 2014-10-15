#!/usr/bin/env python

# Read annotation file, print selected stuff in human-readable format.

# AUTHOR: Tom Skelly (thomas.skelly@fnlcr.nih.gov)

import os
import sys
import optparse
import re        # regular expressions
import cPickle as pickle

from tt_log import logger
import Annotations as anno

VERSION = '20141006.01'

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    if opt.gtfpickle is not None:
        handle = open (opt.gtfpickle, 'r')
        pk = pickle.Unpickler (handle)
        annotList = pk.load()
        handle.close()
    else:
        annotList   = anno.AnnotationList (opt.gtf)

    geneList = annotList.getGene (opt.gene)
    if geneList is None:
        print 'gene %s not found in annotations' % opt.gene
    elif len(geneList) != 1:
        print 'there are %d occurrences of gene %s in annotations' % (len(geneList), opt.gene)
    else:
        geneEnt = geneList[0]

        print 'gene:    ',
        printEnt (geneEnt)

        for transEnt in geneEnt.getChildren():
            print '\ntr:      ',
            printEnt (transEnt)

            for exonEnt in transEnt.getChildren():
                print 'exon:    ',
                printEnt (exonEnt)

    logger.debug('finished')

    return

def printEnt (ent):

    print '%-15s  %9d  %9d  %6d' % (ent.name, ent.start, ent.end, ent.end-ent.start+1)

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <fasta_file> ... ')

    parser.add_option ('--gtf',       help='annotations in gtf format')
    parser.add_option ('--gtfpickle', help='annotations in pickled gtf format')
    parser.add_option ('--gene',      help='gene to print')

    parser.set_defaults (gtf=None,
                         gtfpickle=None,
                         gene=None,
                         )

    opt, args = parser.parse_args()

    return opt, args


if __name__ == "__main__":
    main()
