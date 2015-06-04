#!/usr/bin/env python

# Read a MatchAnnot output file, output entries for selected gene9s)
# only. Optionally reverse the order of exon match lines.

import os
import sys
import optparse
import re

from tt_log import logger

VERSION = '20150331.01'

REGEX_EXONS = re.compile ('(\s*\S+)')

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    geneList = list()
    if opt.genes is not None:
        geneList = opt.genes.split(',')
    gene = None

    if len(args) > 0:
        logger.debug('reading matchAnnot file %s' % args[0])
        handle = open (args[0], 'r')
    else:
        logger.debug('reading matchAnnot data from stdin')
        handle = sys.stdin
    
    exonLines  = list()
    entryLines = list()

    for line in handle:

        if line.startswith ('exon:'):
            exonLines.append(line)                   # save a batch of exon lines
        else:

            if len(exonLines) > 0:                   # if there is a batch pending
                if opt.flip and strand == '-':       # reverse order if requested
                    exonLines = reverseExonList (exonLines)
                entryLines.extend(exonLines)         # add them to the output list
                exonLines = list()

            entryLines.append (line)                 # add non-exon line directly to output list

            if line.startswith ('isoform:'):
                strand = line.split()[-2]

            elif line.startswith ('gene:'):
                gene = line.split()[1]

            elif line.isspace():                     # last line of an entry?
                if len(entryLines) > 0:
                    if opt.genes is None or gene in geneList or entryLines[0].startswith('summary:'):
                        sys.stdout.writelines (entryLines)
                    entryLines =  list()

    if len(entryLines) > 0:
        if entryLines[-1].startswith('summary:'):
            sys.stdout.writelines (entryLines)

    handle.close()

    logger.debug('finished')

    return

def reverseExonList (exonList):
    '''
    Given a list of exon lists, reverse the order of the lines, and
    reverse the begin and end columns in each line.
    '''

    newList =  list()

    for line in exonList:

        # The fields we want to swap are right-justified with a
        # variable number of leading spaces. Split the line into
        # fields including the leading spaces, swap, and rejoin.

        fields = re.split (REGEX_EXONS, line)
        fields = [f for f in fields if f != '']    # split throws in a few null strings: get rid of them
        fields[3:6], fields[6:9] = fields[6:9], fields[3:6]

        newList.append (''.join(fields))

    newList.reverse()

    return newList        

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <SAM_file> ... ', version=VERSION)

    parser.add_option ('--genes',     help='comma-separated list of genes to select (def: all)')
    parser.add_option ('--flip',      help='print exon matches in transcript order', action='store_true')

    opt, args = parser.parse_args()

    return opt, args

if __name__ == "__main__":
    main()
