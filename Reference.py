#!/usr/bin/env python

# Interface to a fasta reference file.

import os
import sys
import re                       # for regular expressions
import math                     # for ceil
import cPickle as pickle

from tt_log import logger

VERSION = '20150623.01'
logger.debug('version %s loaded' % VERSION)

DEF_WINDOW = 23
DEF_THRESH = 0.65

regexChr = re.compile ('>(\S+)')

class Reference (object):
    
    def __init__ (self, filename, altFormat=False):
        '''Create a dict, keyed by chr, from a fasta file.'''

        self.filename = filename
        self.ref      = dict()       # this is the stuff! key=chr value=sequence
        self.chrList  = list()       # preserve original order of chromosomes

        logger.debug('reading reference fasta file %s' % self.filename)

        # Reading fasta lines and appending them one-by-one is very
        # slow (painter's algorithm). Instead, accumulate individual
        # lines in a list, and join them when we've got them all. That
        # turns out to be *much* faster.

        handle = open (self.filename, 'r')

        accum = list()

        for line in handle:

            line = line.strip()

            if line.startswith('>'):

                if len(accum) > 0:
                    self.ref[chr] = ''.join(accum)     # see comment block above
                    accum = list()

                chr = re.match(regexChr, line).group(1)
                if chr in self.ref:
                    raise RuntimeError ('duplicate chromsome %s' % chr)
                self.chrList.append(chr)

            else:
                accum.append(line)

        if len(accum) > 0:              # last one?
            self.ref[chr] = ''.join(accum)

        handle.close()

    @staticmethod
    def fromPickle (filename):
        '''Create a Reference object from a pickle file (alternative to __init__).'''

        logger.debug('reading reference in pickle format from %s' % filename)

        handle = open (filename, 'r')
        pk = pickle.Unpickler (handle)
        ref = pk.load()
        handle.close()

        return ref

    def toPickle (self, filename):

        pickHandle = open (filename, 'w')
        pk = pickle.Pickler (pickHandle, pickle.HIGHEST_PROTOCOL)
        pk.dump (self)
        pickHandle.close()

        logger.debug('wrote reference data to pickle file %s' % filename)

        return

    def chromosomes (self):
        '''
        Generator function for the chromosome names in an Reference
        object. Preserves the order of the original reference file.
        '''

        for chr in self.chrList:
            yield chr

        return

    def bases (self, chr):
        '''Given a chromosome name, return the full sequence for it as a string.'''

        # We should only return a pointer here, so this should be more
        # efficient than the sequence method below, which makes a copy.

        ret = self.ref[chr]
        if ret is not self.ref[chr]:
            raise RuntimeError ('not a ref')
        return ret

    def sequence (self, chr, start=1, howmany=None):

        if howmany is None:
            howmany = len(self.ref[chr]) - start + 1

        return self.ref[chr][start-1:start-1+howmany]      # -1: reference is stored 0-based

    def findPolyAs (self, chr, start, end, strand, window=DEF_WINDOW, thresh=DEF_THRESH):
        '''Find all A/T-rich regions in a specified range.'''

        # At the caller interface, indexes describing genomic regions
        # will be 1-based and inclusive. I.e., start=1, end=100
        # describes the first 100 bases of a chromosome. But
        # internally here, indexes will be pythonic: 0-based and
        # exclusive of the last element. I.e., [0:100] describes those
        # same first 100 bases.

        polyAs = list()
        if end-start+1 < window:
            return polyAs
        
        look4 = 'A' if strand ==  '+' else 'T'
        winStart = start - 1                                              # 0-based index of first base of currently accepted tract
        winEnd   = winStart + window - 1                                  # 0-based index of last base of currently accepted tract
        stopAt   = min(end-1, len(self.ref[chr])-1)                       # 0-based index of last usable base
####        print 'stopAt: %2d  end: %2d' % (stopAt, winEnd)

        while winEnd <= stopAt:

            howmany = self.ref[chr].count(look4, winStart, winEnd+1)
            numNeeded = int(math.ceil( (winEnd-winStart+1) * thresh))     # rounds upward

####            print '---> start: %2d  end: %2d  needed: %2d  howmany: %2d' % (winStart, winEnd, numNeeded, howmany)

            if howmany < numNeeded:

                # Bump by the smallest number of bases that could get
                # us to the total required. E.g.: suppose window=20,
                # thresh=.9, numNeeded=18. If we have 15 bases in the
                # window (howmany=15), we need to bump by at least 3
                # to get to 18.

                winStart += numNeeded - howmany

            else:    # we've found a tract, now extend it

                while winEnd < stopAt:

                    winEnd += 1                                                   # bump window size

                    numNeeded = int(math.ceil( (winEnd-winStart+1) * thresh))     # new numNeeded for larger window
                    if self.ref[chr][winEnd] == look4:
                        howmany += 1                                              # new number of As for larger window
####                    print '     start: %2d  end: %2d  needed: %2d  howmany: %2d' % (winStart, winEnd, numNeeded, howmany)

                    if howmany < numNeeded:                                       # does larger window still work?
                        winEnd -= 1                                               # reject larger window
                        break

                polyAs.append( [winStart+1, winEnd+1, howmany] )                  # returned indexes are 1-based
                winStart = winEnd

            winEnd = winStart + window - 1

        return polyAs


class PhonyRef (Reference):
    '''Reference initialized from a string, to support unitTest.'''

    def __init__ (self, chr, sequence):    

        self.ref = {chr : sequence}


def unitTest ():

    #               0        1         2         3
    #               123456789012345678901234567890
    testcases = ( ['AAAAAAAAAAAAAAAAAAAA',           1, 99, [ [1, 20, 20] ] ],
                  ['AAAAAAAAAAAAAAAAAAAA',           1, 20, [ [1, 20, 20] ] ],
                  ['AAAAAAAAAAAAAAAAAAAA',           1, 19, [ ] ],
                  ['AAAAAAAAAAAAAAAAAAAT',           1, 20, [ [1, 20, 19] ] ],
                  ['TTTAAAAAAAAAAAAAAAAAAAATTTTTTT', 1, 99, [ [2, 23, 20] ] ],
                  ['TTTAAAAAAAAAAAAAAAAAAAATTTTTTT', 6, 99, [ [6, 25, 18] ] ],
                  ['AAAAAAAAAAAAAAAAAAAATTTTTAAAAAAAAAAAAAAAAAAAA',
                                                     1, 99, [ [1, 22, 20], [24, 45, 20] ] ],
                  ['GGAAATGGCCTTATAATAGTTTCCATTGCCTTGTAATTTTTTTCCATTTTTTTCTTTTTA',
                                                     1, 99, [ ] ],
                  )

    chr = 'chr5'
    
    for ix, test in enumerate(testcases):

        phony = PhonyRef(chr, test[0])
        result = phony.findPolyAs (chr, test[1], test[2], '-', window=23, thresh=0.8)

        if result == test[3]:
            print 'test %d passed' % (ix+1)
        else:
            print 'test %d failed:' % (ix+1)
            print result

if __name__ == "__main__":
    unitTest()

