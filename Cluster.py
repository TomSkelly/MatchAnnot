#!/usr/bin/env python

# An IsoSeq cluster object.

import os
import sys
import re                       # for regular expressions
from tt_log import logger

class Cluster (object):

    def __init__ (self, name, flags, chr, start, strand, cigarString, bases):

        self.name   = name
        self.flags  = flags
        self.chr    = chr
        self.start  = start
        self.strand = strand
        self.cigarString = cigarString
        self.bases  = bases
        self.pctGC  = None     # computed and cached below

        
    def percentGC (self, window):
        '''Compute %GC content across the read using sliding window of specified size.'''

        # Note that the self.pctGC list will be shorter than
        # self.bases by the window size. Also note that the %GC
        # reported for position x reflects the window of bases
        # starting at x. When plotting, you may want to adjust the
        # coordinates to center the window on x.

        if self.pctGC is None:

            self.pctGC = list()

            for ix in xrange(len(self.bases)-window):

                gc = self.bases.count ('G', ix, ix+window) \
                   + self.bases.count ('C', ix, ix+window)
                self.pctGC.append (float(gc) / float(window) * 100.0)

        return self.pctGC
        
