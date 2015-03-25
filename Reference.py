#!/usr/bin/env python

# Interface to a fasta reference file.

import os
import sys
import re                       # for regular expressions
import cPickle as pickle

from tt_log import logger

VERSION = '20150318.01'
logger.debug('version %s loaded' % VERSION)

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

    def sequence (self, chr, start=1, howmany=None):

        if howmany is None:
            howmany = len(self.ref[chr]) - start + 1

        return self.ref[chr][start-1:start-1+howmany]      # -1: reference is stored 0-based
