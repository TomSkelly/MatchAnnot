#!/usr/bin/env python

# Look for a list of polyadenylation triggers in a sequence.

import os
import sys
import string
import re                       # for regular expressions

class PolyA (object):
    
    def __init__ (self, motifs=None):

        if motifs is not None:
            self.forward = motifs
        else:
            self.forward = ['AATAAA', 'AAATAA', 'ATAAAA', 'ATTAAA', 'ATAAAT', 'ATAAAG', \
                            'CAATAA', 'TAATAA', 'ATAAAC', 'AAAATA', 'AAAAAA', 'AAAAAG', 'AGTAAA'];

        compltab = string.maketrans ('ACGT', 'TGCA')         # for reverse-complementing reads
        self.reverse = [ x[::-1].translate(compltab) for x in self.forward ]

    def findMotifs (self, bases, strand, reach=None):
        '''
        Look for any of the list of motifs in the first 'reach' bases
        of a sequence string. Motifs are based on forward strand, so
        work backwards if sequence is on reverse strand. Only the
        3'-most occurrence of each motif is found.

        Act as a generator function to return motif+offset for each
        motif found. Offset is from 3' end.
        '''

        motifsFound = dict()
        lenBases = len(bases)

        if strand == '+':

            start = 0 if reach is None else lenBases - reach

            for motif in self.forward:
                where = bases.rfind(motif, start)                 # search from right end
                if where != -1:
                    motifsFound[motif] = len(bases) - where       # distance from 3' end

        else:

            end = lenBases if reach is None else reach

            for ix, motif in enumerate(self.reverse):
                where = bases.find(motif, 0, end+len(motif))      # search from left end
                if where != -1:
                    motifsFound[self.forward[ix]] = where         # distance from 3' end, use forward motif as key


        for motif in self.forward:
            if motif in motifsFound:
                yield motif, motifsFound[motif]

        return
            
