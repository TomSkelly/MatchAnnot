#!/usr/bin/env python

# Process a cigar string.

# Note that this module includes a set of unit tests, which can be
# executed by running this file directly from the command line.

import os
import sys
import re                       # for regular expressions
import math
from tt_log import logger

VERSION = '20150121.01'
logger.debug('version %s loaded' % VERSION)

DEF_Q = 50.0                    # Q score for error-free exon

regexFields   = re.compile('(\d+)([MNDISH])')
regexLeading  = re.compile('^(\d+)S')
regexTrailing = re.compile('(\d+)S$')
regexDigits   = re.compile('(\d+)')
regexBases    = re.compile('([ACGTNacgtn]+)')

class Exon (object):

    def __init__ (self, start, end):
        self.start   = start     # genomic reference start coord
        self.end     = end       # genomic reference end coord
        self.size    = 0         # size of exon in read (see note in 'exons' method')
        self.inserts = 0         # initial values, need to be updated
        self.deletes = 0

    def __lt__ (self, other):
        '''Support sorting of exons based on position then size.'''

        if self.start < other.start:
            return True
        elif self.start > other.start:
            return False
        else:
            return self.end < other.end

    def QScore (self):
        '''Compute Q score for THIS exon.'''

        if not hasattr (self, 'substs'):
            Q = 0.0            # no MD string was supplied, can't compute Q score
        else:

            errors = self.substs + self.inserts + self.deletes

            if (errors == 0):
                Q = DEF_Q
            else:
                length = self.end - self.start + 1
                Q = -10.0 * math.log10( float(errors) / float(length) )

        return Q

class CigarString (object):

    def __init__ (self, string, MD=None):

        self.string  = string
        self.MD      = MD        # MD field from SAM file: may not be present
        self.cfields = re.findall(regexFields, string)

        match = re.search(regexLeading, self.string)
        self.leading = int(match.group(1)) if match is not None else 0

        match = re.search(regexTrailing, self.string)
        self.trailing = int(match.group(1)) if match is not None else 0

        self.exonList = None     # created and cached when/if needed

    def prettyPrint (self):
        '''Separate cigar fields into exons and introns, separated by spaces.'''

        fields  = list()
        filling = False

        for count, op in self.cfields:

            if filling:                         # if we're currently concatenating
                if op in 'SN':                  # if this should be a new field
                    fields.append(count+op)
                    filling = False             # S & N should stand alone
                else:
                    fields[-1] += count+op      # add to existing field
            else:
                fields.append(count+op)         # start a new field
                if op not in 'SN':              # and add subsequent matches to it
                    filling = True

        return ' '.join(fields)

    def genomicLength (self):
        '''Compute the number of bases of the total genomic reference covered by this read.'''

        length = 0

        for count, op in self.cfields:

            # M = bases match (or mismatch) one-for-one: count them
            # N = bases in reference but not in read (e.g., introns in transcript read): count them
            # D = bases in reference but not in read: count them
            # I = bases in read but not in reference: do not count
            # S = bases in read prior to start of alignment: do not count
            # H = bases in read prior to start of alignment: do not count

            if op == 'M' or op == 'N' or op == 'D':
                length += int(count)
                
        return length


    def softclips (self):
        return self.leading, self.trailing

    def fields (self):
        '''Generator function to return (count, op) for each field in cigar string.'''

        for count, op in self.cfields:
            yield int(count), op

        return

    def exons (self, start):
        '''
        Given a cigar string, return a list of Exon objects
        representing exons in the alignment. N's in the cigar string
        separate the exons. The returned objects also include counts
        of insertions and deletions, relative to the reference, in
        each exon. Substitution error count is also present if an MD
        string was supplied.

        The interval start and end coordinates reflect the reference,
        so '100M10I200M' is an exon of length 300 (inserted bases are
        not present in reference), and '100M10D200M' is an exon of
        length 310 (there are 10 bases in the reference that don't
        appear in the string described by the cigar).

        OTOH, self.size reflects the number of bases in the exon as
        found in the actual sequence read. So the '100M10I200M' exon
        mentioned above will have a size of 310.
        '''

        if self.exonList == None:

            exonList = [ Exon(start, start-1) ]    # first exon
            self.exonList = exonList               # cache the eventual list now

            for count, op in self.cfields:

                count = int(count)

                if op == 'M':
                    exonList[-1].end  += count
                    exonList[-1].size += count
                elif op == 'D':
                    exonList[-1].end += count
                    exonList[-1].deletes += count      # bump count of deletions
                elif op == 'I':
                    exonList[-1].size += count         # inserted bases are part of the read length
                    exonList[-1].inserts += count      # bump count of insertions
                elif op == 'N':
                    # start a new exon -- zero-length until we update it
                    exonList.append ( Exon(exonList[-1].end + count + 1, exonList[-1].end + count) )

            # This shouldn't happen. But if we generated a zero-length exon, get rid of it.

            if exonList[-1].start > exonList[-1].end:
                exonList.pop()      # perhaps this should be a fatal error

            self.applyMDString()

            # Adjust for soft-clips.

            exonList[0].start -= self.leading
            exonList[-1].end  += self.trailing

        return self.exonList

    def applyMDString (self):
        '''Use MD string, if present, to fill in substitution error counts by exon.'''

        # Internal method, should only be called from 'exons' method above.

        MD = self.MD
        if MD is None:                 # no MD string supplied
            return

        for exon in self.exonList:
            exon.substs = 0            # creates a new attribute

        exonIx = 0                     # index into list of exons
        strIx  = 0                     # index into MD string
        readIx = 0                     # index into the read
        cursor = self.exonList[0].end - self.exonList[0].start + 1   # # of bases described by MD string so far

        while strIx < len(MD) and not MD[strIx:].isdigit():          # step through the MD string, ...
                                                                     #    quit when only matches remain

            # At this point, we're expecting a field length. But in
            # the case of two adjacent substitutions, some aligners
            # cheat: rather than 123A0C0T, they produce 123ACT. We
            # cope with that by doing nothing. I've added a unit test
            # case for this.

            match = regexDigits.match(MD, strIx)
            if match is not None:
                fieldLen = match.group(1)
                readIx += int(fieldLen)
                strIx += len(fieldLen)

            while cursor <= readIx:             # bump to exon containing the error
                exonIx += 1
                if exonIx >= len(self.exonList):
                    raise RuntimeError ('MD string extends beyond last exon %s' % MD)
                cursor += self.exonList[exonIx].end - self.exonList[exonIx].start + 1

####            print 'strIx: %d  exonIx: %d  readIx: %d  cursor: %d' % (strIx, exonIx, readIx, cursor)

            if MD[strIx] in 'ACGTNacgtn':       # substitution error?
                self.exonList[exonIx].substs += 1
                readIx += 1
                strIx  += 1

            elif MD[strIx] == '^':              # A deletion (e.g., ^AGT)

                strIx += 1
                match = regexBases.match(MD, strIx)
                if match is None:
                    raise RuntimeError ('no bases after ^ in  MD string %s at %d' % (MD, strIx))

                fieldLen = len(match.group(1))
                strIx  += fieldLen
                readIx += fieldLen

            else:
                raise RuntimeError('unexpected field in MD string %s at 5d' % (MD, strIx))

        return

def unitTest ():

    testcases = ( ['100M100N100M', '200', [0,0] ],
                  ['10M4D20M100N10M', '10^ACGT10T9c9', [1,1] ],
                  ['20M100N20M100N20M100N20M', '10t49C4C4C4c0c0c2', [1,0,0,6] ],
                  ['100M100N80M3D20M', '50ACT97TCA27^TTT0T19', [3,4] ],
                  )

    start = 100000

    for ix, test in enumerate(testcases):

        exonList = CigarString(test[0], MD=test[1]).exons(start)
        subs = [exon.substs for exon in exonList]

        if subs == test[2]:
            print 'test %d passed' % (ix+1)
        else:
            print 'test %d failed: ' % (ix+1),
            print subs

    return


if __name__ == "__main__":
    unitTest()

