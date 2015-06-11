#!/usr/bin/env python

# Process a cigar string.

# Note that this module includes a set of unit tests, which can be
# executed by running this file directly from the command line.

import os
import sys
import re                       # for regular expressions
import math
from tt_log import logger

VERSION = '20150611.01'
logger.debug('version %s loaded' % VERSION)

DEF_Q = 50.0                    # Q score for error-free exon

regexFields   = re.compile('(\d+)([MNDISH])')
regexLeading  = re.compile('^(\d+)S')
regexTrailing = re.compile('(\d+)S$')
regexDigits   = re.compile('(\d+)')
regexBases    = re.compile('([ACGTNacgtn]+)')
regexMismatch = re.compile('([0ACGTNacgtn]+)')

class Exon (object):

    def __init__ (self, start, end, offset):
        self.start   = start     # genomic reference start coord
        self.end     = end       # genomic reference end coord
        self.offset  = offset    # offset of start of exon into read
        self.size    = 0         # size of exon in read (see note in 'exons' method')
        self.inserts = 0         # initial values, need to be updated
        self.deletes = 0
        self.substs  = None      # subst error count available only if MD supplied

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

        if self.substs is None:
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

    def __init__ (self, string, start, MD=None):

        self.string  = string
        self.start   = start     # genomic start coordinate

        match = re.search(regexLeading, self.string)
        self.leading = int(match.group(1)) if match is not None else 0
        match = re.search(regexTrailing, self.string)
        self.trailing = int(match.group(1)) if match is not None else 0

        self.MD      = MD        # MD field from SAM file: may not be present
        self.cfields = re.findall(regexFields, string)
        if MD is not None:
            self.expandCigarString()           # add substitution errors to cigar string
        else:
            self.expCfields = self.cfields     # use the original cigar, no subst error counts available

        self.exonList = None     # created and cached when/if needed
        self.variants = None

    def prettyPrint (self):
        '''Separate cigar fields into exons and introns, separated by spaces.'''

        fields  = list()
        filling = False

        for count, op in self.cfields:          # count is a string

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

    def addExon (self, start, offset):
        '''Append a new exon to self.exonList.'''

        if self.exonList == None:
            self.exonList = list()

        self.exonList.append (Exon(start, start-1, offset))

        if self.MD is not None:
            self.exonList[-1].substs = 0             # initial value, needs to be computed later

        return

    def exons (self):
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

            self.variants = list()

            offset = self.leading
            self.addExon (self.start, offset)          # first exon

            exonList = self.exonList                   # make it local to avoid some lookups

            for cf in self.expCfields:

                count = int(cf[0])
                op    = cf[1]                          # cf[2] (refBases) may or may not be present

                if op == 'M':
                    exonList[-1].end  += count
                    exonList[-1].size += count
                elif op == 'D':
                    self.variants.append ( ['D', count, exonList[-1].end+1, offset+exonList[-1].size, cf[2]] )     # +1: end is inclusive
                    exonList[-1].end += count
                    exonList[-1].deletes += count      # bump count of deletions
                elif op == 'I':
                    self.variants.append ( ['I', count, exonList[-1].end+1, offset+exonList[-1].size] )     # no refBases for an insert
                    exonList[-1].size += count         # inserted bases are part of the read length
                    exonList[-1].inserts += count      # bump count of insertions
                elif op == 'X':
                    self.variants.append ( ['X', count, exonList[-1].end+1, offset+exonList[-1].size, cf[2]] )
                    if exonList[-1].substs is None:    # None is initial value, assuming no MD string
                        exonList[-1].substs = 0        # if there's an 'X', we've seen an MD string
                    exonList[-1].substs += count
                    exonList[-1].end    += count
                    exonList[-1].size   += count
                elif op == 'N':
                    # start a new exon -- zero-length until we update it
                    offset += exonList[-1].size
                    self.addExon (exonList[-1].end+count+1, offset)

            # This shouldn't happen. But if we generated a zero-length exon, get rid of it.

            if exonList[-1].start > exonList[-1].end:
                exonList.pop()      # perhaps this should be a fatal error

            # Adjust for soft-clips.

####            exonList[0].start -= self.leading
####            exonList[-1].end  += self.trailing

        return self.exonList

    def expandCigarString (self):
        '''
        Apply the MD string to the cigar string, expanding M
        (match/mismatch) fields into M (match) and X (mismatch) fields.
        '''

        MD = self.MD
        if MD is None:                 # no MD string supplied
            return

        self.expCfields = list()
        count = 0                      # count=0 will cause cfIx to be bumped
        cfIx   = -1                    # index into cfields list (bumped prior to first use)
        strIx  = 0                     # index into MD string
        pendingMatches = 0

        while strIx < len(MD):

            # At this point, we're expecting a field length. But in
            # the case of two adjacent substitutions, some aligners
            # cheat: rather than 123A0C0T, they produce 123ACT. We
            # cope with that by doing nothing. I've added a unit test
            # case for this.

            match = regexDigits.match(MD, strIx)
            if match is not None:
                fieldLen = match.group(1)
                pendingMatches = int(fieldLen)
                strIx  += len(fieldLen)
            else:
                pendingMatches = 0

            while True:               # loop through cfields entries

                if count == 0:
                    cfIx += 1
                    if cfIx >= len(self.cfields):
                        break
                    count = int(self.cfields[cfIx][0])
                    op = self.cfields[cfIx][1]

                if op in 'SIN':                                      # pass straight through
                    self.expCfields.append(self.cfields[cfIx])
                    count = 0                                        # cause a bump

                elif op == 'M' and count <= pendingMatches:          # there's a bigger match (possibly spanning introns) in MD string
                    self.expCfields.append( [count, 'M'] )           # count may have been modified
                    pendingMatches -= count
                    count = 0                                        # cause a bump

                else:
                    break

            # At this point, either:
            #     op is D, MD[] is a deletion, pending=0, or
            #     op is M, MD[] is a substitution, pending >= 0, or
            #     we're at the end of the MD string

            if strIx >= len(MD):
                break                                               # end of MD string

            if op == 'D':

                if pendingMatches > 0:
                    raise RuntimeError ('unexpected pending matches at %d in MD string %s' % (strIx, MD))

                if MD[strIx] != '^':              # A deletion (e.g., ^AGT)
                    raise RuntimeError ('expected ^ at %d in MD string %s' % (strIx, MD))

                strIx += 1                        # bump past the ^
                match = regexBases.match(MD, strIx)
                if match is None:
                    raise RuntimeError ('no bases after ^ in MD string %s at %d' % (MD, strIx))

                refBases = match.group(1)
                fieldLen = len(refBases)
                if fieldLen != count:
                    raise RuntimeError ('mismatching field length at %d in MD string %s' % (strIx, MD))
                strIx  += fieldLen

                self.expCfields.append( [self.cfields[cfIx][0], self.cfields[cfIx][1], refBases] )
                count = 0                                                     # cause a bump

            elif op == 'M':                                                   # an M field with a mismatch in it somewhere

                if pendingMatches > 0:                                        # any remaining matches?
                    self.expCfields.append( [str(pendingMatches), 'M']  )     # new count needs to be a string
                    count -= pendingMatches
                    pendingMatches = 0

                match = regexMismatch.match(MD, strIx)                        # this will find 'T0T0G0A'
                if match is None:
                    raise RuntimeError ('expecting substitution at %d in MD string %s' % (strIx, MD))
                refBases = match.group(1)

                strIx  += len(refBases)                                       # index increment includes the zeroes
                refBases = refBases.replace ('0', '')                         # now get rid of the zeroes
                self.expCfields.append( [str(len(refBases)), 'X', refBases] ) # count of substs does not include zeroes
                count -= len(refBases)

            else:
                raise RuntimeError ('unexpected operation %s in cigar string field %d' % (op, cfIx))

        return

    def eCigAsString (self):
        '''produce a printable string from the expanded cigar fields.'''

        return ' '.join([ '%s%s' % (x[0], x[1]) for x in self.expCfields ])

    def printVariantList (self, bases):

        for var in self.variants:
            type, count, start, offset = var[0:4]
            refBases = var[4] if len(var) == 5 else '-'      # reference base is present unless it's an insert error
            readBases = bases[offset:offset+len(refBases)] if type != 'D' else '-'
            print 'var:      %s  %2d  %9d  %5d  %s  %s' % (type, count, start, offset, refBases, readBases)

        return


def unitTest ():

    testcases = ( ['100M', '50T49', [1], '50M 1X 49M'],
                  ['100M 100N 100M', '99T100', [1,0], '99M 1X 100N 100M'],
                  ['100M 100N 100M 3S', '200', [0,0], '100M 100N 100M 3S'],
                  ['10M4D20M 100N 10M', '10^ACGT10T9c9', [1,1], '10M 4D 10M 1X 9M 100N 1X 9M'],
                  ['20M 100N 20M 100N 20M 100N 20M', '10t49C4C4C4c0c0c2', [1,0,0,6], '10M 1X 9M 100N 20M 100N 20M 100N 1X 4M 1X 4M 1X 4M 1X 1X 1X 2M'],
                  ['3S 100M 100N 80M3D20M', '50ACT97TCA27^TTT0T19', [3,4], '3S 50M 1X 1X 1X 47M 100N 50M 1X 1X 1X 27M 3D 1X 19M'],
                  ['100M10D5M 100N 5M 100N 15M', '100^AAAAATTTTT7ACT15', [0,3,0], '100M 10D 5M 100N 2M 1X 1X 1X 100N 15M'],
                  ['533M1I176M1I66M1I475M1I525M 6725N 62M 10126N 58M 4S', '1837G0G0C0T0C0A1C0T0G0C1T0C0T0G0A0G0C1G0C36', [0,0,19],
                   '533M 1I 176M 1I 66M 1I 475M 1I 525M 6725N 62M 10126N 1X 1X 1X 1X 1X 1X 1M 1X 1X 1X 1X 1M 1X 1X 1X 1X 1X 1X 1X 1M 1X 1X 36M 4S'],
                  )

    start = 100000

    for ix, test in enumerate(testcases):

        cigar = test[0].replace(' ', '',)        # blanks are for readability, get rid of them
        cigObj = CigarString(cigar, start, MD=test[1])
        exonList = cigObj.exons()
        subs = [exon.substs for exon in exonList]

        if subs == test[2]:
            print 'test %d error counts passed' % (ix+1)
        else:
            print 'test %d error counts failed:' % (ix+1),
            print subs

        eCig = cigObj.eCigAsString()

        if eCig == test[3]:
            print 'test %d expanded cigar passed' % (ix+1)
        else:
            print 'test %d expanded cigar failed:' % (ix+1),
            print eCig

    return


if __name__ == "__main__":
    unitTest()

