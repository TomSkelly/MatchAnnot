#!/usr/bin/env python

# Process a cigar string.

import os
import sys
import re                       # for regular expressions
from tt_log import logger

regexFields   = re.compile('(\d+)([MNDISH])')
regexLeading  = re.compile('^(\d+)S')
regexTrailing = re.compile('(\d+)S$')

class Exon (object):

    def __init__ (self, start, end):
        self.start   = start
        self.end     = end
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


class CigarString (object):

    def __init__ (self, string):

        self.string = string
        self.cfields = re.findall(regexFields, string)

        match = re.search(regexLeading, self.string)
        self.leading = int(match.group(1)) if match is not None else 0

        match = re.search(regexTrailing, self.string)
        self.trailing = int(match.group(1)) if match is not None else 0

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
        each exon.

        The interval coordinates reflect the reference, so
        '100M10I200M' is an exon of length 300 (inserted bases are not
        present in reference), and '100M10D200M' is an exon of length
        310 (there are 10 bases in the reference that don't appear in
        the string described by the cigar).
        '''

        exonList = [ Exon(start, start-1) ]    # first exon

        for count, op in self.cfields:

            count = int(count)

            if op == 'M':
                exonList[-1].end += count
            elif op == 'D':
                exonList[-1].end += count
                exonList[-1].deletes += count      # bump count of deletions
            elif op == 'I':
                exonList[-1].inserts += count      # bump count of insertions
            elif op == 'N':
                # start a new exon -- zero-length until we update it
                exonList.append ( Exon(exonList[-1].end + count + 1, exonList[-1].end + count) )

        # This shouldn't happen. But if we generated a zero-length exon, get rid of it.

        if exonList[-1].start > exonList[-1].end:
            exonList.pop()      # perhaps this should be a fatal error

        # Adjust for soft-clips.

        exonList[0].start -= self.leading
        exonList[-1].end  += self.trailing

        return exonList
