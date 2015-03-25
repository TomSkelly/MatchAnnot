#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Smith-Waterman aligner with fixed gap penalty. NOTE: This is not
# going to be fast! It's intended for small jobs like finding
# adapters. We need a C version for the heavy lifting.

# To use, create the aligner object, call setRef and setRead, then
# call fillMatrix, which returns the alignment score. If you like
# that, call alignmentStrings to get printable strings showing the
# alignment.

# An aligner object can be reused for repeated alignments. Ref and
# Read are whatever you last specified.

import sys
import numpy as np
from tt_log import logger

class Aligner (object):

    # maxRef and maxRead can be passsed to the constructor to size the
    # S-W matrix. Otherwise, the matrix will grow in several steps as
    # longer reads and refs are encountered. It will not shrink.

    def __init__ (self, match=1, mismatch=-3, gap=-1, maxRef=1, maxRead=1):

        self._match    = match
        self._mismatch = mismatch
        self._gap      = gap
        self._ref      = None
        self._read     = None
        self._matrix   = np.zeros (shape=(maxRead+1, maxRef+1),
                                   dtype=int)        # read is X axis, ref is Y
        self._matrix[0,:] = 0                        # these values are constant
        self._matrix[:,0] = 0
        self._matrixValid = False

    def setPenalties (self, match=None, mismatch=None, gap=None):
        '''Specify new match, mismatch, and gap penalties.'''

        if match != None:
            self._match = match
        if mismatch != None:
            self._mismatch = mismatch
        if gap != None:
            self._gap = gap
        self._matrixValid = False                    # data in matrix no longer matches current penalties

    def getPenalties (self):
        '''Return the current penalty settings'''

        return (self._match, self._mismatch, self._gap)

    def setRef (self, ref):
        '''Set a new reference, passed as a character string.'''

        self._ref = ref
        self._matrixValid = False                    # data in matrix no longer matches reference

    def setRead (self, read):
        '''Set a new read, passed as a character string.'''

        self._read = read
        self._matrixValid = False

    def fillMatrix (self):
        '''Compute the Smith-Waterman score matrix, and return the best score.'''

        if self._ref  == None:
            raise RuntimeError("reference sequence has not been specified")
        if self._read == None:
            raise RuntimeError("read sequence has not been specified")

        rows = len(self._read)+1                     # read is Y axis
        cols = len(self._ref)+1                      # ref is X axis
        curRows, curCols = self._matrix.shape

        if (curRows < rows or curCols < cols):       # we're gonna need a bigger matrix

            self._matrix = np.zeros (shape=(max(curRows, rows), max(curCols, cols)),
                                     dtype=int)      # dtype=int is *not* the default!

            self._matrix[0,:] = 0                    # these values are constant
            self._matrix[:,0] = 0

            logger.debug ("SW matrix resized to %d x %d" % (self._matrix.shape))

        # Not sure the following bit is any faster, but it makes the
        # code easier to read.

        mx    = self._matrix               # this works as a pointer
        ref   = self._ref
        read  = self._read
        match = self._match
        miss  = self._mismatch
        gap   = self._gap

####        np.set_printoptions (threshold=999999, linewidth=120)

        for row in xrange(1,rows):
            for col in xrange(1,cols):

                if ref[col-1] == read[row-1]:
                    F1 = mx[row-1,col-1] + match
                else:
                    F1 = mx[row-1,col-1] + miss

                F2 = mx[row-1,col]   + gap
                F3 = mx[row,  col-1] + gap

                mx[row,col] = max(F1, F2, F3, 0)

####                print "row: %3d  col: %3d %s %s F1: %3d F2: %3d F3: %3d" \
####                % (row, col, ref[col-1], read[row-1], F1, F2, F3)

####            print mx

        self._matrixValid = True                      # data in matrix is current

        # We need to be careful here: Position [0] of a row of the
        # matrix is a dummy entry; the first base of the reference is
        # in position [1]. Here, we'll stick to saving a column index
        # into the matrix (since bestCol is an internal data
        # item). Later, in allScores and peakPosits, we'll switch to a
        # zero-based indexing scheme.

        self._lastRow = row
        self._bestCol = np.argmax(mx[row,0:cols])     # only check cols which were used
        score = mx[row, self._bestCol]

        return score

    # Once we've created the matrix, we may want to do various things
    # with it, such as creating a displayable representation of the
    # alignment. That is done below, implemented as a separate method
    # so we can avoid doing it if we don't want to (e.g., because the
    # alignment score was too low).

    def alignmentStrings (self, pos=None):
        '''Produce two strings showing the best alignment path through the current S-W matrix.'''

        if not self._matrixValid:
            raise RuntimeError("matrix has not been computed.")

        # Caller can specify a starting psoition for the backtrack, or
        # we can default to self._bestCol, which was set in
        # fillMatrix. Keep in mind that the position, if specified, is
        # that of the *last* base in the alignment. Because the
        # alignment is computed backwards via the walkback, the
        # starting base may not be well-defined.

        row = self._lastRow
        if pos == None:
            col = self._bestCol
        else:
            col = pos

        mx      = self._matrix
        ref     = self._ref
        read    = self._read
        readLow = read.lower()

        alignRef  = []
        alignRead = []

        while row > 0 and col > 0:                    # loop runs backwards, we'll reverse at the end

            m_mch = mx[row-1][col-1];                 # scores for match, insert, delete
            m_ins = mx[row]  [col-1];
            m_del = mx[row-1][col];

            if m_mch >= m_ins and m_mch >= m_del:     # if 'match' scores best ...

                alignRef.append  (ref[col-1])
                if ref[col-1] == read[row-1]:
                    alignRead.append (read[row-1])
                else:
                    alignRead.append (readLow[row-1]) # mismatch gets lower case

                row -= 1
                col -= 1

            elif m_ins >= m_del:                      # else if it's an insert in the ref

                alignRef.append  (ref[col-1])
                alignRead.append ('-')
                col -= 1

            else:                                     # else it's a delete in the read

                alignRef.append  ('-')
                alignRead.append (read[row-1])
                row -= 1

        alignRef.reverse()                            # these were created back-to-front
        alignRead.reverse()

        refString  = ''.join(alignRef)
        readString = ''.join(alignRead)

        return refString, readString

    def allScores (self):
        '''Return a tuple containing alignment scores for every position in the reference.'''

        # We'll return the scores for all true positions, but not for
        # the dummy entry [0]. See the comments in fillMatrix for an
        # explanation.

        # Keep in mind that the score for a given position applies to
        # an alignment which *ends* at that position. There is no
        # unambiguous way to assign a score to a starting position,
        # since backtracks from multiple end positions could arrive at
        # the same start position. In that case, which score do you
        # choose?

        if not self._matrixValid:
            raise RuntimeError("matrix has not been computed.")

        row = self._lastRow
        lastCol = len(self._ref)+1                    # matrix may be larger than current ref
        return tuple(self._matrix[row, 1:lastCol])    # skip the dummy entry [0]

    def peakPosits (self, scores=None, thresh1=None, thresh2=None):
        '''Return a list containing positions of peaks in a list of scores.'''

        if not self._matrixValid:
            raise RuntimeError("matrix has not been computed.")

        # A contiguous run of scores > thresh1 defines a 'range'. In
        # each range, report the highest score as the peak. Don't
        # start a new range until scores drop below thresh2, to avoid
        # calling lots of peaks when scores are bouncing around near
        # thresh1.

        # What constitutes a good score (and hence a good choice of a
        # threshold) is a function of the length of the query read. By
        # default, we take the scores from the current matrix, and can
        # produce a good guess at the thresholds. If the caller
        # supplies the scores, he must also supply at least
        # thresh1. We can make a decent guess at thresh2 from that, if
        # it's not supplied.

        if scores == None:
            scores = self.allScores()
            if thresh1 == None:
                thresh1 = int(len(self._read) * self._match * 0.4)
        elif thresh1 == None:
            raise RuntimeError("if scores are supplied, thresh1 must be supplied.")

        if thresh2 == None:
            thresh2 = thresh1 - 3

        thresh3 = 0

        curMax = 0
        curIx  = 0
        peaks = []

        for ix in xrange(len(scores)):

            if scores[ix] >= thresh1:
                if scores[ix] > curMax:
                    curMax = scores[ix]
                    curIx  = ix
                    thresh3 = int(curMax * 0.8)
####            elif scores[ix] < thresh2:
            elif scores[ix] < thresh3:
                if curMax != 0:
                    peaks.append(curIx)
                    curMax = 0

        if curMax != 0:                # last one?
            peaks.append(curIx)

        return peaks

# Copyright (C) 2011 Genome Research Limited
#
# This library is free software. You can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
