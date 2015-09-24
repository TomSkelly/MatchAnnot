#!/usr/bin/env python

# Read a fasta file of reads, and a fasta file of primers. Find
# multiple occurrences of primers at the end of reads, and chop them
# off. Write a new fasta file.

# Some reads are so scrambled that it isn't obvious what to trim, or
# that they contain any useful sequence at all. A conservative
# approach is taken: If the fix isn't obvious, keep the read as it
# is. Later I may add an option to discard such bad reads.

# We will follow the IsoSeq convention that 5' primers are labelled
# >F0, >F1, ... in the primer fasta file, and 3' primers (including
# RACE) are >R0, >R1, ...  However, given the way SMRTbell sequencing
# and consensus read generation work, the read to be trimmed can be in
# either forward or reverse orientation. So two scenarios must be
# examined:

# Forward:

# F0+  >>>  >>>  >>>  |>>>                          R0+ >>>|  >>>  >>>  >>>
# read ===============|===================================================

# Reverse:
# R0-  <<<  <<<  <<<  |<<<                          F0- <<<|  <<<  <<<  <<<
# read ===============|===================================================

# where the trim point is indicated by |. Note that the SWAligner
# module returns the position of the LAST base of a match, so in the
# F0+ and R0- cases we will need to back up to the trim point.

# We will keep the read in its original sense, and try both
# orientations of each primer, as shown. Note that we do not
# accommodate 'backwards' cases, such as an F0- primer at the front of
# the read.

# For efficiency reasons, only a sub-window of bases at the start and
# end of the read are examined. The window size is a command-line
# parameter.

import os
import sys
import optparse
import re        # regular expressions
import string

from tt_log import logger
import Reference as ref
import SWAligner

VERSION = '20150406.01'

TICK = 100
FASTA_WRAP = 70
MIN_SCORE  = 15
THRESH_VALUE = 0.6
COMPLTAB   = string.maketrans ('ACGTacgt', 'TGCAtgca')         # for reverse-complementing reads

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    if opt.format == 'pickle':
        reads = ref.Reference.fromPickle (opt.input)
    else:
        reads = ref.Reference (opt.input)

    primers = ref.Reference (opt.primers)

    aln = SWAligner.Aligner()
    
    report = None
    if opt.report is not None:
        report = open (opt.report, 'w')

    numReads = 0                                   # counter

    for chr in reads.chromosomes():                # process all the reads

        readIn = reads.sequence(chr)
        aln.setRef (readIn)
        numReads += 1
        if numReads % TICK == 0:
            logger.debug('tick %6d %s' % (numReads, chr))

        if report is not None:
            report.write ('>%s  %d\n' % (chr, len(readIn)))

        for pri in primers.chromosomes():          # for each read, look at all the primers

            forPri = primers.sequence(pri)
            revPri = forPri[::-1].translate(COMPLTAB)

            if pri[0] == 'F':

                trim = tryPrimer ('-', forPri, pri, aln, report)
                if trim > 0:
                    readIn = readIn[ trim : ]         # trim the read
                    aln.setRef (readIn)               # update the aligner's target read if changed

                trim = tryPrimer ('+', revPri, pri, aln, report)
                if trim > 0:
                    readIn = readIn[ : trim ]
                    aln.setRef (readIn)
                    
            elif pri[0] == 'R':

                trim = tryPrimer ('+', forPri, pri, aln, report)
                if trim > 0:
                    readIn = readIn[ : trim ]
                    aln.setRef (readIn)
                    
                trim = tryPrimer ('-', revPri, pri, aln, report)
                if trim > 0:
                    readIn = readIn[ trim : ]
                    aln.setRef (readIn)

            else:
                raise RuntimeError ('primer name %s does not begin with F or R' % pri)

        writeFasta (chr, readIn)

    if report is not None:
        report.close()

    logger.debug('finished')

    return

def tryPrimer (dir, primerIn, primerName, aln, report):

    aln.setRead (primerIn)
    score = aln.fillMatrix()
    trim = 0

    if score >= MIN_SCORE:

        allScores = aln.allScores()
        matchScore = aln.getPenalties()[0]             # score increment for a match
        peakThresh = int(len(primerIn) * matchScore * THRESH_VALUE)
        peakList = list()

        for peak in aln.peakPosits(thresh1=peakThresh):

            peakList.append (peak)

            if report is not None:
                readStr, primerStr = aln.alignmentStrings(pos=peak+1)
                report.write ( '%-10s  %s  %7d  %3d  %s\n' % (primerName, dir, peak, allScores[peak], readStr))
                report.write ( '                             %s\n' % primerStr)

        if len(peakList) > 1:

            if dir == '+':

                trim = peakList[0] + 1               # +1: we want to keep readIn[trim]
                if report is not None:
                    report.write ('clipped %s+ at %d\n' % (primerName, trim))

            else:

                readStr, primerStr = aln.alignmentStrings(pos=peakList[-1]+1)
                backup = len(readStr) - readStr.count('-')
                trim = peakList[-1] - backup + 1
                if report is not None:
                    report.write ('clipped %s- at %d  %d\n' % (primerName, trim, backup))

    return trim

def writeFasta (chr, readIn):

    print '>' + chr

    for ix in xrange(0, len(readIn), FASTA_WRAP):
        print readIn[ix:ix+FASTA_WRAP]

    return

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options]', version=VERSION)

    parser.add_option ('--input',     help='reference file, in format specified by --format')
    parser.add_option ('--format',    help='format of reference file (def: %default)', \
                           type='choice', choices=['standard', 'fasta', 'pickle'])      # 'standard' is the same as 'fasta'
    parser.add_option ('--primers',   help='fasta file containing primers (required)')
    parser.add_option ('--report' ,   help='output report file name (optional)')

    parser.set_defaults (format='fasta',
                         )

    opt, args = parser.parse_args()

    return opt, args


if __name__ == "__main__":
    main()
