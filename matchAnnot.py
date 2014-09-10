#!/usr/bin/env python

# Read a .sam file produced by STAR-aligning IsoSeq isoforms with the
# human genomic reference. Read the GENCODE .gtf annotations file. For
# each isoform, print the genes, transcripts and exons it overlaps.

# Both the SAM and annotation files are assumed to be sorted by chr
# and position.

# AUTHOR: Tom Skelly (thomas.skelly@fnlcr.nih.gov)

import os
import sys
import optparse
import re        # regular expressions

from tt_log import logger
import Annotations as anno
import Best        as best
import Clusters    as cl
import CigarString as cs
import PolyA

VERSION = '20140909.01'

FLAG_NOT_ALIGNED = 0x04         # SAM file flags
FLAG_REVERSE     = 0x10
FLAG_SECONDARY   = 0x100

POLYA_REACH      = 30           # how far from 3' end to look for poly-A motif

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    if opt.clusters is not None:
        clusterList = cl.ClusterList (opt.clusters)

    annotList   = anno.AnnotationList (opt.gtf)
    annotCursor = anno.AnnotationCursor (annotList)

    polyAFinder = PolyA.PolyA()

    regexCluster = re.compile('(c\d+)')      # isoform name format varies, but cnnnn should be in there somewhere
    regexAS = re.compile('(AS:i:\d+)')       # alignment score
    regexUT = re.compile('(uT:A:\d+)')       # mismatch reason (ToDo: Translate this)

    if len(args) > 0:
        logger.debug('reading SAM file %s' % args[0])
        handle = open (args[0], 'r')
    else:
        logger.debug('reading SAM data from stdin')
        handle = sys.stdin

    totReads = 0                  # assorted counters
    totAlign = 0
    totWithGene = 0
    totMulti = 0
    totReverse = 0
    totByScore = [0,0,0,0,0,0]    # indexed by score
    totSplice  = [0,0,0,0,0,0]

    lastPos = dict()              # current position in SAM file by chr (for sort check)

    for line in handle:           # main loop: read the SAM file

        if line.startswith('@'):
            continue

        totReads += 1

        lineFields = line.split('\t')                          # just split it once, not 6 times in list comp
        clusterName, flags, chr, start, cigarString, bases = [lineFields[i] for i in (0,1,2,3,5,9)]
        flags = int(flags)

        if flags & FLAG_NOT_ALIGNED:

            print '\nisoform:  %-16s' % (clusterName)

            if opt.clusters is not None:                       # print cluster (cl:) lines
                clusterID = re.search(regexCluster, clusterName).group(1)
                for FL, cellNo, reads in clusterList.showReads(clusterID):
                    print 'cl:      %2d %-5s ' % (cellNo, FL), '  '.join(['%-16s' % x for x in reads])

            alnScore  = re.search(regexAS, line).group(1)      # alignment score
            alnReason = re.search(regexUT, line).group(1)      # mismatch reason

            print 'result:   %-50s no_alignment_found  %s  %s' % (clusterName, alnReason, alnScore)

            continue

        totAlign += 1

        strand = '-' if (flags & FLAG_REVERSE) else '+'
        cigar = cs.CigarString(cigarString)
        start = int(start)
        end   = start + cigar.genomicLength() - 1;       # -1 to report last base, rather than last+1

        if start < lastPos.get(chr, 0):
            raise RuntimeError ('SAM file is not sorted by position')
        lastPos[chr] = start

        exons = cigar.exons (start)

        print '\nisoform:  %-16s    %9d                    %9d         %-5s  %s  %6d   %dS+%dS' \
            % (clusterName, start, end, chr, strand, end-start, cigar.leading, cigar.trailing),
        if flags & FLAG_SECONDARY:
            print ' multimap',
            totMulti += 1
        print

        print 'cigar:    %s' % cigarString

        if opt.clusters is not None:                     # print cluster (cl:) lines
            clusterID = re.search(regexCluster, clusterName).group(1)
            for FL, cellNo, reads in clusterList.showReads(clusterID):
                print 'cl:      %2d %-5s ' % (cellNo, FL), '  '.join(['%-16s' % x for x in reads])

        foundPolyA = False
        print 'polyA:  ',                                # print 'polyA:' line
        for motif, offset in polyAFinder.findMotifs (bases, strand, POLYA_REACH):
            print ' %s: %4d' % (motif, offset),
            foundPolyA = True
        print

        # Now let's do the genes...

        # Bump cursor past any genes which end prior to the start of
        # the current read. We're done looking at them.

        annotCursor.advance (chr, start)

        bestHit = best.Best()

        # Loop through genes this cluster overlaps. Try the aligned
        # strand first, but if no joy, try the other strand, since
        # IsoSeq can get it backwards sometimes.

        for str2try in (strand, '-' if strand == '+' else '+'):        # try aligned strand first

            for curGene in annotCursor.getOverlappingGenes (chr, start, end, str2try):

                print 'gene:     %-16s    %9d            %6d             %9d  %5d     %s' \
                    % (curGene.name, curGene.start, curGene.start-start, curGene.end, curGene.end-end, curGene.strand),
                if str2try != strand:
                    print '  rev',
                print

                bestTran, bestScore = matchTranscripts (exons, curGene)     # match this cluster to all transcripts of gene
                bestHit.update (bestScore, [curGene, bestTran])             # best transcript of best gene so far?

            if bestHit.which > 1:                                           # if decent match found, don't try the other strand
                break

        if bestHit.value is None:
            print 'result:   %-50s no_genes_found' % (clusterName)
        else:

            bestGene, bestTran = bestHit.which
            print 'result:   %-50s  %-20s  %-24s  ex: %2d  sc: %d' % (clusterName, bestGene.name, bestTran.name, len(exons), bestHit.value),

            if bestGene.strand != strand:                                   # if best hit was on other strand from alignment
                print 'rev',
                totReverse += 1

            if bestHit.value >= 3:
                delta5 = bestTran[0].start - exons[0].start              # 5' delta
                delta3 = bestTran[-1].end  - exons[-1].end               # 3' delta
                print ' 5-3: %5d %5d' % (delta5, delta3),
            print

            totWithGene += 1
            totByScore[bestHit.value] += 1
            if foundPolyA:
                totSplice[bestHit.value] += 1

    print  '\nsummary: version %s\n' % VERSION

    if opt.clusters is not None:
        for cellNo, cell in clusterList.showCells():
            print 'summary:   cell %d = %s' % (cellNo, cell)
        print

    print 'summary: %7d isoforms read' % totReads
    print 'summary: %7d isoforms aligned, of which %d were multiply mapped' % (totAlign, totMulti)
    print 'summary: %7d isoforms hit at least one gene, of which %d were on opposite strand' % (totWithGene, totReverse)
    print

    for score in xrange(5,-1,-1):
        print 'summary: %7d isoforms scored %d, of which %6d had splice termination motif' \
            % (totByScore[score], score, totSplice[score])

    logger.debug('finished')

def matchTranscripts (exons, gene):
    '''
    Given the list of annotated transcripts for a gene of interest,
    compare our isoform's exons to the exons of each transcript. Look for
    case(s) where there is a one-for-one match.

    The logic here is a bit cumbersome, because we need to loop over the
    transcripts twice: first, to find the highest-scorer, then again to
    print them.
    '''

    strings = list()         # overlaps are computed in the first pass. save them (as strings) to print in pass 2.
    scores  = list()         # scores saved for the print step

    bestScore = best.Best()
    bestTrunc = best.Best(reverse=True)

    for ix, tran in enumerate (gene):

        tranExons = tran.children                                     # the list of exons for this transcript

        over1, over2 = findOverlaps (exons, tranExons)
        strings.append(overlap2string(over1))                         # save it for the print step as a string

        if not isMatch (over1):                                       # do leading exons (at least) match?
            score = 1
        elif len(over1) != len(over2):                                # do *all* exons match (i.e., overlap)?
            score = 2
        elif not internalMatch (exons, tranExons):                    # do their sizes match?
            score = 3
        else:
            score = 4                                                 # may get promoted later
            trunc = abs(exons[0].start - tranExons[0].start) \
                +   abs(exons[-1].end  - tranExons[-1].end)           # leading and trailing exon truncation amount
            bestTrunc.update (trunc, ix)                              # keep track of score-4 with the smallest UTR size difference

        bestScore.update (score, ix)
        scores.append(score)                                          # save for the print step

    if bestScore.value == 4:                                          # promote the best 4 (if there is one) to a 5
        scores[bestTrunc.which] = 5
        bestScore.update (5, bestTrunc.which)

    # Now we print the results:

####    for ix in xrange(gene.numChildren()):
####        tran = gene[ix]

    for ix, tran in enumerate (gene):

        print 'tr:       %-20s  sc: %d  ex: %2d  %5d  id: %-20s   %s' \
            % (tran.name, scores[ix], tran.numChildren(), tran.length, tran.ID, strings[ix])

####        if scores[ix] >= 3:
        showCoords (exons, tran)

    return gene[bestScore.which], scores[bestScore.which]     # return best transcript and score


def findOverlaps (list1, list2):
    '''
    Given two lists of intervals, find the overlaps between them.  The
    lists contain objects which have 'start' and 'end' attributes --
    they need not be of the same type, or suclasses of each other.

    We return two arrays, one for each input list. There is an entry
    in the output array corresponding to each input list entry. The
    output entry is an array of indices of the entries in the other
    list which overlap this input entry.

    Note that we could get by here with generating only one of the two
    overlap lists. It is easy enough to derive one from the
    other. E.g.:

    Given:  [0] [1] [] [2,3] [] [4] [6]
    Invert: [0] [1] [3] [3] [5] [] [6]

    But that wouldn't save much. And this way the interface is nice
    and symmetrical.
    '''

    over1 = [ [] for x in list1 ]
    over2 = [ [] for x in list2 ]
    pos2  = 0

    for pos1 in xrange(len(list1)):     # for each list1 entry, find all list2 which overlap it

        while pos2 < len(list2) and list2[pos2].end < list1[pos1].start:
            pos2 += 1                   # we're done with this list2 entry

        ix = pos2                       # start from current list2 entry and look forward
        while ix < len(list2) and list2[ix].start < list1[pos1].end:
            over1[pos1].append(ix)
            over2[ix].append(pos1)
            ix += 1

        pos1 += 1

    return over1, over2


def overlap2string (list1):
    '''Create a printable string representation of an overlap list.'''

    # As of version 20140903.01, exons in *printed* output are
    # numbered 1..N, rather than 0..N-1.

    fields = list()

    for sublist in list1:
        strList =  [str(x+1) for x in sublist]                # join doesn't do ints
        fields.append ('[' + ','.join(strList) + ']')

    return ' '.join(fields)


def isMatch (list1):
    '''Test whether an interval overlap list hits every interval.'''

    # Note that if the supplied list is shorter than the list it was
    # compared to, it will pass. I.e., the following situation is a
    # pass:

    #    list-1:   [0] [1] [2] [3] [4]
    #    list-2:   [0] [1] [2] [3] [4] [] [] []

    for ix in xrange(len(list1)):
        if len(list1[ix]) != 1 or list1[ix][0] != ix:
            return False

    return True


def internalMatch (list1, list2):
    '''
    Given two lists of exon intervals which overlap one-for-one,
    determine whether their coordinates match exactly, EXCEPT for the
    start of the first exon and the end of the last exon. That
    variation is currently thought to be caused (at the 3' end, at
    least) by polyadenylation at more sites than the annotations
    record.
    '''
    
    for ix in xrange(1,len(list1)):                    # check starts, skip exon 0
        if list1[ix].start != list2[ix].start:
            return False

    for ix in xrange(len(list1)-1):                    # check ends, skip last exon
        if list1[ix].end != list2[ix].end:
            return False

    return True


def canPrint (over):
    '''Determine whether a set of overlapping exons can be meaningfully printed.'''

    last = -1                              # internally, exons are numbered 0..N-1
    hit  = False

    for group in over:

        if len(group) > 1:                 # can't do [3,4]
            return False

        if  len(group) == 1:
            if group[0] <= last:
                return False               # can't do [4] [4]
            last = group[0]
            hit = True                     # we found at least one overlap

    return hit

def showCoords (readExons, tranExons):
    '''
    Print coordinates of matched exons, given their interval
    lists. Note that this only works when the exons match
    one-for-one. The first list is assumed to also contain counts of
    insertions and deletions.
    '''

    overR, overT = findOverlaps (readExons, tranExons)      # already done in matchTranscripts, but easier to recompute than to save
    if not canPrint (overR):                                # check that there are no multi-exon overlaps
        return

    ixR = 0
    ixT = 0

    while ixR < len(overR) and ixT < len(overT):

        if len(overR[ixR]) == 0 and len(overT[ixT]) == 0:

            if readExons[ixR].start < tranExons[ixT].start:

                printReadExon (ixR, readExons[ixR])
                ixR += 1

            else:  

                printTranExon (ixT, tranExons[ixT])
                printStartStop (tranExons, tranExons[ixT])
                ixT += 1

        elif len(overR[ixR]) == 0:

            printReadExon (ixR, readExons[ixR])
            ixR += 1

        elif len(overT[ixT]) == 0:

            printTranExon (ixT, tranExons[ixT])
            printStartStop (tranExons, tranExons[ixT])
            ixT += 1

        elif overR[ixR][0] == ixT and overT[ixT][0] == ixR:

            printMatchingExons (ixR, ixT, readExons[ixR], tranExons[ixT])
            printStartStop (tranExons, tranExons[ixT])
            ixR += 1
            ixT += 1

        else:
            raise RuntimeError ('improper overlap: %s <-> %s' % (overlap2string(overR), overlap2string(overT)))

    while ixR < len(overR):

        exonR = readExons[ixR]
        printReadExon (ixR, exonR)
        ixR += 1

    while ixT < len(overT):

        exonT = tranExons[ixT]
        printTranExon (ixT, exonT)
        printStartStop (tranExons, exonT)
        ixT += 1


def printMatchingExons (ixR, ixT, exonR, exonT):

    print 'exon:                %2d  %2d   %9d  %9d  %5d  %9d  %9d  %5d      len: %4d %4d  ins: %2d  del: %2d' \
        % (ixR+1, ixT+1, \
               exonR.start, exonT.start, exonT.start-exonR.start, \
               exonR.end,   exonT.end,   exonT.end-exonR.end, \
               exonR.end-exonR.start+1,  exonT.end-exonT.start+1, \
               exonR.inserts, exonR.deletes),  # comma: line continued below


def printReadExon (ixR, exonR):

    print 'exon:                %2d   .   %9d          .      .  %9d          .      .      len: %4d    .  ins: %2d  del: %2d' \
        % (ixR+1, exonR.start, exonR.end, exonR.end-exonR.start+1, exonR.inserts, exonR.deletes)        # no comma: EOL here

def printTranExon (ixT, exonT):

    print 'exon:                 .  %2d           .  %9d      .          .  %9d      .      len:    . %4d' \
        % (ixT+1, exonT.start, exonT.end, exonT.end-exonT.start+1),  # comma: line continued below

def printStartStop (tranExons, exonT):

    if hasattr (tranExons, 'startcodon'):
        if exonT.start <= tranExons.startcodon and exonT.end >= tranExons.startcodon:
            print '  start: %4d' % (tranExons.startcodon - exonT.start),        # distance from start of exon
    if hasattr (tranExons, 'stopcodon'):
        if exonT.start <= tranExons.stopcodon and exonT.end >= tranExons.stopcodon:
            print '  stop:  %4d' % (tranExons.stopcodon  - exonT.end),          # distance from end of exon
    print

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <SAM_file> ... ')

    parser.add_option ('--gtf',      help='annotations in gtf format (required)')
    parser.add_option ('--clusters', help='cluster_report.csv file name (optional)')

    parser.set_defaults (gtf=None,
                         clusters=None,
                         )

    opt, args = parser.parse_args()

    return opt, args


if __name__ == "__main__":
    main()
