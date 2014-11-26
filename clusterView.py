#!/usr/bin/env python

# Read annotation file, plot exons for each transcript of a specified gene.

import os
import sys
import optparse
import re        # regular expressions

from tt_log import logger

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import Annotations as anno
import Best        as best
import Cluster     as cl
import CigarString as cs

VERSION = '20141125.01'

DEF_OUTPUT = 'exons.png'        # default plot filename
DEF_YSCALE = 1.0                # default  Y-axis scale factor

FIG_WIDTH = 14
FIG_HEIGHT_PER_TRANS = 0.2      # figure height depends on the number of rows
MAX_LABELS = 20                 # how many labels fit across the X axis

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    # Find all the exons in all the transcripts for the gene, put them
    # in a list.

    tranList = list()                                      # list of Transcript objects
    exonList = list()                                      # list of Exon objects

    if opt.gtf is not None:
        getGeneFromAnnotation (opt, tranList, exonList)    # lists will be changed
    if opt.matches is not None:
        getGeneFromMatches (opt, tranList, exonList)       # lists will be changed
    if len(exonList) == 0:
        raise RuntimeError ('no exons found for gene %s in annotation or match files' % opt.gene)

    forwardStrand = '-' if opt.flip else '+'
    if exonList[0].strand == forwardStrand:
        exonList.sort(key=lambda x: x.start)               # sort the list by start position
        blocks = assignBlocks (opt, exonList)              # assign each exon to a block
    else:
        exonList.sort(key=lambda x: x.end, reverse=True)   # sort the list by decreasing end position
        blocks = assignBlocksReverse (opt, exonList)       # assign each exon to a block -- backwards

    tranNames = orderTranscripts (tranList)      # set Y-axis coords for the transcripts

    plt.figure (figsize=(FIG_WIDTH, int(FIG_HEIGHT_PER_TRANS*len(tranList)*opt.yscale)))

    plotExons (exonList, blocks)                 # plot the exons

    plotStartStop (tranList, blocks)             # start/stop  codons to plot

    plotBoundaries (tranNames, blocks)           # plot block boundaries as vertical lines

    if opt.title is not None:
        plt.title (opt.title)
    else:
        plt.title ('transcripts for gene %s' % opt.gene)

    plt.savefig (opt.output)
    plt.close()

    if opt.details is not None:
        printDetails (opt, blocks, exonList)

    logger.debug('finished')

    return

def getGeneFromAnnotation (opt, tranList, exonList):
    '''Add to lists of transcripts and exons: annotations for gene of interest.'''

    if opt.gtf == None:
        return tranList, exonList

    if opt.format == 'pickle':
        annotList   = anno.AnnotationList.fromPickle (opt.gtf)
    elif opt.format == 'alt':
        annotList   = anno.AnnotationList (opt.gtf, altFormat=True)
    else:     # standard format
        annotList   = anno.AnnotationList (opt.gtf)

    allGenes = annotList.getGeneDict()
    if opt.gene not in allGenes:
        raise RuntimeError ('gene %s is not in the annotation file' % opt.gene)
    geneList = allGenes[opt.gene]       # a list of Annotation objects
    if len(geneList) > 1:
        logger.warning('gene %s appears %d times in annotations, first occurrence plotted' \
                           % (opt.gene, len(geneList)))
    myGene = geneList[0]

    for tran in myGene.getChildren():                       # tran is an Annotation object

        myTran = Transcript(tran.name, annot=True)

        if hasattr(tran, 'startcodon'):
            myTran.startcodon = tran.startcodon
        if hasattr(tran, 'stopcodon'):
            myTran.stopcodon = tran.stopcodon

        for exon in tran.getChildren():                     # exon is an Annotation object
            myExon = Exon(myTran, exon.name, exon.start, exon.end, exon.strand)
            exonList.append (myExon)
            myTran.exons.append(myExon)

        tranList.append (myTran)

    return tranList, exonList

def getGeneFromMatches (opt, tranList, exonList):
    '''Add to lists of transcripts and exons: clusters which matched gene of interest.'''

    if opt.matches == None:
        return tranList, exonList

    regex = re.compile ('(c\d+)')                                      # regex for cluster ID
    omits = [] if opt.omit is None else opt.omit.split(',')            # clusters which must not be included
    shows = [] if opt.show is None else opt.show.split(',')            # clusters which must be included

    localList = list()                                                 # temporary list of clusters
    totClusters = 0

    clusterDict = cl.ClusterDict.fromPickle (opt.matches)

    for cluster in clusterDict.getClustersForGene(opt.gene):           # cluster is Cluster object

        totClusters += 1

        match = re.search (regex, cluster.name)
        if match is not None and match.group(1) in shows:              # if this is a force-include cluster
            localList.append ( [cluster, 'f999999p999999'] )           # fake sort key to push it to the front
        elif match is None or match.group(1) not in omits:
            full, partial = cluster.getFP()
            sortKey = 'f%06dp%06d' % (full, partial)                   # single key includes full and partial counts
            localList.append ( [cluster, sortKey] )

    if opt.howmany is not None:
        localList.sort(key=lambda x: x[1], reverse=True)               # sort by full/partial counts
        localList = localList[:opt.howmany]                            # keep the top N entries (which will include the forces)

    for ent in localList:

        cluster = ent[0]
        myTran = Transcript(cluster.name, score=cluster.bestScore)

        cigar = cs.CigarString(cluster.cigarString)
        for exonNum, exon in enumerate(cigar.exons (cluster.start)):   # exon is a cs.Exon object
            exonName = '%s/%d' % (myTran.name, exonNum)                # exons don't have names: make one up
            myExon = Exon(myTran, exonName, exon.start, exon.end, cluster.strand)
            exonList.append (myExon)
            myTran.exons.append(myExon)

        tranList.append (myTran)

    logger.debug('kept %d of %d clusters' % (len(localList), totClusters))

    return tranList, exonList

def assignBlocks (opt, exonList):
    '''
    Assign exons to blocks, separated by sequence which is intronic in
    all transcripts. exonList is assumed to be sorted by ascending
    start position.
    '''

    adjust  = 0
    blockNo = 0
    exonIx  = 0
    blocks = list()

    while exonIx < len(exonList):               # can't use enumerate here, it's a double loop

        blockStart = exonList[exonIx].start     # block start = start of first exon in block
        blockEnd   = exonList[exonIx].end       # initial value, updated in the loop below
        blockStartIx = exonIx

        while exonIx < len(exonList) and exonList[exonIx].start <= blockEnd:
            myExon = exonList[exonIx]
            if myExon.end > blockEnd:
                blockEnd = myExon.end
            myExon.block = blockNo
            myExon.tran.blocks.add(blockNo)    # transcript has an exon in this block
            myExon.adjStart = myExon.start - blockStart + adjust
            exonIx += 1

        adjust += blockEnd - blockStart + 1
        blocks.append(Block(blockStart, blockEnd, adjust))
        blockNo += 1

    return blocks

def assignBlocksReverse (opt, exonList):
    '''
    Like assignblocks, but for the reverse strand, ordering blocks
    from the 5' end of the transcript. exonList is assumed to be
    sorted by decreasing exon end position.
    '''

    # I did this as a separate mirror image of assignBlocks, rather
    # than clutter the scenery with lots of forward/reverse checks.

    adjust  = 0
    blockNo = 0
    exonIx  = 0
    blocks = list()

    while exonIx < len(exonList):               # can't use enumerate here, it's a double loop

        blockStart = exonList[exonIx].end       # block start = end of last exon in block
        blockEnd   = exonList[exonIx].start     # initial value, updated in the loop below
        blockStartIx = exonIx

        while exonIx < len(exonList) and exonList[exonIx].end >= blockEnd:

            myExon = exonList[exonIx]
            if myExon.start < blockEnd:
                blockEnd = myExon.start
            myExon.block = blockNo
            myExon.tran.blocks.add(blockNo)    # transcript has an exon in this block
            myExon.adjStart = blockStart - myExon.end + adjust
            exonIx += 1

        adjust += blockStart - blockEnd + 1
        blocks.append(Block(blockStart, blockEnd, adjust))
        blockNo += 1

    return blocks

def orderTranscripts (tranList):
    '''
    Order the transcripts (i,e., assign each a Y coordinate) so similar
    transcripts are close to each other.
    '''

    # The measure of similarity used here is block occupancy: The
    # distance between two transcripts is the number of blocks where
    # one transcript has exons, and the other doesn't. How many exons
    # there are, or how similar they are in length, is not looked at.

    # The ordering is done using a greedy nearest-neighbor
    # heuristic. To do it optimally turns it into a Traveling Salesman
    # problem.

    tranNames = list()
    curTran = tranList[0]            # arbitrarily start with the first transcript
    tranIx = 0

    while True:                                     # loop until break below

        tranNames.append(curTran.name)              # needed for yticks call
        curTran.tranIx = tranIx

        bestTran = best.Best(reverse=True)
        for myTran in tranList:                     # find the next closest transcript
            if myTran.tranIx is None:               # if transcript hasn't been indexed yet
                diff = len(curTran.blocks.symmetric_difference(myTran.blocks))
                bestTran.update(diff, myTran)

        if bestTran.which is None:                  # every transcript has its index: we're done
            break

        curTran = bestTran.which
        tranIx += 1

    return tranNames

def plotExons (exonList, blocks):
    '''Plot exons.'''

    for myExon in exonList:

        exonSize = myExon.end - myExon.start + 1
        adjStart = myExon.adjStart

        if myExon.tran.annot:                       # exon from annotation file?
            color = 'green'
            blocks[myExon.block].annot = True       # this block includes annotation
        else:
            color = 'purple' if myExon.tran.score == 5 else 'blue'

        plt.hlines (myExon.tran.tranIx+1, adjStart, adjStart+exonSize, linewidth=3, colors=color, zorder=1)

    return

def plotStartStop (tranList, blocks):
    '''Add start/stop codons to plot.'''

    for tran in tranList:
        if tran.annot:                             # only annotations know about start/stops
            if hasattr(tran, 'startcodon'):
                plotCodon (tran, tran.startcodon, blocks, 'chartreuse')
            if hasattr(tran, 'stopcodon'):
                plotCodon (tran, tran.stopcodon, blocks, 'red')

    return

def plotCodon (tran, posit, blocks, color):
    '''Add a codon mark to the plot.'''

    for blk in blocks:

        if blk.start <= posit and blk.end >= posit or \
                blk.start >= posit and blk.end <= posit:      # check in both strand directions

            xPos = blk.boundary - abs(blk.end-posit)
            plt.scatter (xPos, tran.tranIx+1, s=25, c=color, marker='v', zorder=2)
            break

    return

def plotBoundaries (tranNames, blocks):
    '''Plot exon boundaries as vertical lines.'''

    ticksStart  = list()     # x-axis tick positions in phony space
    ticksEnd    = list()
    labelsStart = list()     # x-axis labels are real genomic coordinates of block start
    labelsEnd   = list()

    tickSep = blocks[-1].boundary / MAX_LABELS     # separation required to avoid label overlap

    frompos = 0
    for bound in blocks:

        plt.vlines (bound.boundary, 0, len(tranNames)+1, linewidth=0.5, colors='red')    # block boundary

        if not bound.annot:            # if block does not include annotation exon(s), make it blush
            plt.axvspan(frompos, bound.boundary, facecolor='pink', alpha=0.8)

        if len(ticksStart) == 0 or frompos - ticksStart[-1] > tickSep:     # add tick only if there's room for the label
            ticksStart.append(frompos)
            labelsStart.append(str(bound.start))
        if len(ticksEnd) == 0 or bound.boundary - ticksEnd[-1] > tickSep:
            ticksEnd.append(bound.boundary)
            labelsEnd.append(str(bound.end))

        frompos = bound.boundary       # new block start

    plt.grid(axis='y')                     # turn on horizontal dotted lines
    plt.xlim (0, blocks[-1].boundary)
    plt.ylim (len(tranNames)+1, 0)
    plt.xticks (ticksStart, labelsStart, fontsize='xx-small')
    plt.yticks (xrange(1,len(tranNames)+2), tranNames, fontsize='xx-small')

    # These lines came from an obscure posting on StackOverflow. They
    # create a second X axis at the bottom of the figure. I haven't
    # got a clue how they work, but they do.

    x2 = plt.twiny()
    x2.set_frame_on(True)
    x2.patch.set_visible(False)
    x2.xaxis.set_ticks_position('bottom')
    x2.xaxis.set_label_position('bottom')
    x2.spines['bottom'].set_position(('outward', 20))

    plt.axes(x2)

    plt.grid(axis='y')
    plt.xlim (0, blocks[-1].boundary)
    plt.ylim (len(tranNames)+1, 0)
    plt.xticks (ticksEnd, labelsEnd, fontsize='xx-small')
    plt.yticks (xrange(1,len(tranNames)+2), tranNames, fontsize='xx-small')

    return

def printDetails (opt, blocks, exonList):
    '''Print numerical details on block occupancy.'''

    exonIx = 0
    lastBoundary = 0

    handle = open (opt.details, 'w')

    for ix, blk in enumerate(blocks):

        blockSize = blk.boundary - lastBoundary
        lastBoundary = blk.boundary
        handle.write ('\nblock %d (%d)  %d  %d:\n' % (ix, blockSize, blk.boundary, blk.start))

        while exonIx < len(exonList) and exonList[exonIx].block == ix:
            myExon = exonList[exonIx]
            exonSize = myExon.end-myExon.start+1
            handle.write ('    %-32s  %s  %9d  %9d  %5d  %5d  %d\n' \
                              % (myExon.name, myExon.strand, myExon.start, myExon.end, \
                                     exonSize, myExon.adjStart, myExon.adjStart+exonSize))
            exonIx += 1

    handle.close()

    return

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options]')

    parser.add_option ('--gtf',     help='annotations file, in format specified by --format (optional)')
    parser.add_option ('--format',  help='format of annotation file: standard, alt, pickle (def: %default)', \
                           type='choice', choices=['standard', 'alt', 'pickle'])
    parser.add_option ('--matches', help='pickle file from matchAnnot.py (optional)')
    parser.add_option ('--gene',    help='gene to plot (required)')
    parser.add_option ('--omit',    help='clusters to ignore, e.g., c1234,c2345,c3456')
    parser.add_option ('--show',    help='clusters to force shown, even if underpopulated')
    parser.add_option ('--howmany', help='how many clusters to plot (def: all)', type='int')
    parser.add_option ('--output',  help='output plot file name (def: %default)')
    parser.add_option ('--flip',    help='reverse plot orientation (def: mRNA 5\' on left)', action='store_true')
    parser.add_option ('--yscale',  help='amount by which to scale Y axis (def: %default)', type='float')
    parser.add_option ('--details', help='output file name for details in text format (def: no output)')
    parser.add_option ('--title',   help='title for top of figure')

    parser.set_defaults (format='standard',
                         output=DEF_OUTPUT,
                         yscale=DEF_YSCALE,
                         )

    opt, args = parser.parse_args()

    return opt, args


class Transcript (object):
    '''Just a struct actually, containing data about a transcript.'''

    def __init__ (self, name, score=None, annot=False):
        
        self.name   = name
        self.score  = score
        self.annot  = annot            # transcript comes from annotations?
        self.tranIx = None             # y-axis coordinate of transcript
        self.blocks = set()            # blocks where this transcript has exon(s)
        self.exons  = list()           # Exon objects for this transcript


class Exon (object):
    '''Struct containing data about an exon.'''

    def __init__ (self, tran, name, start, end, strand):

        self.tran     = tran           # Transcript object containing this exon
        self.name     = name
        self.start    = start
        self.end      = end
        self.strand   = strand
        self.block    = None           # block number where this exon resides
        self.adjStart = None           # start of exon in phony x-axis coordinates


class Block (object):
    '''Struct for plot block.'''

    # A plot block is a vertical span representing a range contiguous
    # bases. Plot blocks are separated by vertical lines representing
    # regions of the reference, of unspecified length, which contain
    # no exons.

    # one implication of that scheme is that the x axis of the plot is
    # meaningless: it represents neither genomic nor RNA sequence range.

    def __init__ (self, start, end, boundary):

        self.start    = start          # actual genomic start coord
        self.end      = end            # actual genomic end coord
        self.boundary = boundary       # right-hand boundary x-coord in phony space
        self.annot    = False          # block contains annotation exons?


if __name__ == "__main__":
    main()
