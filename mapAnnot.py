#!/usr/bin/env python

# Read annotation file, print counts, plot transcript sizes;

import os
import sys
import optparse
import re        # regular expressions

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import Annotations as anno

DEF_OUTPUT = 'transcript_sizes.png'
DEF_XMAX   = 4000
DEF_TITLE = 'Transcript lengths in annotation file'

SIZE_BINS = (1, 3, 5, 8, 9999999)
SIZE_COLORS = ('orange', 'red', 'yellow', 'green', 'blue')
SIZE_LEGENDS = ('1 exon', '2/3 exons', '4/5 exons', '6/7/8 exons', '8 exons')

def main ():

    opt, args = getParms()

    gtf = args[0]

    annotList = anno.AnnotationList (gtf)

    tranSizes = list()
    for ix in xrange(len(SIZE_BINS)):
        tranSizes.append (list())

    for chr in annotList.chromosomes():
        for strand in annotList.strands(chr):
            for geneEnt in annotList.geneList(chr, strand).nextChild():
                for tranEnt in geneEnt.nextChild():

                    if tranEnt.length > 100000:
                        print '%-5s  %s  %2d  %5d  %s' % (chr, strand, tranEnt.numChildren(), tranEnt.length, tranEnt.name)
                
                    for ix, size in enumerate(SIZE_BINS):
                        if tranEnt.numChildren() <= size:
                            tranSizes[ix].append (tranEnt.length)
                            break

    plt.figure (figsize=(12, 6))
    counts, bins, patches = plt.hist(tranSizes, bins=80, range=(0,opt.xmax), rwidth=0.8, color=SIZE_COLORS, histtype='barstacked', label=SIZE_LEGENDS)
    plt.legend(loc='best', prop={'size':10})
    plt.xlabel('transcript length')
    plt.ylabel('number of transcripts')

    if opt.title is not None:
        plt.suptitle(opt.title)

    plt.savefig (opt.output)
    plt.close()

    print counts
    print bins
    print patches


def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options] <fasta_file> ... ')

    parser.add_option ('--output',               help='Output file name (def: %default)')
    parser.add_option ('--xmax',     type='int', help='Maximum transcript length (x-axis) to plot (def: %default)')
    parser.add_option ('--title',                help='Title for top of figure (def: %default)')

    parser.set_defaults (output=DEF_OUTPUT,
                         xmax=DEF_XMAX,
                         title=DEF_TITLE,
                         )

    opt, args = parser.parse_args()

    return opt, args



if __name__ == "__main__":
    main()
