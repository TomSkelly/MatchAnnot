#!/usr/bin/env python

# Interface to an annotation .gtf file (e.g., GENCODE).

import os
import sys
import re                       # for regular expressions
from tt_log import logger

class Annotation (object):

    def __init__ (self, start, end, strand, name):
        '''
        Annotation file data about a gene, transcript, or exon.
        '''

        self.start  = start
        self.end    = end
        self.strand = strand
        self.name   = name
        self.children = None
        
    def addChild (self, child):
        '''
        Add a child Annotation to the list of children of this
        Annotation object. Not that exons in the - strand are listed
        in descending position order in the annotation file. We want
        them always in ascending order, hence the start position check
        done here. It also constitutes a sanity check on the order of
        the incoming data.
        '''

        if self.children == None:
            self.children = [child]
        elif child.start >= self.children[-1].start:
            self.children.append(child)
        elif child.start <= self.children[0].start:
            self.children.insert(0, child)
        else:
####            raise RuntimeError('out of sequence: %s' % child.name)
            logger.debug('out of sequence: %s' % child.name)
            self.children.append(child)
            self.children.sort()

    def __getitem__ (self, ix):
        '''Support indexing for an Annotation object by returning children[ix].'''
        return self.children[ix]

    def __len__ (self):
        '''Support indexing for an Annotation object by returning len(children).'''
        return len(self.children)

    def nextChild (self):
        '''Generator function to return children one by one.'''

        if self.children is None:
            return
        for child in self.children:
            yield child
        return

    def numChildren (self):
        if self.children is None:
            return 0
        return len(self.children)


class AnnotationList (object):
    
    def __init__ (self, filename):
        '''
        self.annot is a dict keyed by chr.  Each dict entry is a list
        of dicts, one per gene. Each such dict includes a list of
        children, which are dicts with the same structure as the
        parent. The children in a gene entry are transcripts. The
        children in a transcript entry are exons.
        '''

        logger.debug('reading annotations from %s' % filename)

        self.filename = filename
        self.annot    = dict()       # this is the stuff

        regexGene = re.compile ('gene_name \"([^\"]+)\"\; ')
        regexTran = re.compile ('transcript_name \"([^\"]+)\"\; ')
        regexTID  = re.compile ('transcript_id \"([^\"]+)\"\; ')
        regexExon = re.compile ('exon_number \"?(\d+)\"?\; ')         # some files have quotes around exon number, some don't

        numGenes = 0
        numTrans = 0
        numExons = 0

        handle = open (filename, 'r')

        for line in handle:

            if line.startswith('#'):               # skip comment line
                continue

            chr, source, type, start, end, score, strand, frame, attrs = line.strip().split('\t')
            start = int(start)
            end   = int(end)

            chrEnt = self.annot.setdefault (chr, Annotation(0, 0, '+', chr))    # dummy top entry for chr

            if type == 'gene':

                numGenes += 1
                geneName = re.search(regexGene, attrs).group(1)
                geneEnt = Annotation (start, end, strand, geneName)
                chrEnt.addChild(geneEnt)

            elif type == 'transcript':

                numTrans += 1
                geneName = re.search(regexGene, attrs).group(1)
                tranName = re.search(regexTran, attrs).group(1)
                tranID   = re.search(regexTID,  attrs).group(1)       # transcript id looks like: ENST00000456328.2

                if geneName != geneEnt.name:
                    raise RuntimeError ('gene name %s != %s in transcript %s' % (geneName, geneEnt.name, tranName))

                tranEnt = Annotation (start, end, strand, tranName)
                tranEnt.ID     = tranID                               # only transcripts have ID and length attributes
                tranEnt.length = 0
                geneEnt.addChild(tranEnt)

            elif type == 'exon':

                numExons += 1
                geneName = re.search(regexGene, attrs).group(1)
                tranName = re.search(regexTran, attrs).group(1)
                exonNum  = int(re.search(regexExon, attrs).group(1))
                exonName = '%s/%d' % (tranName, exonNum)     # exons don't have names: make one up

                if geneName != geneEnt.name:
                    raise RuntimeError ('gene name %s != %s in transcript %s' % (geneName, geneEnt.name, tranName))
                if tranName != tranEnt.name:
                    raise RuntimeError ('transcript name %s != %s in exon %d' % (tranName, tranEnt.name, exonNum))
                if exonNum  != tranEnt.numChildren() + 1:
                    raise RuntimeError ('transcript name %s exons out of sequence' % (tranName))

                tranEnt.length += end - start + 1            # add this exon to total transcript length

                exonEnt = Annotation (start, end, strand, exonName)
                tranEnt.addChild(exonEnt)

            elif type == 'start_codon':
                tranEnt.startcodon = start
            elif type == 'stop_codon':
                tranEnt.stopcodon  = start

        handle.close()

        logger.debug('read %d genes, %d transcripts, %d exons' % (numGenes, numTrans, numExons))

    def chromosomes (self):
        '''Generator function for the chromosome names in an AnnotationList object.'''

        for chr in sorted(self.annot.keys()):
            yield chr

        return

    def geneList (self, chr):
        '''Return a top-level Annotation object containing the genes for a chr.'''

        if chr not in self.annot:
            raise RuntimeError ('chromosome %s not in annotation' % chr)
        return self.annot[chr]


class AnnotationCursor (object):
    '''
    Class supporting the walking of an AnnotationList.
    '''

    def __init__ (self, annotList):

        self.annotList = annotList
        self.posit = dict( [ (chr, 0) for chr in annotList.chromosomes() ] )

    def advance (self, chr, start):
        '''
        Move the cursor past genes which end prior to the specified
        start position on specified chr.
        '''

        geneList = self.annotList.geneList (chr)
        numGenes = geneList.numChildren()
        curPos = self.posit[chr]
        while curPos < numGenes and geneList[curPos].end < start:
            curPos += 1
        self.posit[chr] = curPos

    def getOverlappingGenes (self, chr, start, end, strand):
        '''
        Generator function returns genes from the specified strand,
        starting from current cursor position, up to point where gene
        starts after end position of specified range. Does NOT change
        the current cursor position, since the next query may want
        some of the same genes. (To move the cursor, call advance.)
        '''

        # Here we need to deal with the case where gene A completely overlaps gene B:

        #         10000                                   20000
        #  gene A |-------------------------------------------|
        #  gene B      |----------------|
        #  read 1                           |-----------|
        #  read 2                                                  |-----------|

        # $pos will stick at 10000 until read 2 is encountered. But we
        # don't want to report gene B for read 1.


        geneList = self.annotList.geneList (chr)
        numGenes = geneList.numChildren()
        curPos = self.posit[chr]

        while curPos < numGenes and geneList[curPos].start <= end:
            
            curGene = geneList[curPos]

            if curGene.strand == strand and curGene.end > start:
                yield curGene

            curPos += 1

        return
