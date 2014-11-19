#!/usr/bin/env python

# Interface to an annotation .gtf file (e.g., GENCODE).

import os
import sys
import re                       # for regular expressions
import cPickle as pickle

from tt_log import logger

VERSION = '20141104.01'
logger.debug('version %s loaded' % VERSION)

class Annotation (object):

    def __init__ (self, start, end, strand, name):
        '''Annotation file data about a gene, transcript, or exon.'''

        self.start  = start
        self.end    = end
        self.strand = strand
        self.name   = name
        self.children = None
        
    def __getitem__ (self, ix):
        '''Support indexing for an Annotation object by returning children[ix].'''
        return self.children[ix]

    def __len__ (self):
        '''Support indexing for an Annotation object by returning len(children).'''
        return len(self.children)

    def updateStartEnd (self, start, end):
        '''
        With alternative format GTF annotation files, the start and
        end coordinates of genes and transcripts must be inferred from
        the exon entries. This method adjusts them given new information.
        '''

        if self.start > start:
            self.start = start
        if self.end   < end:
            self.end   = end

    def addChild (self, child):
        '''
        Add a child Annotation to the list of children of this
        Annotation object. Note that exons in the - strand are listed
        in descending position order in the annotation file. We want
        them always in ascending order, hence the start position check
        done here. It also constitutes a sanity check on the order of
        the incoming data.
        '''

        # The sort at the bottom of this block looks a little
        # scary. But the 500 sorts required to process the
        # out-of-sequence entries on chrY in the GENCODE19 gtf file
        # take only ~150 msec. It would be nice to sort only once, at
        # the end -- but the Annotation object doesn't know when the
        # end has arrived. It has to keep itself consistent after
        # every call to addChild.

        if self.children == None:
            self.children = [child]
        elif child.start >= self.children[-1].start:
            self.children.append(child)
        elif child.start <= self.children[0].start:
            self.children.insert(0, child)
        else:
####            raise RuntimeError('out of sequence: %s' % child.name)
####            logger.debug('out of sequence: %s' % child.name)
            self.children.append(child)
            self.children.sort(key=lambda x: x.start)
            
    def getChildren (self):
        '''Generator function to return children, which are Annotation objects themselves, one by one.'''

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
    
    def __init__ (self, filename, altFormat=False):
        '''
        self.annot is a dict keyed by chr.  Each dict entry is a
        top-level Annotation object for that chr. An Annotation object
        contains a list of next-level Annotation objects (its
        'children'). The children of a chr entry are genes.  The
        children of a gene entry are transcripts. The children of a
        transcript entry are exons.
        '''

        # When matchAnnot got wider use, it emerged that the format
        # for tags in GTF files is not standardized. The default path
        # here processes a GENCODE-like file, which includes separate
        # entries for genes, transcripts and exons in a hierarchical
        # arrangement. The altFormat path accepts a file where genes
        # and transcripts are identified implicitly by showing up in
        # exon entries.

        # I could have implemented this several ways. One logical
        # choice would have been to create two subclasses of the
        # AnnotationList class. But there are no obvious benefits to
        # that extra generality, since the two classes would differ
        # only in the constructor, and the caller's interface would be
        # a bit more clumsy. So all I've done is call one of two
        # methods out of the constructor, based on a flag passed by
        # the caller.

        # Note that there is a thrid way to create an AnnotationList
        # object: See AnnotationList.fromPickle below.

        self.filename = filename
        self.annot    = dict()       # this is the stuff! key=chr value=top-level Annotation object for chr
        self.geneDict = None         # lookup table by gene name (created and cached when needed)

        if not altFormat:
            self.initFromStandard()
        else:
            self.initFromAlt()

    def initFromStandard (self):

        logger.debug('reading annotations in standard format from %s' % self.filename)

        regexGene = re.compile ('gene_name \"([^\"]+)\"\;')
        regexTran = re.compile ('transcript_name \"([^\"]+)\"\;')
        regexTID  = re.compile ('transcript_id \"([^\"]+)\"\;')
        regexExon = re.compile ('exon_number \"?(\d+)\"?\;')         # some files have quotes around exon number, some don't

        numGenes = 0
        numTrans = 0
        numExons = 0

        handle = open (self.filename, 'r')

        for line in handle:

            if line.startswith('#'):               # skip comment line
                continue

            chr, source, type, start, end, score, strand, frame, attrs = line.strip().split('\t')
            start = int(start)
            end   = int(end)

#           chrEnt = self.annot.setdefault (chr, Annotation(0, 0, '+', chr))    # dummy top entry for chr
#               Used to do this as above. But that creates an Annotation object
#               every time, then throws it away if chr already exists!

            if chr not in self.annot:
                self.annot[chr] = Annotation(0, 0, '+', chr)       # dummy top entry for chr
            chrEnt = self.annot[chr]

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

        return

    def initFromAlt (self):
        '''Initialize an AnnotationList from a GTF annotation file in Alternative format.'''

        # Alternative format annotation files include entries for
        # exons only. There are no explicit entries for genes or
        # transcripts. Gene and transcript information must be
        # inferred from the exon entries. Specifically,
        # gene/transcript start and end coordinates will be the lower
        # and upper bounds of the exons they contain.

        logger.debug('reading annotations in alternate format from %s' % self.filename)

        regexGene = re.compile ('gene_name \"*([^\"]+)\"*\;')
        regexTran = re.compile ('transcript_name \"*([^\"]+)\"*\;')
        regexTID  = re.compile ('transcript_id \"*([^\"]+)\"*\;')
        regexExon = re.compile ('exon_number \"*(\d+)\"*\;')         # some files have quotes around exon number, some don't

        numGenes = 0
        numTrans = 0
        numExons = 0

        geneEnt = None
        tranEnt = None

        handle = open (self.filename, 'r')

        for line in handle:

            if line.startswith('#'):               # skip comment line
                continue

            chr, source, type, start, end, score, strand, frame, attrs = line.strip().split('\t')
            start = int(start)
            end   = int(end)

            if chr not in self.annot:
                self.annot[chr] = Annotation(0, 0, '+', chr)       # dummy top entry for chr
            chrEnt = self.annot[chr]

            if type == 'exon':

                match = re.search(regexGene, attrs)
                if match is None:
                    raise RuntimeError ('no gene_name field in %s' % line)
                geneName = match.group(1)

                if geneEnt is None or geneName != geneEnt.name:
                    geneEnt = Annotation (start, end, strand, geneName)
                    chrEnt.addChild (geneEnt)
                    numGenes += 1
                else:
                    geneEnt.updateStartEnd (start, end)       # expand start/end coords of current gene
                    
                # Capture transcript_name and transcript_id if they
                # both exist, otherwise take whichever we find.

                matchName = re.search(regexTran, attrs)
                matchID   = re.search(regexTID,  attrs)
                if matchID is not None:
                    tranID   = matchID.group(1)
                    tranName = matchName.group(1) if matchName is not None else tranID
                elif matchName is not None:
                    tranName = matchName.group(1)
                    tranID   = tranName
                else:
                    raise RuntimeError ('no transcript_name/transcript_id field in %s' % line)

                if tranEnt is None or tranName != tranEnt.name:
                    tranEnt = Annotation (start, end, strand, tranName)
                    tranEnt.ID     = tranID                   # only transcripts have ID and length attributes
                    tranEnt.length = 0
                    geneEnt.addChild(tranEnt)
                    numTrans += 1
                else:
                    tranEnt.updateStartEnd (start, end)       # expand start/end coords of current transcript

                tranEnt.length += end - start + 1             # add this exon to total transcript length
                    
                matchExon = re.search(regexExon, attrs)
                if matchExon is not None:
                    exonNum = int(matchExon.group(1))
                    if exonNum != tranEnt.numChildren() + 1:
                        raise RuntimeError ('transcript name %s exons out of sequence' % (tranName))
                else:
                    exonNum = tranEnt.numChildren() + 1

                exonName = '%s/%d' % (tranName, exonNum)     # exons don't have names: make one up

                exonEnt = Annotation (start, end, strand, exonName)
                tranEnt.addChild(exonEnt)
                numExons += 1

            elif type == 'start_codon':
                tranEnt.startcodon = start
            elif type == 'stop_codon':
                tranEnt.stopcodon  = start

        handle.close()

        logger.debug('read %d genes, %d transcripts, %d exons' % (numGenes, numTrans, numExons))

        return

    @staticmethod
    def fromPickle (filename):
        '''Create an AnnotationList object from a pickle file (alternative to __init__).'''

        logger.debug('reading annotations in pickle format from %s' % filename)

        handle = open (filename, 'r')
        pk = pickle.Unpickler (handle)
        annotList = pk.load()
        handle.close()

        return annotList

    def toPickle (self, filename):

        pickHandle = open (filename, 'w')
        pk = pickle.Pickler (pickHandle, pickle.HIGHEST_PROTOCOL)
        pk.dump (self)
        pickHandle.close()

        logger.debug('wrote annotation data to pickle file %s' % filename)

        return

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

    def getGeneDict (self):
        '''Create and cache dict: key=gene name, value=Annotation object for gene.'''

        if self.geneDict is None:

            logger.debug('creating gene name lookup table')

            self.geneDict = dict()

            for chr in self.chromosomes():
                for gene in self.annot[chr].getChildren():
                    self.geneDict.setdefault(gene.name, []).append(gene)

        return self.geneDict

    def getGene (self, geneName):

        gd = self.getGeneDict()
        if geneName in gd:
            return gd[geneName]
        else:
            return None


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
