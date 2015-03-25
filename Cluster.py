#!/usr/bin/env python

# Classes for managing IsoSeq clusters, including annotation matches.

import os
import sys
import re                       # for regular expressions
import cPickle as pickle

from tt_log import logger

VERSION = '20150109.01'
logger.debug('version %s loaded' % VERSION)

regexFP = re.compile('f(\d+)p(\d+)')      # finds full and partial read counts in cluster ID

class Cluster (object):

    def __init__ (self, name, flags, chr, start, strand, cigar, bases):

        self.name   = name
        self.flags  = flags
        self.chr    = chr
        self.start  = start
        self.strand = strand
        self.cigar  = cigar
        self.bases  = bases
        self.pctGC  = None     # computed and cached below
        self.bestGene  = None
        self.bestTran  = None
        self.bestScore = None

    def best (self, bestGene, bestTran, bestScore):
        '''Add match findings.'''
        
        self.bestGene  = bestGene
        self.bestTran  = bestTran
        self.bestScore = bestScore

        return

    def getFP (self):
        '''Return counts of full and partial reads making up the cluster.'''

        match = re.search (regexFP, self.name)
        if match is None:
            raise RuntimeError ('full/partial counts not found in cluster name %s' % self.name)
        return int(match.group(1)), int(match.group(2))

    def percentGC (self, window):
        '''Compute (and cache) %GC content across the read using sliding window of specified size.'''

        # Note that the self.pctGC list will be shorter than
        # self.bases by the window size. Also note that the %GC
        # reported for position x reflects the window of bases
        # starting at x. When plotting, you may want to adjust the
        # coordinates to center the window on x.

        if self.pctGC is None:

            self.pctGC = list()

            for ix in xrange(len(self.bases)-window):

                gc = self.bases.count ('G', ix, ix+window) \
                   + self.bases.count ('C', ix, ix+window)
                self.pctGC.append (float(gc) / float(window) * 100.0)

        return self.pctGC
        
class ClusterDict (object):

    def __init__ (self):

        self.clusterDict = dict()
        self.geneDict = None

    @staticmethod
    def fromPickle (filename):
        '''Create a ClusterDict object from a pickle file (alternative to __init__).'''

        handle = open (filename, 'r')
        pk = pickle.Unpickler (handle)
        clusterDict = pk.load()
        handle.close()

        logger.debug('read %d clusters in pickle format from %s' % (len(clusterDict), filename))

        return clusterDict

    def __len__ (self):
        return len(self.clusterDict)

    def addCluster (self, cluster):
        '''Add a cluster object to the dictionary.'''

        self.clusterDict[cluster.name] = cluster

        return

    def getGeneDict (self):
        '''Create and cache mapping from gene name to clusters which match it.'''

        if self.geneDict is None:
            self.geneDict = dict()
            for cluster in self.clusterDict.values():
                if cluster.bestGene is not None:
                    self.geneDict.setdefault(cluster.bestGene.name, []).append(cluster)
    
        return self.geneDict

    def getClustersForGene (self, gene):
        '''Generator function to return clusters for specified gene.'''

        gd = self.getGeneDict()
        if gene not in gd:
            return

        for cluster in gd[gene]:
            yield cluster

        return

    def toPickle (self, filename):

        self.geneDict =  None         # no need to pickle this, it can be recreated

        pickHandle = open (filename, 'w')
        pk = pickle.Pickler (pickHandle, pickle.HIGHEST_PROTOCOL)
        pk.dump (self)
        pickHandle.close()

        logger.debug('wrote %d clusters to pickle file %s' % (len(self.clusterDict), filename))

        return
