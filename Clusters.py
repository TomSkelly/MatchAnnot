#!/usr/bin/env python

# Interface to cluster_report.csv file produced by IsoSeq.

import os
import sys
import re                       # for regular expressions
from tt_log import logger

class ClusterList (object):
    
    def __init__ (self, filename):
        '''
        '''

        logger.debug('reading clusters from %s' % filename)

        self.filename = filename
        self.clusters = dict()       # this is the stuff
        self.cells    = dict()       # key=cell long name  value=cell number
        self.numCells = 0
        self.numClusters = 0

        regexCCS = re.compile ('_CCS$')

        handle = open (filename, 'r')
        handle.readline()            # skip header line

        for line in handle:

            self.numClusters += 1

            clusterID, readName, FL = line.strip().split()
            cell, ZMW, coords = readName.split('/')
            coords = re.sub (regexCCS, '', coords)          # get rid of '_CCS' at end of read range
            shortName = ZMW + '|' + coords

            if cell not in self.cells:                      # have we seen this cell before?
                self.numCells += 1
                self.cells[cell] = self.numCells            # if not, give it a number
            cellNo = self.cells[cell]
            
            clusterEnt = self.clusters.setdefault(clusterID, {}).setdefault(FL, {}).setdefault(cellNo, [])
            clusterEnt.append (shortName)

        handle.close()

        logger.debug('read %d clusters from %d cells' % (self.numClusters, self.numCells))

    def showReads (self, clusterID):
        '''Generator function returns all reads for a spccified cluster.'''

        clusterEnt = self.clusters[clusterID]

        for FL in sorted(clusterEnt.keys()):
            for cellNo in sorted(clusterEnt[FL].keys()):
                yield (FL, cellNo, clusterEnt[FL][cellNo])      # last item is a list of reads

        return

    def showCells (self):
        '''Generator function returns long cell names and their indexes, in index order.'''

        for cell in sorted(self.cells, key=lambda cell: self.cells[cell]):
            yield self.cells[cell], cell

        return
