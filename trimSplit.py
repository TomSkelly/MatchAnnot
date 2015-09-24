#!/usr/bin/env python

# Given a fasta file of reads-of-insert, qsub N trimPrimers jobs, each
# of which will process a subset of the file. Combine the outpts when
# all are done.

import os
import sys
import optparse
import re        # regular expressions
import datetime

from tt_log import logger

VERSION = '20150403.01'

DEF_NJOBS = 16
DEF_TMPDIR = 'tmp'

STARTWAIT = 30                       # seconds to delay start of qsub'd jobs

def main ():

    logger.debug('version %s starting' % VERSION)

    opt, args = getParms()

    makeTempDir (opt.tmpdir)

    nSeqs = countSeqs (opt.input)
    logger.debug('%s contains %d sequences' % (opt.input, nSeqs))

    seqsPerJob = (nSeqs + opt.njobs - 1) / opt.njobs
    logger.debug('each of %d jobs will process %d sequences' % (opt.njobs, seqsPerJob))

    chunkList = makeFastaChunks (opt, nSeqs, seqsPerJob)

    for chunk in chunkList:
        chunk.makeScript()
        chunk.submitScript()

    submitFinalJobs (opt, chunkList)

    logger.debug('finished')

    return

def makeTempDir (dir):

    if os.path.isdir (dir):
        logger.warning('WARNING: temp directory %s already exists' % dir)
    else:
        os.makedirs (dir)
        
    return

def countSeqs (filename):
    '''Run grep -c ">" on a fasta file to count the sequences it contains.'''

    command = 'grep -c ">" %s' % filename
    popen_file = os.popen(command)
    response = popen_file.read().strip()
    rc       = popen_file.close()
    if rc is not None:
        logger.error("command failed, rc=%d" % rc)
        raise RuntimeError
    if not response.isdigit():
        logger.error("grep -c returned:" % response)
        raise RuntimeError
    
    return int(response)

def makeFastaChunks (opt, nSeqs, seqsPerJob):

    curPos = 1
    thisRec = 0
    chunkList = list()

    fastaIn = open (opt.input, 'r')
    line = fastaIn.readline()
    if not line.startswith('>'):
        raise RuntimeError ('first fasta line is not a header: %s' % line)

    while (line):                               # outer loop reads the whole input file

        endPos = min(curPos+seqsPerJob-1, nSeqs)
        fastaChunk = Chunk (opt, curPos, endPos)
        logger.debug('writing chunk %s' % fastaChunk.inputChunkName)
        with open (fastaChunk.inputChunkName, 'w') as chunkOut:

            while (line):                       # inner loop writes one chunk

                if line.startswith('>'):

                    thisRec += 1
                    if thisRec > seqsPerJob:
                        chunkList.append(fastaChunk)
                        thisRec = 0
                        break                   # break inner loop and exit with block, closing output file

                    curPos += 1

                chunkOut.write (line)
                line = fastaIn.readline()

    if thisRec > 0:
        chunkList.append(fastaChunk)

    fastaIn.close()

    return chunkList

def submitFinalJobs (opt, chunkList):

    chunkFiles = ['%s \\\n' % chk.trimmedChunkName for chk in chunkList]

    sh = list()
    sh.append('#!/bin/bash\n\n')
    sh.append('set -o errexit\n')
    sh.append('set -o nounset\n\n')

    sh.append('cat \\\n')
    sh.extend(chunkFiles)
    sh.append(' > %s\n' % opt.output)

    if opt.report is not None:
        reportFiles = ['%s \\\n' % chk.reportChunkName for chk in chunkList]
        sh.append('\ncat \\\n')
        sh.extend(reportFiles)
        sh.append(' > %s\n' % opt.report)

    finalScriptName = '%s/trim_final.sh' % opt.tmpdir
    handle =  open (finalScriptName, 'w')
    handle.writelines (sh)
    handle.close()

    deps = ':'.join ([chk.jobno for chk in chunkList])

    cmd = list()
    cmd.append('qsub')
    cmd.append('-N trim_final')       # job name
    cmd.append('-o trim_final.out')   # output file
    cmd.append('-j oe')               # combine stdout and stderr
    cmd.append('-l nodes=1:ppn=1,walltime=4:00:00')    # resources required
    cmd.append('-d . ')               # working directory (strangely, ./ is not the default)
    cmd.append('-r n')                # do NOT attempt to restart on failure
    cmd.append('-V')                  # export all environment variables to job
    cmd.append('-W umask=0002')       # make logs rw-rw-r--
    cmd.append('-m n')                # don't send any mail
    cmd.append('-W depend=afterok:%s' % deps)
    cmd.append(finalScriptName)       # script to run

    command = ' '.join(cmd)
    logger.debug ('running %s' % command)
    
    popen_file = os.popen(command)
    response = popen_file.read().strip()
    rc = popen_file.close()
    if rc is not None:
        logger.error('command failed, rc=%d' % rc)
        raise RuntimeError

    logger.debug ('jobno is %s' % response)

    return response

def getParms ():                       # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage='%prog [options]', version=VERSION)

    parser.add_option ('--input',   help='fasta file to be trimmed (required)')
    parser.add_option ('--primers', help='fasta file of primers (required)')
    parser.add_option ('--output',  help='output fasta file (required)')
    parser.add_option ('--report',  help='output alignments report text file (optional)')
    parser.add_option ('--njobs',   help='number of jobs to submit (def: %default)', type='int')
    parser.add_option ('--tmpdir',  help='temporary directory (def: %default)')

    parser.set_defaults (njobs=DEF_NJOBS,
                         tmpdir=DEF_TMPDIR,
                         )

    opt, args = parser.parse_args()

    return opt, args


class Chunk (object):
    '''Manage a chunk, including keeping all the file names consistent and in one place.'''

    JOBNO_PATTERN = re.compile('^(\d+)')  # the number part of job number

    def __init__ (self, opt, start, end):

        self.opt   = opt
        self.start = start
        self.end   = end
        self.jobno = None     # filled in when submitted
        self.jobName          = 'trim.%06d'                        % (start)
        self.scriptName       = '%s/trim.%06d_%06d.sh'             % (opt.tmpdir, start, end)
        self.scriptOutput     = '%s/trim.%06d_%06d.out'            % (opt.tmpdir, start, end)
        self.inputChunkName   = '%s/input_chunk.%06d_%06d.fasta'   % (opt.tmpdir, start, end)
        self.trimmedChunkName = '%s/trimmed_chunk.%06d_%06d.fasta' % (opt.tmpdir, start, end)
        if opt.report is None:
            self.reportChunkName  = None
        else:
            self.reportChunkName  = '%s/report_chunk.%06d_%06d.txt'  % (opt.tmpdir, start, end)

    def makeScript (self):

        sh = list()
        sh.append('#!/bin/bash\n\n')
        sh.append('set -o errexit\n')
        sh.append('set -o nounset\n\n')
        sh.append('~/work/MatchAnnot/trimPrimers.py \\\n')
        sh.append('    --input %s \\\n' % self.inputChunkName)
        sh.append('    --primers %s \\\n' % self.opt.primers)

        if self.opt.report is not None:
            sh.append('    --report %s \\\n' % self.reportChunkName)
            
        sh.append('    > %s\n' % self.trimmedChunkName)

        handle =  open (self.scriptName, 'w')
        handle.writelines (sh)
        handle.close()

        return

    def submitScript (self):

        # Dependent job submission will fail if parent has already
        # completed. So delay all job startups by a short amount of time.

        startAt = datetime.datetime.now() + datetime.timedelta(0, STARTWAIT)
        startAtStr = startAt.strftime('%Y%m%d%H%M.%S')

        cmd = list()
        cmd.append('qsub')
        cmd.append('-N %s' % self.jobName)           # job name
        cmd.append('-o %s' % self.scriptOutput)      # output file
        cmd.append('-j oe')               # combine stdout and stderr
        cmd.append('-l nodes=1:ppn=1,walltime=4:00:00')    # resources required
        cmd.append('-a %s' % startAtStr)  # delay start, see above
        cmd.append('-d . ')               # working directory (strangely, ./ is not the default)
        cmd.append('-r n')                # do NOT attempt to restart on failure
        cmd.append('-V')                  # export all environment variables to job
        cmd.append('-W umask=0002')       # make logs rw-rw-r--
        cmd.append('-m n')                # don't send any mail
        cmd.append(self.scriptName)       # script to run

        command = ' '.join(cmd)
        logger.debug ('running %s' % command)
        
        popen_file = os.popen(command)
        response = popen_file.read().strip()
        rc = popen_file.close()
        if rc is not None:
            logger.error('command failed, rc=%d' % rc)
            raise RuntimeError

        match = re.match (Chunk.JOBNO_PATTERN, response)
        if match is None:
            logger.error("invalid job sequence number: %s" % jobSeqStr)
            raise RuntimeError

        response = match.group(1)
        logger.debug ('jobno is %s' % response)
        self.jobno = response
        return response

if __name__ == "__main__":
    main()
