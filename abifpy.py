#!/usr/bin/env python

"""Python module for reading .ab1 trace files"""

import datetime
import re
import struct

# dictionary for deciding which values to extract and contain in self.meta
EXTRACT = {
            'DATA1':'raw1',
            'DATA2':'raw2',
            'DATA3':'raw3',
            'DATA4':'raw4',
            'DySN1':'dyename',
            'GTyp1':'polymer',
            'FWO_1':'baseorder',
            'MODL1':'model', 
            'PLOC2':'tracepeaks',
            'SMPL1':'sample', 
            'TUBE1':'well',
          }     

_BYTEFMT = {
            1: ">{0}b", 2: ">{0}s", 4: ">{0}h", 5: ">{0}i", 
            7: ">{0}f", 10: ">h2b", 11: ">4b", 12: ">2i2b",
            18: ">{0}s", 19: ">{0}s"
           }

__version__ = '0.3.6'

class Trace(object):
    """Class representing trace file"""
    def __init__(self, inFile, trimming=False):        
        with open(inFile) as source:
            self._raw = source.read()
            try:
                if not self._raw[:4] == 'ABIF':
                    raise IOError('Input is not a valid trace file')
            except IOError:
                source = None
                raise
            else:
                # header data structure:
                # file type, file, version, tag name, tag number, element type code,
                # element size, number of elements, data size, data offset, handle,
                # file type, file version
                fmtHead = '>4sH4sIHHIII'
                # dictionary for containing file metadata
                self.meta = {}
                # dictionary for containing extracted directory data
                self.tags = {}
                self.meta['id'] = inFile.replace('.ab1','').replace('.fsa','')
                self.trimming = trimming
                # values contained in file header
                self._header = struct.unpack(fmtHead, self._raw[:30])
                # file format version
                self.version = self._header[1]

                # build dictionary of data tags and metadata
                for entry in self._parse():
                    key = entry.tagName + str(entry.tagNum)
                    self.tags[key] = entry
                    # only extract data from tags we care about
                    if key in EXTRACT:
                        # e.g. self.meta['well'] = 'B6'
                        self.meta[EXTRACT[key]] = self.get_data(key)
    
    def _parse(self):
        """Generator for directory contents."""
        # directory data structure:
        # file type, file, version, tag name, tag number, 
        # element type code, element size, number of elements
        # data size, data offset, handle
        fmtHead = '>4sIHHIIII'
        headElemSize = self._header[5]
        headElemNum = self._header[6]
        headOffset = self._header[8]
        index = 0
        
        while index < headElemNum:
            start = headOffset + index * headElemSize
            finish = headOffset + (index + 1) * headElemSize
            # added directory offset to tuple
            # to handle directories with data size <= 4 bytes
            dirContent =  struct.unpack(fmtHead, self._raw[start:finish]) + (start,)
            yield _TraceDir(dirContent, self._raw)
            index += 1
   
    def get_data(self, key):
        return self.tags[key].tagData

    def show_tag(self, key):
        print 'tag name:', self.tags[key].tagName
        print 'tag number:', self.tags[key].tagNum
        print 'element code:', self.tags[key].elemCode
        print 'element size:', self.tags[key].elemSize
        print 'number of element:', self.tags[key].elemNum
        print 'data size:', self.tags[key].dataSize
        print 'data offset:', self.tags[key].dataOffset
        print 'data handle:', self.tags[key].dataHandle
        print 'tag offset:', self.tags[key].tagOffset

    def seq(self, ambig=False):
        """Returns sequence contained in the trace file."""
        data = self.get_data('PBAS2')

        if not ambig:
            seq = re.sub("K|Y|W|M|R|S", 'N', data)
        else:
            seq = data

        if self.trimming:
            return self.trim(seq)
        else:
            return seq

    def qual(self, char=True):
        """Returns a list of read quality values.
        
        Keyword argument:
        char -- True: returns ascii representation of phred values, False: returns phred values
        """
        data = self.get_data('PCON2')
        qualList = []

        if not char:    
            for value in data:
                qualList.append(ord(value))
        else:
            for value in data:
                if ord(value) > 93:
                    value = chr(93)
                qualList.append(chr(ord(value) + 33))
        # return value is list to make it compatible with the char option
        # and with the export() function (for writing .qual files)
        if self.trimming:
            return self.trim(qualList)
        else:
            return qualList

    def seqrecord(self):
        """Returns a SeqRecord object of the trace file."""
        try:
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            biopython = True
        except ImportError:
            biopython = False

        if biopython:
            seq = self.seq()
            return SeqRecord(Seq(seq), id=self.meta['id'], name="",
                             description=self.meta['sample'],
                             letter_annotations = {"phred_quality":self.qual(char=False)})
        else:
            print 'Biopython was not detected. No SeqRecord was created.'
            return None

    def export(self, outFile="", qual='fasta'):
        """Writes the trace file sequence to a fasta file.
        
        Keyword argument:
        outFile -- output file name (detault 'tracefile'.fa)
        qual -- 'fasta': write fasta file, 'qual': write qual file, 'fastq': write fastq file

        """
        
        if outFile == "":
            fileName = self.meta['id']
            if qual == 'fasta':
                fileName += '.fa'
            elif qual == 'qual':
                fileName += '.qual'
            elif qual == 'fastq':
                fileName += '.fq'
            else:
                raise ValueError('Invalid file format: {0}.'.format(qual))
        else:
            fileName = outFile
        
        if qual == 0:
            contents = '>{0} {1}\n{2}\n'.format(
                        fileName, self.meta['sampleid'], self.seq())
        elif qual == 1:
            contents = '>{0} {1}\n{2}\n'.format(
                        fileName, self.meta['sampleid'], ' '.join(self.qual(char=False)))
        elif qual == 2:
            contents = '@{0} {1}\n{2}\n+{0} {1}\n{3}\n'.format(
                        fileName, self.meta['sampleid'], self.seq(), ''.join(self.qual()))

        with open(fileName, 'w') as outFile:
            outFile.writelines(contents)

    def trim(self, seq, cutoff=0.05):
        """Trims the sequence using Richard Mott's modified trimming algorithm.
        
        Keyword argument:
        seq -- sequence to be trimmed
        cutoff -- probability cutoff value

        Trimmed bases are determined from their segment score, ultimately
        determined from each base's quality values. If a segment score is
        negative, the entire segment will be trimmed.
        
        """
        # set flag for trimming
        take = False
        # set segment size for calculating segment score
        segment = 20
        trimStart = 0
        trimFinish = len(seq)
        
        if len(seq) <= segment:
            raise ValueError('Sequence can not be trimmed because it is shorter than the trim segment size')
        else:
            # calculate score for trimming
            scoreList =  []
            # actually the same as seq.qual(), but done to avoid infinite recursion
            qualHack = (ord(x) for x in self.get_data('PCON2'))
            
            for qual in qualHack:
                # calculate probability back from formula used
                # to calculate phred qual values
                prob = 10 ** (qual/-10.0)
                scoreList.append(cutoff - prob)
            # algorithm for obtaining trimming position values
            for index in xrange(len(seq)-segment):
                # calculate segment score
                # if segment score is negative, trim the region out
                segmentScore = reduce(lambda x,y:x+y,
                                 scoreList[index:index+segment])
                if not take and segmentScore > 0:
                    trimStart = index
                    take = True
                elif take and segmentScore <= 0:
                    # if segment length is longer than remaining bases
                    # take up everything until the end
                    if index + segment <= len(seq):
                        trimFinish = index + segment - 1
                        break
                    else:
                        break
            return seq[trimStart:trimFinish]

class _TraceDir(object):
    """Class representing directory content."""
    def __init__(self, tagEntry, rawData):
        self.tagName = tagEntry[0]
        self.tagNum = tagEntry[1]
        self.elemCode = tagEntry[2]
        self.elemSize = tagEntry[3]
        self.elemNum = tagEntry[4]
        self.dataSize = tagEntry[5]
        self.dataOffset = tagEntry[6]
        self.dataHandle = tagEntry[7]
        self.tagOffset = tagEntry[8]

        # if data size is <= 4 bytes, data is stored inside the directory
        # so offset needs to be changed
        if self.dataSize <= 4:
            self.dataOffset = self.tagOffset + 20

        self.tagData = self._unpack(rawData)

        # no need to use tuple if len == 1
        if self.elemCode not in [10, 11, 12, 1024] and len(self.tagData) == 1:
            self.tagData = self.tagData[0]

        # account for different data types
        if self.elemCode == 10:
            self.tagData = datetime.date(*(self.tagData))
        elif self.elemCode == 11:
            self.tagData = datetime.time(*(self.tagData))
        elif self.elemCode == 18:
            self.tagData = self.tagData[1:]
        elif self.elemCode == 19:
            self.tagData = self.tagData[:-1]
            
    def _unpack(self, rawData):
        start = self.dataOffset
        finish = self.dataOffset + self.dataSize
        num = self.elemNum
        data = rawData[start:finish]

        if self.elemCode in _BYTEFMT:
            return struct.unpack(_BYTEFMT[self.elemCode].format(num), data)
        else:
            return None
