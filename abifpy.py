#!/usr/bin/env python

"""Python module for reading .ab1 trace files"""

import datetime
import re
import struct
import sys

# dictionary for deciding which values to extract and contain in self.meta
TAGS = {
            'DATA1':'raw1',
            'DATA2':'raw2',
            'DATA3':'raw3',
            'DATA4':'raw4',
            'DySN1':'dyename',
            'GTyp1':'polymer',
            'FWO_1':'baseorder',
            'MODL1':'model', 
            'SMPL1':'samplename', 
            'TUBE1':'well',
       } 

__version__ = '0.35'

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
                self.dir = {}
                self.meta['id'] = inFile.replace('.ab1','').replace('.fsa','')
                self.trimming = trimming
                # values contained in file header
                self._header = struct.unpack(fmtHead, self._raw[:30])
                # file format version
                self.version = self._header[1]

                # build dictionary of data tags and metadata
                for entry in self._parse():
                    key = entry.tagName + str(entry.tagNum)
                    self.dir[key] = entry
                    # only extract data from tags we care about
                    if key in TAGS:
                        # e.g. self.meta['well'] = 'B6'
                        self.meta[TAGS[key]] = self.getData(key)
    
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
            yield TraceDir(dirContent, self._raw)
            index += 1
   
    def getData(self, key):
        return self.dir[key].tagData

    def seq(self):
        """Returns sequence contained in the trace file."""
        data = self.getData('PBAS2')
        seq = re.sub("K|Y|W|M|R|S",'N',data)

        if self.trimming:
            return self.trim(seq)
        else:
            return seq

    def qual(self, char=True):
        """Returns a list of read quality values.
        
        Keyword argument:
        char -- True: returns ascii representation of phred values, False: returns phredvalues
        """
        data = self.getData('PCON2')
        qualList = []

        if not char:    
            for value in data:
                qualList.append(str(ord(value)))
        else:
            for value in data:
                if ord(value) > 93:
                    value = chr(93)
                qualList.append(chr(ord(value) + 33))
        # return value is list to make it compatible with the char option
        # and with the write() function (for writing .qual files)
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
                      description=self.meta['sampleid'])
        else:
            print 'Biopython was not detected. No SeqRecord was created.'
            return None

    def export(self, outFile="", qual=0):
        """Writes the trace file sequence to a fasta file.
        
        Keyword argument:
        outFile -- output file name (detault 'tracefile'.fa)
        qual -- 0: write fasta file, 1: write qual file, 2: write fastq file

        """
        
        if outFile == "":
            fileName = self.meta['id']
            if qual == 0:
                fileName += '.fa'
            elif qual == 1:
                fileName += '.qual'
            elif qual == 2:
                fileName += '.fq'
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
            qualHack = (ord(x) for x in self.getData('PCON2'))
            
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

class TraceDir(object):
    """Class representing directory content."""
    def __init__(self, tagEntry, rawData):
        self.tagName = tagEntry[0]
        self.tagNum = tagEntry[1]
        self.elemCode = tagEntry[2]
        self.elemSize = tagEntry[3]
        self.elemNum = tagEntry[4]
        self.dirSize = tagEntry[5]
        self.dirOffset = tagEntry[6]
        self.dirHandle = tagEntry[7]
        self.tagOffset = tagEntry[8]

        # if data size is <= 4 bytes, data is stored inside the directory
        # so offset needs to be changed
        if self.dirSize <= 4:
            self.dirOffset = self.tagOffset + 20

        if self.elemCode == 1:
            fmt = ">{0}b".format(self.elemNum)
            self.tagData = struct.unpack(fmt, 
                    rawData[self.dirOffset:self.dirOffset+self.dirSize])[0]
        elif self.elemCode == 2:
            fmt = ">{0}s".format(self.elemNum)
            self.tagData = struct.unpack(fmt, 
                    rawData[self.dirOffset:self.dirOffset+self.dirSize])[0]
        elif self.elemCode == 4:
            fmt = ">{0}h".format(self.elemNum)
            self.tagData = struct.unpack(fmt, 
                    rawData[self.dirOffset:self.dirOffset+self.dirSize])
        elif self.elemCode == 5:
            fmt = ">{0}l".format(self.elemNum)
            self.tagData = struct.unpack(fmt, 
                    rawData[self.dirOffset:self.dirOffset+self.dirSize])[0]
        elif self.elemCode == 7:
            fmt = ">{0}f".format(self.elemNum)
            self.tagData = struct.unpack(fmt,
                    rawData[self.dirOffset:self.dirOffset+self.dirSize])[0]
        elif self.elemCode == 10:
            fmt = ">hbb"
            year, month, date = struct.unpack(fmt,
                    rawData[self.dirOffset:self.dirOffset+self.dirSize])
            self.tagData = datetime.date(year, month, date)
        elif self.elemCode == 11:
            fmt = ">4b"
            hour, minute, second, hsecond = struct.unpack(fmt,
                    rawData[self.dirOffset:self.dirOffset+self.dirSize])
            self.tagData = datetime.time(hour, minute, second, hsecond)
        elif self.elemCode == 12:
            fmt = ">iibb"
            self.tagData = struct.unpack(fmt,
                    rawData[self.dirOffset:self.dirOffset+self.dirSize])[0]
        elif self.elemCode == 18:
            fmt = ">{0}s".format(self.elemNum-1)
            self.tagData = struct.unpack(fmt,
                    rawData[self.dirOffset+1:self.dirOffset+self.dirSize])[0].strip()
        elif self.elemCode == 19:
            fmt = ">{0}s".format(self.elemNum-1)
            self.tagData = struct.unpack(fmt,
                    rawData[self.dirOffset:self.dirOffset+self.dirSize-1])[0]
        elif self.elemCode == 1024:
            self.tagData = 'not available'
        else:
            self.tagData = None
    
