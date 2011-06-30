#!/usr/bin/env python

"""Python module for reading .ab1 trace files"""

import datetime
import re
import struct

# dictionary for deciding which values to extract and contain in self.meta
EXTRACT = {
            'SMPL1':'sample', 
            'TUBE1':'well',
            'DySN1':'dye',
            'GTyp1':'polymer',
            'MODL1':'model', 
            'DATA1':'raw1',
            'DATA2':'raw2',
            'DATA3':'raw3',
            'DATA4':'raw4',
            'PLOC2':'tracepeaks',
            'FWO_1':'baseorder',
          }     

# dictionary for unpacking tag values
_BYTEFMT = {
            1: 'b',     # byte
            2: 's',     # char
            3: 'H',     # word
            4: 'h',     # short
            5: 'i',     # long
            6: '2i',    # rational, legacy unsupported
            7: 'f',     # float
            8: 'd',     # double
            10: 'h2B',  # date
            11: '4B',   # time
            12: '2i2b', # thumb
            13: 'B',    # bool
            14: '2h',   # point, legacy unsupported
            15: '4h',   # rect, legacy unsupported
            16: '2i',   # vPoint, legacy unsupported
            17: '4i',   # vRect, legacy unsupported
            18: 's',    # pString
            19: 's',    # cString
            20: '2i',   # Tag, legacy unsupported
           }

__version__ = '0.4'

class Trace(object):
    """Class representing trace file"""
    def __init__(self, inFile, trimming=False):        
        with open(inFile, 'rb') as source:
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
            # dictionary for containing file metadata
            self.meta = {}
            # dictionary for containing extracted directory data
            self.tags = {}
            self.meta['id'] = inFile
            self.trimming = trimming
            # values contained in file header
            self._header = struct.unpack('>4sH4sI2H3I', self._raw[:30])
            # file format version
            self.version = self._header[1]

            # build dictionary of data tags and metadata
            for entry in self._parse_tag():
                key = entry.tagName + str(entry.tagNum)
                self.tags[key] = entry
                # only extract data from tags we care about
                if key in EXTRACT:
                    # e.g. self.meta['well'] = 'B6'
                    self.meta[EXTRACT[key]] = self.get_data(key)
    
    def __str__(self):
        """Prints data associated with the file."""
        summary = ['Trace file: {0}'.format(self.meta['id'])]
        summary.append('Sequence name: {0}'.format(self.meta['sample']))
        summary.append('Sequence:\n{0}'.format(self.seq()))
        summary.append('Quality values:\n{0}'.format(''.join(self.qual())))
        summary.append('Sample well: {0}'.format(self.meta['well']))
        summary.append('Dye: {0}'.format(self.meta['dye']))
        summary.append('Polymer: {0}'.format(self.meta['polymer']))

        return '\n'.join(summary)
    
    def _parse_tag(self):
        """Generator for directory contents."""
        # directory data structure:
        # file type, file, version, tag name, tag number, 
        # element type code, element size, number of elements
        # data size, data offset, handle
        headElemSize = self._header[5]
        headElemNum = self._header[6]
        headOffset = self._header[8]
        index = 0
        
        while index < headElemNum:
            start = headOffset + index * headElemSize
            finish = headOffset + (index + 1) * headElemSize
            # added directory offset to tuple
            # to handle directories with data size <= 4 bytes
            dirContent =  struct.unpack('>4sI2H4I', self._raw[start:finish]) + (start,)
            yield _TraceDir(dirContent, self._raw)
            index += 1

    def get_data(self, key):
        """Returns data stored in a tag."""
        return self.tags[key].tagData

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

        if not char:    
            qualList = [ord(value) for value in data]
        else:
            qualList = [chr(ord(value) + 33) for value in data]
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

    def export(self, outFile="", fmt='fasta'):
        """Writes the trace file sequence to a fasta file.
        
        Keyword argument:
        outFile -- output file name (detault 'tracefile'.fa)
        fmt -- 'fasta': write fasta file, 'qual': write qual file, 'fastq': write fastq file

        """
        if outFile == "":
            fileName = self.meta['id']
            if fmt == 'fasta':
                fileName += '.fa'
            elif fmt == 'qual':
                fileName += '.qual'
            elif fmt == 'fastq':
                fileName += '.fq'
            else:
                raise ValueError('Invalid file format: {0}.'.format(fmt))
        else:
            fileName = outFile
        
        if fmt == 'fasta':
            contents = '>{0} {1}\n{2}\n'.format(
                        self.meta['id'], 
                        self.meta['sample'], 
                        self.seq())
        elif fmt == 'qual':
            contents = '>{0} {1}\n{2}\n'.format(
                        self.meta['id'], 
                        self.meta['sample'], 
                        ' '.join(map(str, self.qual(char=False))))
        elif fmt == 'fastq':
            contents = '@{0} {1}\n{2}\n+{0} {1}\n{3}\n'.format(
                        self.meta['id'], 
                        self.meta['sample'], 
                        self.seq(), ''.join(self.qual()))

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
            # calculate probability back from formula used
            # to calculate phred qual values
            scoreList = [cutoff - (10 ** (qual/-10.0)) for qual in self.qual(char=False)]
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

    def __str__(self):
        """Prints data associated with a tag."""
        summary = ['Tag name: {0}'.format(self.tagName)]
        summary.append('Tag number: {0}'.format(self.tagNum))
        summary.append('Element code: {0}'.format(self.elemCode))
        summary.append('Element size: {0}'.format(self.elemSize))
        summary.append('Number of element(s): {0}'.format(self.elemNum))
        summary.append('Data size: {0}'.format(self.dataSize))
        summary.append('Data offset: {0}'.format(self.dataOffset))
        summary.append('Data handle: {0}'.format(self.dataHandle))
        summary.append('Tag offset: {0}'.format(self.tagOffset))
       
        return '\n'.join(summary)

    def _unpack(self, raw):
        """Returns tag data"""
        if self.elemCode in _BYTEFMT:
            
            # because ">1s" unpacks differently from ">s"
            num = '' if self.elemNum == 1 else str(self.elemNum)
            fmt = '>' + num +  _BYTEFMT[self.elemCode]
            start = self.dataOffset
            finish = self.dataOffset + self.dataSize

            data = struct.unpack(fmt, raw[start:finish])
            
            # no need to use tuple if len(data) == 1
            if self.elemCode not in [10, 11] and len(data) == 1:
                data = data[0]

            # account for different data types
            if self.elemCode == 10:
                return datetime.date(*data)
            elif self.elemCode == 11:
                return datetime.time(*data)
            elif self.elemCode == 13:
                return bool(data)
            elif self.elemCode == 18:
                return data[1:]
            elif self.elemCode == 19:
                return data[:-1]
            else:
                return data
        else:
            return None
