#!/usr/bin/env python

"""Python module for reading .ab1 trace files."""

import datetime
import struct

# dictionary for deciding which values to extract and contain in self.data
EXTRACT = {
            'TUBE1': 'well',
            'DySN1': 'dye',
            'GTyp1': 'polymer',
            'MODL1': 'model', 
            'RUND1': 'run start date',
            'RUND2': 'run finish date',
            'RUND3': 'data collection start date',
            'RUND4': 'data collection finish date',
            'RUNT1': 'run start time',
            'RUNT2': 'run finish time',
            'RUNT3': 'data collection start time',
            'RUNT4': 'data collection finish time',
            'DATA1': 'raw1',
            'DATA2': 'raw2',
            'DATA3': 'raw3',
            'DATA4': 'raw4',
            'PLOC2': 'tracepeaks',
            'FWO_1': 'baseorder',
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

# header structure
_HEADFMT = '>4sH4sI2H3I'

# directory data structure
_DIRFMT = '>4sI2H4I'

__version__ = '0.5'


class Trace(object):
    """Class representing trace file."""
    def __init__(self, inFile, trimming=False):        
        self._handle = open(inFile, 'rb')
        try:
            self._handle.seek(0)
            if not self._handle.read(4) == 'ABIF':
                raise IOError('Input is not a valid trace file')
        except IOError:
            self._handle = None
            raise
        else:
            # header data structure:
            # file type, file, version, tag name, tag number, element type code,
            # element size, number of elements, data size, data offset, handle,
            # file type, file version
            # dictionary for containing file metadata
            self.data = {}
            # dictionary for containing extracted directory data
            self.tags = {}
            self.trimming = trimming
            # values contained in file header
            self._handle.seek(0)
            header = struct.unpack(_HEADFMT, 
                     self._handle.read(struct.calcsize(_HEADFMT)))
            # file format version
            self.version = header[1]

            # build dictionary of data tags and metadata
            for entry in self._parse_header(header):
                key = entry.tagName + str(entry.tagNum)
                self.tags[key] = entry
                # only extract data from tags we care about
                if key in EXTRACT:
                    # e.g. self.data['well'] = 'B6'
                    self.data[EXTRACT[key]] = self.get_data(key)

            self.id = inFile.replace('.ab1', '')
            self.name = self.get_data('SMPL1')
            self.seq = self.get_data('PBAS2')
            self.qualVal = [ord(value) for value in self.get_data('PCON2')]
            self.qual = ''.join([chr(value + 33) for value in self.qualVal])

            if trimming:
                self.seq, self.qual, self.qualVal = map(self.trim, 
                                                        [self.seq, self.qual,
                                                        self.qualVal])

            self._handle.close()
    
    def __repr__(self):
        """Represents data associated with the file."""
        if len(self.seq) > 10:
            seq = "{0}...{1}".format(self.seq[:5], self.seq[-5:])
            qualVal = "[{0}, ..., {1}]".format(
                      repr(self.qualVal[:5])[1:-1], 
                      repr(self.qualVal[-5:])[1:-1])
        else:
            seq = self.seq
            qualVal = self.qualVal

        return "{0}({1}, qualVal:{2}, id:{3}, name:{4})".format(
                self.__class__.__name__, repr(seq), qualVal,
                repr(self.id), repr(self.name))
    
    def _parse_header(self, header):
        """Generator for directory contents."""
        # header structure:
        # file signature, file version, tag name, tag number, 
        # element type code, element size, number of elements
        # data size, data offset, handle
        headElemSize = header[5]
        headElemNum = header[6]
        headOffset = header[8]
        index = 0
        
        while index < headElemNum:
            start = headOffset + index * headElemSize
            # added directory offset to tuple
            # to handle directories with data size <= 4 bytes
            self._handle.seek(start)
            dirEntry =  struct.unpack(_DIRFMT, 
                        self._handle.read(struct.calcsize(_DIRFMT))) + (start,)
            index += 1
            yield _TraceDir(dirEntry, self._handle)

    def get_data(self, key):
        """Returns data stored in a tag."""
        return self.tags[key].tagData

    def seq_remove_ambig(self, seq):
        """Replaces extra ambiguous bases with 'N'."""
        seq = self.seq

        import re
        return re.sub("K|Y|W|M|R|S", 'N', seq)

    def export(self, outFile="", fmt='fasta'):
        """Writes the trace file sequence to a fasta file.
        
        Keyword argument:
        outFile -- output file name (detault 'tracefile'.fa)
        fmt -- 'fasta': write fasta file, 'qual': write qual file, 'fastq': write fastq file

        """
        if outFile == "":
            fileName = self.id
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
                        self.id, 
                        self.name, 
                        self.seq)
        elif fmt == 'qual':
            contents = '>{0} {1}\n{2}\n'.format(
                        self.id, 
                        self.name, 
                        ' '.join(map(str, self.qualVal)))
        elif fmt == 'fastq':
            contents = '@{0} {1}\n{2}\n+{0} {1}\n{3}\n'.format(
                        self.id, 
                        self.name, 
                        self.seq, ''.join(self.qual))

        with open(fileName, 'w') as outFile:
            outFile.writelines(contents)

    def trim(self, seq, cutoff=0.05):
        """Trims the sequence using Richard Mott's modified trimming algorithm.
        
        Keyword argument:
        seq -- sequence to be trimmed
        cutoff -- probability cutoff value

        Trimmed bases are determined from their segment score, ultimately
        determined from each base's quality values. 
        
        More on:
        http://www.phrap.org/phredphrap/phred.html
        http://www.clcbio.com/manual/genomics/Quality_trimming.html
        """
        # set flag for trimming
        start = False
        # set minimum segment size
        segment = 20
        trimStart = 0
        
        if len(seq) <= segment:
            raise ValueError('Sequence can not be trimmed because \
                             it is shorter than the trim segment size')
        else:
            # calculate probability back from formula used
            # to calculate phred qual values
            scoreList = [cutoff - (10 ** (qual/-10.0)) for 
                         qual in self.qualVal]

            # calculate cummulative scoreList
            # if cummulative value < 0, set to 0
            # first value is set to 0 (assumption: trimStart is always > 0)
            runningSum = [0]
            for i in xrange(1, len(scoreList)):
                num = runningSum[-1] + scoreList[i]
                if num < 0:
                    runningSum.append(0)
                else:
                    runningSum.append(num)
                    if not start:
                        # trimStart = value when cummulative starts to be > 0
                        trimStart = i
                        start = True

            # trimFinish = index of the highest cummulative value,
            # marking the segment with the highest cummulative score 
            trimFinish = runningSum.index(max(runningSum)) 

            return seq[trimStart:trimFinish]

class _TraceDir(object):
    """Class representing directory content."""
    def __init__(self, tagEntry, handle):
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

        self.tagData = self._unpack(handle)

    def __repr__(self):
        """Represents data associated with a tag."""
        summary = ['tagName: {0}'.format(repr(self.tagName))]
        summary.append('tagNumber: {0}'.format(repr(self.tagNum)))
        summary.append('elemCode: {0}'.format(repr(self.elemCode)))
        summary.append('elemSize: {0}'.format(repr(self.elemSize)))
        summary.append('elemNum: {0}'.format(repr(self.elemNum)))
        summary.append('dataSize: {0}'.format(repr(self.dataSize)))
        summary.append('dataOffset: {0}'.format(repr(self.dataOffset)))
        summary.append('dataHandle: {0}'.format(repr(self.dataHandle)))
        summary.append('tagOffset: {0}'.format(repr(self.tagOffset)))
        summary.append('tagData: {0}'.format(repr(self.tagData)))
       
        return '\n'.join(summary)

    def _unpack(self, handle):
        """Returns tag data"""
        if self.elemCode in _BYTEFMT:
            
            # because ">1s" unpacks differently from ">s"
            num = '' if self.elemNum == 1 else str(self.elemNum)
            fmt = "{0}{1}{2}".format('>', num, _BYTEFMT[self.elemCode])
            start = self.dataOffset
    
            handle.seek(start)
            data = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))
            
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
