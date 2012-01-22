#!/usr/bin/env python
#
# abifpy.py
# python module for reading abi trace files
# http://github.com/bow/abifpy

"""Python module for reading .ab1 trace files."""

import datetime
import struct
from os.path import splitext, basename

from sys import version_info

__all__ = ['Trace']

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

__version__ = '0.9'


# to handle py3 IO
def py3_get_string(byte):
    if version_info[0] < 3:
        return byte
    else:
        return byte.decode()

def py3_get_byte(string):
    if version_info[0] < 3:
        return string
    else:
        return string.encode()

class Trace(object):
    """Class representing trace file."""
    def __init__(self, in_file, trimming=False):        
        self._handle = open(in_file, 'rb')
        try:
            self._handle.seek(0)
            if not self._handle.read(4) == py3_get_byte('ABIF'):
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
                key = entry.tag_name + str(entry.tag_num)
                self.tags[key] = entry
                # only extract data from tags we care about
                if key in EXTRACT:
                    # e.g. self.data['well'] = 'B6'
                    self.data[EXTRACT[key]] = self.get_data(key)

            self.id = self._get_file_id(in_file)
            self.name = self.get_data('SMPL1')
            self.seq = self.get_data('PBAS2')
            self.qual = ''.join([chr(ord(value) + 33) for value in self.get_data('PCON2')])
            self.qual_val = [ord(value) for value in self.get_data('PCON2')]

            if trimming:
                self.seq, self.qual, self.qual_val = map(self.trim, 
                                                        [self.seq, self.qual,
                                                        self.qual_val])

    def __repr__(self):
        """Represents data associated with the file."""
        if len(self.seq) > 10:
            seq = "{0}...{1}".format(self.seq[:5], self.seq[-5:])
            qual_val = "[{0}, ..., {1}]".format(
                      repr(self.qual_val[:5])[1:-1], 
                      repr(self.qual_val[-5:])[1:-1])
        else:
            seq = self.seq
            qual_val = self.qual_val

        return "{0}({1}, qual_val:{2}, id:{3}, name:{4})".format(
                self.__class__.__name__, repr(seq), qual_val,
                repr(self.id), repr(self.name))
    
    def _parse_header(self, header):
        """Generator for directory contents."""
        # header structure:
        # file signature, file version, tag name, tag number, 
        # element type code, element size, number of elements
        # data size, data offset, handle
        head_elem_size = header[5]
        head_elem_num = header[6]
        head_offset = header[8]
        index = 0
        
        while index < head_elem_num:
            start = head_offset + index * head_elem_size
            # added directory offset to tuple
            # to handle directories with data size <= 4 bytes
            self._handle.seek(start)
            dir_entry =  struct.unpack(_DIRFMT, 
                        self._handle.read(struct.calcsize(_DIRFMT))) + (start,)
            index += 1
            yield _TraceDir(dir_entry, self._handle)

    def _get_file_id(self, in_file):
        """Returns filename without extension."""
        return splitext(basename(in_file))[0]

    def close(sel):
        """Closes the Trace file object."""
        self._handle.close()
    

    def get_data(self, key):
        """Returns data stored in a tag."""
        return self.tags[key].tag_data

    def seq_remove_ambig(self, seq):
        """Replaces extra ambiguous bases with 'N'."""
        import re
        seq = self.seq
        return re.sub("K|Y|W|M|R|S", 'N', seq)

    def export(self, out_file="", fmt='fasta'):
        """Writes the trace file sequence to a fasta file.
        
        Keyword argument:
        out_file -- output file name (detault 'tracefile'.fa)
        fmt -- 'fasta': write fasta file, 'qual': write qual file, 'fastq': write fastq file

        """
        if out_file == "":
            file_name = self.id
            if fmt == 'fasta':
                file_name += '.fa'
            elif fmt == 'qual':
                file_name += '.qual'
            elif fmt == 'fastq':
                file_name += '.fq'
            else:
                raise ValueError('Invalid file format: {0}.'.format(fmt))
        else:
            file_name = out_file
        
        if fmt == 'fasta':
            contents = '>{0} {1}\n{2}\n'.format(
                        self.id, 
                        self.name, 
                        self.seq)
        elif fmt == 'qual':
            contents = '>{0} {1}\n{2}\n'.format(
                        self.id, 
                        self.name, 
                        ' '.join(map(str, self.qual_val)))
        elif fmt == 'fastq':
            contents = '@{0} {1}\n{2}\n+{0} {1}\n{3}\n'.format(
                        self.id, 
                        self.name, 
                        self.seq, ''.join(self.qual))

        with open(file_name, 'w') as out_file:
            out_file.writelines(contents)

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
        trim_start = 0
        
        if len(seq) <= segment:
            raise ValueError('Sequence can not be trimmed because \
                             it is shorter than the trim segment size')
        else:
            # calculate probability back from formula used
            # to calculate phred qual values
            score_list = [cutoff - (10 ** (qual/-10.0)) for 
                         qual in self.qual_val]

            # calculate cummulative score_list
            # if cummulative value < 0, set to 0
            # first value is set to 0 (assumption: trim_start is always > 0)
            running_sum = [0]
            for i in range(1, len(score_list)):
                num = running_sum[-1] + score_list[i]
                if num < 0:
                    running_sum.append(0)
                else:
                    running_sum.append(num)
                    if not start:
                        # trim_start = value when cummulative starts to be > 0
                        trim_start = i
                        start = True

            # trim_finish = index of the highest cummulative value,
            # marking the segment with the highest cummulative score 
            trim_finish = running_sum.index(max(running_sum)) 

            return seq[trim_start:trim_finish]

class _TraceDir(object):
    """Class representing directory content."""
    def __init__(self, tag_entry, handle):
        self.tag_name = py3_get_string(tag_entry[0])
        self.tag_num = tag_entry[1]
        self.elem_code = tag_entry[2]
        self.elem_size = tag_entry[3]
        self.elem_num = tag_entry[4]
        self.data_size = tag_entry[5]
        self.data_offset = tag_entry[6]
        self.data_handle = tag_entry[7]
        self.tag_offset = tag_entry[8]

        # if data size is <= 4 bytes, data is stored inside the directory
        # so offset needs to be changed
        if self.data_size <= 4:
            self.data_offset = self.tag_offset + 20

        self.tag_data = self._unpack(handle)

    def __repr__(self):
        """Represents data associated with a tag."""
        summary = ['tag_name: {0}'.format(repr(self.tag_name))]
        summary.append('tag_number: {0}'.format(repr(self.tag_num)))
        summary.append('elem_code: {0}'.format(repr(self.elem_code)))
        summary.append('elem_size: {0}'.format(repr(self.elem_size)))
        summary.append('elem_num: {0}'.format(repr(self.elem_num)))
        summary.append('data_size: {0}'.format(repr(self.data_size)))
        summary.append('data_offset: {0}'.format(repr(self.data_offset)))
        summary.append('data_handle: {0}'.format(repr(self.data_handle)))
        summary.append('tag_offset: {0}'.format(repr(self.tag_offset)))
        summary.append('tag_data: {0}'.format(repr(self.tag_data)))
       
        return '\n'.join(summary)

    def _unpack(self, handle):
        """Returns tag data"""
        if self.elem_code in _BYTEFMT:
            
            # because ">1s" unpacks differently from ">s"
            num = '' if self.elem_num == 1 else str(self.elem_num)
            fmt = "{0}{1}{2}".format('>', num, _BYTEFMT[self.elem_code])
            start = self.data_offset
    
            handle.seek(start)
            data = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))
            
            # no need to use tuple if len(data) == 1
            if self.elem_code not in [10, 11] and len(data) == 1:
                data = data[0]

            # account for different data types
            if self.elem_code == 2:
                return py3_get_string(data)
            elif self.elem_code == 10:
                return datetime.date(*data)
            elif self.elem_code == 11:
                return datetime.time(*data)
            elif self.elem_code == 13:
                return bool(data)
            elif self.elem_code == 18:
                return py3_get_string(data[1:])
            elif self.elem_code == 19:
                return py3_get_string(data[:-1])
            else:
                return data
        else:
            return None
