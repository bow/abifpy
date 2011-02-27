#!/usr/bin/env python

"""Python module for reading .ab1 trace files"""

import re
import struct

# dictionary for deciding which values to extract and contain in self.meta
TAGS = {
            'HCFG3':'instrument', 
            'SMPL1':'sampleid', 
            'TUBE1':'well',
       } 

__version__ = '0.3'

class Trace(object):
    """Class representing trace file"""
    def __init__(self, in_file, trimming=False):        
        try:
            with open(in_file) as source:
                self._raw = source.read()
            if not self._raw[:4] == 'ABIF':
                raise IOError('Input is not a valid .ab1 trace file')
        except IOError as (strerror):
            print "IOError: {0}".format(strerror)
        else:
            # header data structure:
            # file type, file, version, tag name, tag number, element type code,
            # element size, number of elements, data size, data offset, handle,
            # file type, file version
            fmt_head = '>4sH4sIHHIII'
            # dictionary for containing file metadata
            self.meta = {}
            # dictionary for containing extracted directory data
            self.data = {}
            self.meta['id'] = in_file.replace('.ab1','')
            self.trimming = trimming
            # values contained in file header
            self._header = struct.unpack(fmt_head, self._raw[:30])
            # file format version
            self.version = self._header[1]

            # build dictionary of data tags
            for entry in self._gen_dir():
                self.data[entry[0] + str(entry[1])] = entry

            # retrieve attributes from tags
            for item in self.data.keys():
                # only extract data from tags we care about
                if item in TAGS:
                    # e.g. self.meta['well'] = 'B6'
                    self.meta[TAGS[item]] = self.get_data(item)
    
    def _gen_dir(self):
        """Generator for directory contents."""
        # directory data structure:
        # file type, file, version, tag name, tag number, 
        # element type code, element size, number of elements
        # data size, data offset, handle
        fmt_head = '>4sIHHIIII'
        head_elemsize = self._header[5]
        head_elemnum = self._header[6]
        head_offset = self._header[8]
        index = 0
        
        while index < head_elemnum:
            start = head_offset + index * head_elemsize
            finish = head_offset + (index + 1) * head_elemsize
            # added directory offset to tuple
            # to handle directories with data size <= 4 bytes
            yield struct.unpack(fmt_head, self._raw[start:finish]) + (start,)
            index += 1

    def get_data(self, dir_entry):
        """Extracts data from directories in the file."""
        tag_name, tag_no, elem_code, elem_size, elem_no, dir_size, dir_offset, dir_handle, data_offset = self.data[dir_entry]

        # if data size is <= 4 bytes, data is stored inside the directory
        # so offset needs to be changed
        if dir_size <= 4:
            dir_offset = data_offset + 20

        if elem_code == 2:
            fmt = str(dir_size) + 's'
            data = struct.unpack(fmt, 
                    self._raw[dir_offset:dir_offset+dir_size])[0]
        elif elem_code == 18:
            fmt = str(dir_size-1) + 's'
            data = struct.unpack(fmt,
                    self._raw[dir_offset+1:dir_offset+dir_size])[0]
        elif elem_code == 19:
            fmt = str(dir_size-1) + 's'
            data = struct.unpack(fmt,
                    self._raw[dir_offset:dir_offset+dir_size-1])[0]
        else:
            data = None
        return data
    
    def seq(self):
        """Returns sequence contained in the trace file."""
        data = self.get_data('PBAS2')
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
        data = self.get_data('PCON2')
        qual_list = []

        if not char:    
            for value in data:
                qual_list.append(str(ord(value)))
        else:
            for value in data:
                if ord(value) > 93:
                    value = chr(93)
                qual_list.append(chr(ord(value) + 33))
        # return value is list to make it compatible with the char option
        # and with the write() function (for writing .qual files)
        if self.trimming:
            return self.trim(qual_list)
        else:
            return qual_list

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

    def write(self, out_file="", qual=0):
        """Writes the trace file sequence to a fasta file.
        
        Keyword argument:
        out_file -- output file name (detault 'tracefile'.fa)
        qual -- 0: write fasta file, 1: write qual file, 2: write fastq file

        """
        
        if out_file == "":
            file_name = self.meta['id']
            if qual == 0:
                file_name += '.fa'
            elif qual == 1:
                file_name += '.qual'
            elif qual == 2:
                file_name += '.fq'
        else:
            file_name = out_file
        
        if qual == 0:
            contents = '>{0} {1}\n{2}\n'.format(
                        file_name, self.meta['sampleid'], self.seq())
        elif qual == 1:
            contents = '>{0} {1}\n{2}\n'.format(
                        file_name, self.meta['sampleid'], ' '.join(self.qual(char=False)))
        elif qual == 2:
            contents = '@{0} {1}\n{2}\n+{0} {1}\n{3}\n'.format(
                        file_name, self.meta['sampleid'], self.seq(), ''.join(self.qual()))

        with open(file_name, 'w') as out_file:
            out_file.writelines(contents)

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
        trim_start = 0
        trim_finish = len(seq)
        
        if len(seq) <= segment:
            print "Sequence length for {0} is shorter than trim segment size ({1}). Sequence not trimmed.".format(self.meta['id'], segment)
        else:
            # calculate score for trimming
            score_list =  []
            # actually the same as seq.qual(), but done to avoid infinite recursion
            qual_hack = [ord(x) for x in self.get_data('PCON2')]
            
            for qual in qual_hack:
                # calculate probability back from formula used
                # to calculate phred qual values
                prob = 10 ** (qual/-10.0)
                score_list.append(cutoff - prob)
            # algorithm for obtaining trimming position values
            for index in xrange(len(seq)-segment):
                # calculate segment score
                # if segment score is negative, trim the region out
                segment_score = reduce(lambda x,y:x+y,
                                 score_list[index:index+segment])
                if not take and segment_score > 0:
                    trim_start = index
                    take = True
                elif take and segment_score <= 0:
                    # if segment length is longer than remaining bases
                    # take up everything until the end
                    if index + segment <= len(seq):
                        trim_finish = index + segment - 1
                        break
                    else:
                        break
            return seq[trim_start:trim_finish]
