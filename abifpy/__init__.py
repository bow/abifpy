#!/usr/bin/env python

"""Python module for reading .ab1 trace files"""

import re
import struct

# directory data structure:
# tag name, tag number, element type code, element size, number of elements
# data size, data offset, handle
FMT_DIR = '>4sIHHIIII'
# for header +file type, file version
FMT_HEAD = '>4sH4sIHHIII'
# only retrieve the data we care about:
# plate barcode, machine model, sequence, quality values sample id, well number
TAGS = {'HCFG3':'instrument', 'PBAS2':'seq', 'PCON2':'qual', 'SMPL1':'sampleid', 'TUBE1':'well'}


class Trace(object):
    """Class representing trace file"""
    def __init__(self, in_file, trimming=False, all_tags=False):        
        try:
            with open(in_file) as source:
                self._data = source.read()
            if not self._data[:4] == 'ABIF':
                raise IOError('Input is not a valid .ab1 trace file')
        except IOError as (strerror):
            print "IOError: {0}".format(strerror)
        else:
            self._header = struct.unpack(FMT_HEAD, self._data[:30])
            self.version = self._header[1]
            self.tags = {}
            self.id = in_file.replace('.ab1','')

            # build dictionary of tags that we care about if all_tags=False
            # otherwise get all tags
            for entry in self._gen_dir(self._header):
                if not all_tags:
                    if (entry[0] + str(entry[1])) in TAGS:
                        self.tags[entry[0] + str(entry[1])] = entry
                else:
                    self.tags[entry[0] + str(entry[1])] = entry

            # retrieve attributes from tags
            for item in self.tags:
                self._extract_dir(self.tags[item])

            if trimming:
                self.trim()
    
    def _gen_dir(self, head_entry):
        """Generator for directory contents."""
        head_elemsize = head_entry[5]
        head_elemnum = head_entry[6]
        head_offset = head_entry[8]
        
        index = 0
        while index < head_elemnum:
            start = head_offset + index * head_elemsize
            finish = head_offset + (index + 1) * head_elemsize
            # added directory offset to tuple
            # to handle directories with data size <= 4 bytes
            yield struct.unpack(FMT_DIR, self._data[start:finish]) + (start,)
            index += 1

    def _extract_dir(self, dir_entry):
        """Extracts data from directories in the file."""
        # if data size is <= 4 bytes, data is stored inside the directory
        # so offset needs to be changed
        tag_name, tag_no, elem_code, elem_size, elem_no, dir_size, dir_offset, dir_handle, data_offset = dir_entry

        if dir_size <= 4:
            dir_offset = data_offset + 20

        if elem_code == 2:
            fmt = str(dir_size) + 's'
            data = struct.unpack(fmt, 
                    self._data[dir_offset:dir_offset+dir_size])[0]
            if tag_name == 'PCON':
                self.qual = self._get_qual(data)
            elif tag_name == 'PBAS':
                # replaces semi-ambiguous DNA characters with 'N'
                self.seq = re.sub("K|Y|W|M|R|S",'N',data)
       
        elif elem_code == 18:
            fmt = str(dir_size-1) + 's'
            data = struct.unpack(fmt,
                    self._data[dir_offset+1:dir_offset+dir_size])[0]
            if tag_name == 'SMPL':
                self.sampleid = data
            elif tag_name == 'TUBE':
                self.well = data
       
        elif elem_code == 19:
            fmt = str(dir_size-1) + 's'
            data = struct.unpack(fmt,
                    self._data[dir_offset:dir_offset+dir_size-1])[0]
            if tag_name == 'HCFG':
                self.instrument = data
    
    def _get_qual(self, qual):
        """Returns a list of read quality values."""
        qual_list = []
        for value in qual:
            qual_list.append(ord(value))
        return qual_list

    def _get_score(self, cutoff):
        """Returns a list of read probability scores."""
        score_list = []
        for qual in self.qual:
            # calculate probability back from formula used
            # to calculate phred qual values
            prob = 10 ** (qual/-10.0)
            score_list.append(cutoff - prob)
        return score_list

    def seqrecord(self):
        """Returns a SeqRecord object of the trace file."""
        try:
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            biopython = True
        except ImportError:
            biopython = False

        if biopython:
            return SeqRecord(Seq(self.seq), id=self.id, name="",
                      description=self.sampleid)
        else:
            print 'Biopython was not detected. No SeqRecord was created.'
            return None

    def write(self, output=""):
        """Writes the trace file sequence to a fasta file.
        
        Keyword argument:
        output -- output file name (detault 'tracefile'.fa)

        """
        if output == "":
            output = self.id + '.fa'

        with open(output, 'rw') as out_file:
            contents = '>{0} {1}\n{2}\n'.format(
                        self.id, self.sampleid, self.seq)
            out_file.writelines(contents)

    def trim(self, segment=20, cutoff=0.05):
        """Trims the sequence using Richard Mott's modified trimming algorithm.
        
        Keyword argument:
        cutoff -- probability cutoff value
        window -- segment size for calculating segment score

        Trimmed bases are determined from their segment score, ultimately
        determined from each base's quality values. If a segment score is
        negative, the entire segment will be trimmed.
        
        """
        # set flag for trimming
        take = False
        trim_start = 0
        trim_finish = len(self.seq)
        # obtain scores of each base (based on qual values)
        self.score = self._get_score(cutoff)
        
        if len(self.score) <= segment:
            print "Sequence length for {0} is shorter than trim segment size (20). Sequence not trimmed.".format(self.id)
        else:
            for index in xrange(len(self.seq)-segment):
                # calculate segment score
                # if segment score is negative, trim the region out
                segment_score = reduce(lambda x,y:x+y,
                                 self.score[index:index+segment])
                # algorithm for obtaining trimming position values
                if not take and segment_score > 0:
                    trim_start = index
                    take = True
                elif take and segment_score <= 0:
                    # if segment length is longer than remaining bases
                    # take up everything until the end
                    if index + segment <= len(self.seq):
                        trim_finish = index + segment - 1
                        break
                    else:
                        break
                
            self.seq = self.seq[trim_start:trim_finish]
            self.qual = self.qual[trim_start:trim_finish]
            self.score = self.score[trim_start:trim_finish]
