#!/usr/bin/env python

import struct

# directory data structure:
# tag name, tag number, element type code, element size, number of elements
# data size, data offset, handle
FMT_DIR = '>4sIHHIIII'
# for header +file type, file version
FMT_HEAD = '>4sH4sIHHIII'
# only retrieve the data we care about:
# plate barcode, color data, machine model, sample tracking id, sequence, 
# peak locations, max quality value, trim probability
TAGS = ('CTID', 'DATA', 'HCFG', 'LIMS', 'PBAS', 'PLOC', 'phQL', 'phTR')


class Trace(object):

    """Class representing trace file"""

    def __init__(self, tfile):        
        with open(tfile) as source:
            self._data = source.read()

        if self._data[:4] == 'ABIF':
            self._header = struct.unpack(FMT_HEAD, self._data[:30])
            self.version = self._header[1]
            self.elemsize = self._header[5]
            self.elemnum = self._header[6]
            self.offset = self._header[8]
            self.tags = {}

            for entry in self.gen_dir():
                if entry[0] in TAGS:
                    self.tags[entry[0] + str(entry[1]).zfill(3)] = entry
        else:
            del self._data
            print "File error. Make sure the file is a proper .ab1 file."
    
    def gen_dir(self):
        index = 0
        while index < self.elemnum:
            start = self.offset + index * self.elemsize
            finish = self.offset + (index + 1) * self.elemsize
            yield struct.unpack(FMT_DIR, self._data[start:finish])
            index += 1


#class Seq(Trace):

# write parser for file header
# write parser for directory entries

# open file
## parse file header
## parse directory entries
# get sequence
## find directory entries related to sequence information
# get miscellaneous data
# write sequence
