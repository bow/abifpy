#!/usr/bin/env python

import struct

fmt_head = '>4sH4sIHHIII'
fmt_dir = '>4sIHHIIII'

class Trace(object):

    """Class representing trace file"""

    def __init__(self, tfile):        
        with open(tfile) as source:
            self._data = source.read()

        if self._data[:4] == 'ABIF':
            self._header = struct.unpack(fmt_head, self._data[:30])
            self.version = self._header[1]
            self.elemsize = self._header[5]
            self.elemnum = self._header[6]
            self.offset = self._header[8]
            self.tags = {}

            for i in self.gen_dir():
                if (i[:4], 

        else:
            del self._data
            print "File error. Make sure the file is a proper .ab1 file."
    
    def gen_dir(self):
        index = 0
        while index < self.elemnum:
            start = self.offset + index * self.elemsize
            finish = self.offset + (index + 1) * self.elemsize
            yield struct.unpack(fmt_dir, self._data[start:finish])
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
