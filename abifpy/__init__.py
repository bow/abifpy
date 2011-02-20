#!/usr/bin/env python

import struct

# directory data structure:
# tag name, tag number, element type code, element size, number of elements
# data size, data offset, handle
FMT_DIR = '>4sIHHIIII'
# for header +file type, file version
FMT_HEAD = '>4sH4sIHHIII'
# only retrieve the data we care about:
# plate barcode, machine model, sequence, quality values sample id, well number
TAGS = {'CTID1':'plateid', 'HCFG3':'instrument', 'PBAS2':'seq', 'PCON2':'qual', 'SMPL1':'sampleid', 'TUBE1':'well'}


class Trace(object):

    """Class representing trace file"""

    def __init__(self, tfile, all_tags=False):        
        with open(tfile) as source:
            self._data = source.read()

        if self._data[:4] == 'ABIF':
            self._header = struct.unpack(FMT_HEAD, self._data[:30])
            self.version = self._header[1]
            self._elemsize = self._header[5]
            self._elemnum = self._header[6]
            self._offset = self._header[8]
            self.tags = {}

            # build dictionary of tags that we care about if all_tags=False
            # otherwise get all tags
            for entry in self._gen_dir():
                if not all_tags:
                    if (entry[0] + str(entry[1])) in TAGS:
                        self.tags[entry[0] + str(entry[1])] = entry
                else:
                    self.tags[entry[0] + str(entry[1])] = entry

            # retrieve attributes from tags
            for item in self.tags:
                self._decode_dir(self.tags[item])
        else:
            del self._data
            print "File error. Make sure the file is a proper .ab1 file."
    
    # generator for directories
    def _gen_dir(self):
        index = 0
        while index < self._elemnum:
            start = self._offset + index * self._elemsize
            finish = self._offset + (index + 1) * self._elemsize
            yield struct.unpack(FMT_DIR, self._data[start:finish]) + (start,)
            index += 1

    # method to decode contents in a directory structure
    def _decode_dir(self, 
                    (tag_name, tag_no, elem_code, elem_size, elem_no,
                     dir_size, dir_offset, dir_handle, data_offset)):
        if elem_code == 2:
            fmt = str(elem_no) + 's'
            data = struct.unpack(fmt, 
                    self._data[dir_offset:dir_offset+dir_size])[0]
            if tag_name == 'PCON':
                self.qual = data
            elif tag_name == 'PBAS':
                self.seq = data
       
        elif elem_code == 18:
            # if data size is <= 4 byte, data is stored inside the directory
            # so offset needs to be changed
            if dir_size <= 4:
                dir_offset = data_offset + 20
            fmt = str(elem_no-1) + 's'
            data = struct.unpack(fmt,
                    self._data[dir_offset+1:dir_offset+elem_no])[0]
            if tag_name == 'SMPL':
                self.sampleid = data
            elif tag_name == 'TUBE':
                self.well = data
       
        elif elem_code == 19:
            fmt = str(elem_no-1) + 's'
            data = struct.unpack(fmt,
                    self._data[dir_offset:dir_offset+dir_size-1])[0]
            if tag_name == 'CTID':
                self.plateid = data
            elif tag_name == 'HCFG':
                self.instrument = data
     
        else:
            pass    

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
