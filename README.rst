======
ABIFPY
======

------------------------------------------
Python module for reading .ab1 trace files
------------------------------------------

abifpy is a python module that extracts sequence and various other data from Applied Biosystem's, Inc. format (ABIF) file. The module was written based on the `official spec`_ released by Applied Biosystem.

Usage
=====

::

    $ python
    >>> import abifpy
    >>> yummy = abifpy.Trace('tracefile.ab1')

By default, only these data are extracted::

    yummy.seq         # string of untrimmed sequence data as called by the basecaller
    yummy.qual        # list of quality values after basecalling
    yummy.sampleid    # string of sample ID entered before the run
    yummy.well        # string of well ID
    yummy.instrument  # string of sequencing machine model
    yummy.id          # string of trace file name

Additionally, these attributes can also be accessed::

    yummy.tags        # dictionary of extracted tags with values corresponding to directory structure
    yummy._header     # tuple of extracted header values
    yummy._data       # string representation of file contents

You can invoke the ``all_tags=True`` option when instantiating the class to get all tags available. These tags can then be viewed with the ``yummy.tags`` attribute. Be warned that the module can only read data from the extracted tags above. If you want to make sense of the extra tags, refer to the `official spec`_. 

The module can be used with Biopython if it is installed. It provides the following methods::

    yummy.seqrecord()   # returns a SeqRecord object of the trace file
    yummy.write()       # writes a fasta file of the trace file

License
=======

Copyright (c) 2011 by Wibowo Arindrarto

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

.. _official spec: http://www.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf
