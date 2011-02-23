======
ABIFPY
======

------------------------------------------
Python module for reading .ab1 trace files
------------------------------------------

abifpy is a python module that extracts sequence and various other data from
Applied Biosystem's, Inc. format (ABIF) file. The module was written based on
the `official spec`_ released by Applied Biosystems.

Usage
=====

::

    $ python
    >>> import abifpy
    >>> yummy = abifpy.Trace('tracefile.ab1')

Or if you want to perform base trimming directly: ::
    
    >>> yummy = abifpy.Trace('tracefile.ab1', trimming=True)

The module can be used with Biopython if it is installed. It provides the
following methods: ::

    yummy.write(out_file="", qual=0)       
    # writes a fasta (qual=0), qual (qual=1), or fastq (qual=1) file from the trace file
    # default output is tracefile.fa

    yummy.seq()
    # returns a string of nucleotide sequence as called by the basecaller

    yummy.qual(char=True)
    # returns a list of ascii characters of phred quality values (offset 33)
    # if char=False, the phred quality values is returned instead

    yummy.trim()        
    # trims the sequence using Richard Mott's algorithm (used in phred)
    # can be used for trimming quality values returned by yummy.qual() as well
    
    yummy.seqrecord()   
    # returns a SeqRecord object of the trace file

    yummy.get_dir()
    # returns a metadata stored in the file, accepts keys from yummy.tags (see below)
    # half-cooked method, not yet capable of extracting the entire file metadata

The file metadata (e.g. sample well, sequencing instrument) can be looked up in
the ``self.meta`` dictionary: ::

    yummy.meta['sampleid']    # string of sample ID entered before the run
    yummy.meta['well']        # string of well ID
    yummy.meta['instrument']  # string of sequencing machine model
    yummy.meta['id']          # string of trace file name

Keys for ``yummy.meta`` are the values of ``abifpy.TAGS``, except for ``'id'``.

Additionally, these attributes can also be accessed: ::

    yummy.tags        # dictionary of tags with values of directory contents
    yummy._header     # tuple of extracted header values
    yummy._raw        # string representation of file contents

You can get the metadata not contained in ``yummy.meta`` by using ``yummy.get_dir()``
with one of the keys in ``yummy.tags`` as the argument, e.g.::

    >>> yummy.get_dir('GTyp1')
    'POP7'

Be warned though that this method is half-cooked. Sometimes it returns the value you want,
other times it will only return ``None``. For more info on the file metadata, refer to 
the `official spec`_. 

Installation
============

Just add the abifpy directory to your ``$PYTHONPATH`` (in ``.bashrc`` to make it persistent).

License
=======

Copyright (c) 2011 by Wibowo Arindrarto

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

.. _official spec: http://www.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf
