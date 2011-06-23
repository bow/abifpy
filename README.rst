======
ABIFPY
======

------------------------------------------
Python module for reading .ab1 trace files
------------------------------------------

abifpy is a python module that extracts sequence and various other data from
Applied Biosystem's, Inc. format (ABIF) file. The module was written based on
the `official spec`_ released by Applied Biosystems.

The module provides the following items: ::

*class* abifpy.Trace(inFile)
    Class representing the trace file ``inFile``.

*class* abifpy.TraceDir(tagEntry, rawData)
    Class representing directory data in the trace file. ``tagEntry`` is
    the a tuple for unpacked directory data and ``rawData`` is the contents
    of the trace file. You would not normally need this unless you are
    playing around with the file metadata.

seq()
    Returns a string of nucleotide sequence from the trace file.

seqrecord()   
    Returns a SeqRecord object of the trace file (if Biopython is installed).

qual([char=True])
    Returns a list of ascii characters of phred quality values (offset 33)
    if ``char=False``, the phred quality values is returned instead.

trim(sequence[, cutoff=0.05])        
    Returns a trimmed sequence using Richard Mott's algorithm (used in phred)
    with the probability cutoff of 0.05, can be used for trimming quality
    values returned by ``qual()`` as well.
    
getData(key)
    Returns metadata stored in the file, accepts keys from ``dir`` (see below).

export([outFile="", qual=0])       
    Writes a fasta (``qual=0``), qual (``qual=1``), or fastq (``qual=1``) file
    from the trace file. Default output is ``tracefile.fa``.

TAGS
    Dictionary for determining which metadata are extracted.

meta
    Dictionary that contains the file metadata. The keys are values of ``TAGS``,
    except for ``id`` which is the trace file name.

dir
    Dictionary of tags with values of data directory class instance. Keys are tag name and 
    tag number, concatenated. Use ``getData()`` to access values in each ``dir`` entry.

Usage
=====

::

    $ python
    >>> import abifpy
    >>> yummy = abifpy.Trace('tracefile.ab1')

Or if you want to perform base trimming directly::
    
    >>> yummy = abifpy.Trace('tracefile.ab1', trimming=True)

Sequencing results can be accessed as string or as a Biopython SeqRecord object::

    >>> yummy.seq()
    'GCCAAGGTGCAGACTTCCATCT'
    >>> yummy.seqrecord()
    SeqRecord(seq=Seq('GCCAAGGTGCAGACTTCCATCT', Alphabet()), id='tracefile1', name='', description='tracefile1 seq', dbxrefs=[])

If trimming was not performed when instantiating, you can still do it afterwards::
    
    >>> yummy.trim(yummy.seq())

The quality values itself can be trimmed as well::

    >>> yummy.trim(yummy.qual())

Viewing the trace file metadata is easy::

    >>> yummy.meta['sampleid']
    'TOPO_clone1_F'
    >>> yummy.meta['well']
    'B6'
    >>> yummy.meta['model']
    '3730'

Metadata not contained in ``meta`` can be viewed using ``getData()``
with one of the keys in ``dir`` as the argument, e.g.::

    >>> yummy.getData('PTYP1')
    '96-well'

For more info on the meaning of these tags and the file metadata, consult the `official spec`_. 

Installation
============

Just add the abifpy directory to your ``$PYTHONPATH`` (in ``.bashrc`` to make it persistent).

License
=======

abifpy is licensed under the MIT License.

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
