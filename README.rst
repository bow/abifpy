======
ABIFPY
======

------------------------------------------
Python module for reading .ab1 trace files
------------------------------------------

abifpy is a python module that extracts sequence and various other data from
Applied Biosystem's, Inc. format (ABIF) file. The module was written based on
the `official spec`_ released by Applied Biosystems.

The module provides the following items:

*class* Trace(inFile)
    Class representing the trace file ``inFile``.

Trace Objects
=============

seq
    String of base-called nucleotide sequence stored in the file.

qual
    String of phred quality characters of the base-called sequence.

qualVal
    List of phred quality values of the base-called sequence.

id
    String of the sequence file name.

name
    String of the sample name entered prior to sequencing.

trim(sequence[, cutoff=0.05])        
    Returns a trimmed sequence using Richard Mott's algorithm (used in phred)
    with the probability cutoff of 0.05. Can be used on ``seq``, ``qual``, and
    ``qualVal``.
    
get_data(key)
    Returns metadata stored in the file, accepts keys from ``tags`` (see below).

export([outFile="", fmt='fasta'])       
    Writes a fasta (``fmt='fasta'``), qual (``fmt='qual'``), or 
    fastq (``fmt='fastq'``) file from the trace file. Default format is ``fasta``.

close()
    Closes the Trace file object.

seq_remove_ambig(seq)
    Replaces extra ambigous base characters (K, Y, W, M, R, S) with 'N'. Accepts ``seq``
    for input.

EXTRACT
    Dictionary for determining which metadata are extracted.

data
    Dictionary that contains the file metadata. The keys are values of ``EXTRACT``.

tags
    Dictionary of tags with values of data directory class instance. Keys are tag name and 
    tag number, concatenated. Use ``get_data()`` to access values in each ``tags`` entry.

Usage
=====

::

    $ python
    >>> import abifpy
    >>> yummy = abifpy.Trace('tracefile.ab1')

Or if you want to perform base trimming directly::
    
    >>> yummy = abifpy.Trace('tracefile.ab1', trimming=True)

Sequence can be accessed with the ``seq`` attribute. Other attributes of note
are ``qual`` for phred quality characters, ``qualVal`` for phred quality values,
``id`` for sequencing trace file name, and ``name`` for the sample name::

    >>> yummy.seq
    'CCAAGGTGCAGACTTCCATCT'
    >>> yummy.qual
    '3824DESHSSSSS:DSSSSSS'
    >>> yummy.qualVal
    [18, 23, 17, 19, 35, 36, 50, 39, 50, 50, 50, 50, 50, 25, 35, 50, 50, 50, 50, 50, 50]
    >>> yummy.id
    'tracefile'
    >>> yummy.name
    'TOPO_clone1_fwd'

If trimming was not performed when instantiating, you can still do it afterwards::
    
    >>> yummy.trim(yummy.seq)

The quality values itself can be trimmed as well::

    >>> yummy.trim(yummy.qual)

Viewing the trace file metadata is easy. Use the values from ``EXTRACT``
as the keys in ``data``::

    >>> yummy.data['well']
    'B6'
    >>> yummy.data['model']
    '3730'
    >>> yummy.data['run start date']
    datetime.date(2009, 12, 10)

metadata not contained in ``data`` can be viewed using ``get_data()``
with one of the keys in ``tags`` as the argument, e.g.::

    >>> yummy.get_data('PTYP1')
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
