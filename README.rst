======
ABIFPY
======

------------------------------------------
Python module for reading .ab1 trace files
------------------------------------------

abifpy is a python module that extracts sequence and various other data from Applied Biosystem's, Inc. format (ABIF) file. The module is made based on the `official spec`_ released by Applied Biosystem.

Usage
=====

::

 $ python
 >>> import abifpy
 >>> yummy = abifpy.Trace('tracefile.ab1')

By default, only these tags are extracted::

    yummy.seq         # sequence data as called by the basecaller, untrimmed
    yummy.qual        # quality values after basecalling
    yummy.sampleid    # sample id entered before the run
    yummy.well        # well coordinate of the sample (A1-K12, if using a 96-well plate)
    yummy.instrument  # sequencing machine model

You can invoke the ``all_tags=True`` option when instantiating the class to get all tags available. These tags can then be viewed using as the ``yummy.tags`` attribute. Be warned that the module are only able to read data from the five tags above. If you want to make sense of the extra tags, all the information is in the `official spec`_. 


License
=======

Copyright (c) 2011 by Wibowo Arindrarto

Permission is hereby granted, free of charge, to any person obtaining a copyof this software and associated documentation files (the "Software"), to dealin the Software without restriction, including without limitation the rightsto use, copy, modify, merge, publish, distribute, sublicense, and/or sellcopies of the Software, and to permit persons to whom the Software isfurnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included inall copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS ORIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THEAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHERLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS INTHE SOFTWARE.

.. _official spec: http://www.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf
