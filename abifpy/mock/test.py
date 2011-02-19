#!/usr/bin/env python

import struct
import abifpy

x = abifpy.Trace('trace1.ab1')

print 'qual', x.qual
print 'seq', x.seq
print 'sample id', x.sampleid
print 'well', x.well
print 'plate id', x.plateid
print 'instrument', x.instrument
