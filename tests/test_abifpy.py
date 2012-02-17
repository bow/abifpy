#!/usr/bin/env python
#
# test_abifpy.py
# python2 unit tests for abifpy
# http://github.com/bow/abifpy

import unittest
import datetime
import abifpy

class TestAbif(unittest.TestCase):

    def __init__(self, filename):
        unittest.TestCase.__init__(self, methodName='runTest')
        self.filename = filename
        self.abif = abifpy.Trace(self.filename)
        self.trimmed_seq = self.abif.trim(self.abif.seq)
        self.untrimmed_seq = self.abif.seq

    def shortDescription(self):
        return "Testing %s" % self.filename

    def runTest(self):
        self.file_type()
        self.tag_data()
        self.tag_data_len()
        self.tag_data_type()
        self.trim_is_shorter() 

    def file_type(self):
        print "\nChecking file type..."
        self.abif._handle.seek(0)
        self.assertEqual(self.abif._handle.read(4), 'ABIF')

    def tag_data(self):
        print "Checking tag data parsing..."
        for key in self.abif.tags:
            # should be assertNone(), not available in py2.6
            if self.abif.tags[key].elem_code != 1024:
                self.assertNotEqual(self.abif.get_data(key), None)

    def tag_data_len(self):
        print "Checking tag data lengths..."
        for key in self.abif.tags:
            code = self.abif.tags[key].elem_code
            data = self.abif.get_data(key)
            
            # 10 & 11 returns datetime, 12 is not clear, 1024 is ignored
            # only check for strings, arrays, and numbers
            if code not in [10, 11, 12, 1024]:
                # account for null character in pString and cString
                mod = 1 if code in [18, 19] else 0
                # if data is int/float, len is always 1
                obtained = len(data) if not isinstance(data, (int, float, bool)) else 1
                expected = self.abif.tags[key].elem_num - mod
                self.assertEqual(obtained, expected)
                
    def tag_data_type(self):
        print "Checking tag data types..."
        for key in self.abif.tags:
            code = self.abif.tags[key].elem_code
            data = self.abif.get_data(key)

            # user data should return None
            if code == 1024:
                self.assertEqual(data, None)
            # check for string return type in tags 2, 18, 19
            elif code in [2, 18, 19]:
                self.assertTrue(isinstance(data, basestring))
            # check for datetime return tag
            elif code == 10:
                self.assertTrue(isinstance(data, datetime.date))
            # check for datetime return tag
            elif code == 11:
                self.assertTrue(isinstance(data, datetime.time))
            elif code == 13:
                self.assertTrue(isinstance(data, bool))
            # check for number return types
            # some tags' data are still in a tuple of numbers, so will have to
            # iterate over them
            elif isinstance(data, tuple):
                for item in data:
                    self.assertTrue(isinstance(item, (int, float)))
            # otherwise just check for type directly
            else:
                self.assertTrue(isinstance(data, (int, float)))

    def trim_is_shorter(self):
        print "Checking trimmed sequence length..."
        self.assertTrue(len(self.trimmed_seq) < len(self.untrimmed_seq))

class TestAbifFake(unittest.TestCase):

    def __init__(self, filename):
        unittest.TestCase.__init__(self, methodName='runTest')
        self.filename = filename

    def shortDescription(self):
        return "Testing %s" % self.filename

    def runTest(self):
        self.fake_file_type()

    def fake_file_type(self):
        print "\nIOError is raised if file is not ABIF..."
        self.assertRaises(IOError, abifpy.Trace, self.filename)

class TestAbifEmpty(unittest.TestCase):

    def __init__(self, filename):
        unittest.TestCase.__init__(self, methodName='runTest')
        self.filename = filename
        self.abif = abifpy.Trace(self.filename)

    def shortDescription(self):
        return "Testing %s" % self.filename

    def runTest(self):
        self.short_sequence_untrimmed()

    def short_sequence_untrimmed(self):
        print "\nValueError is raised if sequence length is shorter than 20."
        self.assertRaises(ValueError, self.abif.trim, self.abif.seq) 


abif_real = ['3730.ab1', '3100.ab1', '310.ab1',]
abif_fake = ['fake.ab1',]
abif_empty = ['empty.ab1',]

def run_suite():
    suite = unittest.TestSuite([TestAbif(n) for n in abif_real])
    suite.addTests(unittest.TestSuite([TestAbifFake(n) for n in abif_fake]))
    suite.addTests(unittest.TestSuite([TestAbifEmpty(n) for n in abif_empty]))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(run_suite())
