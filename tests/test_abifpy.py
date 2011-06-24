import unittest
import datetime
import abifpy


class TestAbif(object):

    def tearDown(self):
        del self.abif
        del self.trimmedSeq
        del self.untrimmedSeq

    def testFileType(self):
        self.assertEqual(self.abif._raw[:4], 'ABIF')

    def testAllTagsParsed(self):
        for key in self.abif.tags:
            # should be assertNone(), not available in py2.6
            if self.abif.tags[key].elemCode != 1024:
                assert (self.abif.get_data(key) is not None)

    def testTagDataLen(self):
        for key in self.abif.tags:
            code = self.abif.tags[key].elemCode
            data = self.abif.get_data(key)
            
            # 10 & 11 returns datetime, 12 is not clear, 1024 is ignored
            if code not in [10, 11, 12, 1024]:
                # account for null character in pString and cString
                mod = 1 if code in [18, 19] else 0
                # if data is int/float, len is always 1
                obtained = len(data) if type(data).__name__ not in ['int', 'float'] else 1
                expected = self.abif.tags[key].elemNum - mod
                
                self.assertEqual(obtained, expected)

    def testTrimIsShorter(self):
        self.assertTrue(len(self.trimmedSeq) <= len(self.untrimmedSeq))

    def testTrimIsSubset(self):
        self.assertTrue(self.trimmedSeq in self.untrimmedSeq)


class TestAbif3730(TestAbif, unittest.TestCase):

    def setUp(self):
        self.abif = abifpy.Trace('3730.ab1')
        self.trimmedSeq = self.abif.trim(self.abif.seq())
        self.untrimmedSeq = self.abif.seq()


class TestAbif3100(TestAbif, unittest.TestCase):

    def setUp(self):
        self.abif = abifpy.Trace('3100.ab1')
        self.trimmedSeq = self.abif.trim(self.abif.seq())
        self.untrimmedSeq = self.abif.seq()


class TestAbif310(TestAbif, unittest.TestCase):

    def setUp(self):
        self.abif = abifpy.Trace('310.ab1')
        self.trimmedSeq = self.abif.trim(self.abif.seq())
        self.untrimmedSeq = self.abif.seq()


class TestAbifFake(unittest.TestCase):

    def testFakeFileType(self):
        source = 'fake.ab1'
        self.assertRaises(IOError, abifpy.Trace, source)


class TestAbifEmpty(unittest.TestCase):

    def testShortSequenceNotTrimmed(self):
        self.abif = abifpy.Trace('empty.ab1')
        self.assertRaises(ValueError, self.abif.trim, self.abif.seq()) 


if __name__ == '__main__':
    unittest.main()
