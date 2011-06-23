import unittest
import abifpy

class TestAbif(object):

    def tearDown(self):
        del self.abif
        del self.trimmedSeq
        del self.untrimmedSeq

    def testFileType(self):
        self.assertEqual(self.abif._raw[:4], 'ABIF')

    def testAllTagsParsed(self):
        for key in self.abif.dir:
            # should be assertNone(), not available in py2.6
            assert (self.abif.getData(key) is not None)

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
