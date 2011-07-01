import unittest
import datetime
import abifpy


class TestAbif(object):

    def tearDown(self):
        del self.abif
        del self.trimmedSeq
        del self.untrimmedSeq

    def testFileType(self):
        assert self.abif._raw[:4] == 'ABIF'

    def testTagData(self):
        for key in self.abif.tags:
            # should be assertNone(), not available in py2.6
            if self.abif.tags[key].elemCode != 1024:
                assert self.abif.get_data(key) != None

    def testTagDataLen(self):
        for key in self.abif.tags:
            code = self.abif.tags[key].elemCode
            data = self.abif.get_data(key)
            
            # 10 & 11 returns datetime, 12 is not clear, 1024 is ignored
            # only check for strings, arrays, and numbers
            if code not in [10, 11, 12, 1024]:
                # account for null character in pString and cString
                mod = 1 if code in [18, 19] else 0
                # if data is int/float, len is always 1
                obtained = len(data) if not isinstance(data, (int, float, bool)) else 1
                expected = self.abif.tags[key].elemNum - mod
                
                assert obtained == expected
                
    def testTagDataType(self):
        for key in self.abif.tags:
            code = self.abif.tags[key].elemCode
            data = self.abif.get_data(key)

            # user data should return None
            if code == 1024:
                assert data == None
            # check for string return type in tags 2, 18, 19
            elif code in [2, 18, 19]:
                assert isinstance(data, basestring)
            # check for datetime return tag
            elif code == 10:
                assert isinstance(data, datetime.date)
            # check for datetime return tag
            elif code == 11:
                assert isinstance(data, datetime.time)
            elif code == 13:
                assert isinstance(data, bool)
            # check for number return types
            # some tags' data are still in a tuple of numbers, so will have to
            # iterate over them
            elif isinstance(data, tuple):
                for item in data:
                    assert isinstance(item, (int, float))
            # otherwise just check for type directly
            else:
                assert isinstance(data, (int, float))

    def testTrimIsShorter(self):
        assert len(self.trimmedSeq) < len(self.untrimmedSeq)


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
