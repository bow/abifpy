import unittest
import datetime
import abifpy

class TestAbif(unittest.TestCase):

    def __init__(self, filename):
        unittest.TestCase.__init__(self, methodName='runTest')
        self.filename = filename
        self.abif = abifpy.Trace(self.filename)
        self.trimmedSeq = self.abif.trim(self.abif.seq())
        self.untrimmedSeq = self.abif.seq()

    def runTest(self):
        self.file_type()
        self.tag_data()
        self.tag_data_len()
        self.tag_data_type()
        self.trim_is_shorter() 

    def file_type(self):
        """File type is ABIF."""
        self.assertEqual(self.abif._raw[:4], 'ABIF')

    def tag_data(self):
        """All file tags returns data."""
        for key in self.abif.tags:
            # should be assertNone(), not available in py2.6
            if self.abif.tags[key].elemCode != 1024:
                self.assertNotEqual(self.abif.get_data(key), None)

    def tag_data_len(self):
        """All file tags returns the right length of data."""
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
                
                self.assertEqual(obtained, expected)
                
    def tag_data_type(self):
        """All file tag data is of the correct type."""
        for key in self.abif.tags:
            code = self.abif.tags[key].elemCode
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
        """Trimming shortens sequence."""
        self.assertTrue(len(self.trimmedSeq) < len(self.untrimmedSeq))

class TestAbifFake(unittest.TestCase):

    def test_fake_file_type(self):
        """IOError is raised if file is not ABIF."""
        source = 'tests/fake.ab1'
        self.assertRaises(IOError, abifpy.Trace, source)

class TestAbifEmpty(unittest.TestCase):

    def test_short_sequence_untrimmed(self):
        """ValueError is raised if sequence length is shorter than 20."""
        self.abif = abifpy.Trace('tests/empty.ab1')
        self.assertRaises(ValueError, self.abif.trim, self.abif.seq()) 


abif_real = ['tests/3730.ab1', 'tests/3100.ab1', 'tests/310.ab1',]
abif_fake = ['tests/fake.ab1',]
abif_empty = ['tests/empty.ab1',]

def run_suite():
    return unittest.TestSuite([TestAbif(n) for n in abif_real])

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(run_suite())
