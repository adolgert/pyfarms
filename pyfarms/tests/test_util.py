import tempfile
from unittest import TestCase

import pyfarms.util as util

class TestUtil(TestCase):
    def test_nonexist(self):
        with self.assertRaises(RuntimeError):
            util.check_filename("nonexistent", "hopefully not existent")

    def test_exist(self):
        f=tempfile.NamedTemporaryFile()
        util.check_filename(f.name, "this module file.")


if __name__=='__main__':
    a=TestUtil()
    a.test_nonexist()
    a.test_exist()
    