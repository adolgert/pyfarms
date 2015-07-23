import tempfile
from unittest import TestCase

import util

class TestUtil(TestCase):
    def test_nonexist(self):
        with self.assertRaises(RuntimeError):
            util.check_filename("nonexistent")

    def test_exist(self):
        util.check_filename("test_util.py")


if __name__=='__main__':
    a=TestUtil()
    a.test_nonexist()
    a.test_exist()
    