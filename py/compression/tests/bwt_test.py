import unittest

import bwt

class TestBWT(unittest.TestCase):
    def test_bwt(self):
        out = bwt.bwt('banana')
        self.assertEqual(out, 'annb$aa')

    def test_ibwt(self):
        input = 'annb$aa'
        out = bwt.ibwt(input)
        self.assertEqual(out, 'banana')

    def test_identity(self):
        encoding = bwt.bwt("baseball")
        decoding = bwt.ibwt(encoding)
        self.assertEqual(decoding, "baseball")

if __name__ == '__main__':
    unittest.main()
