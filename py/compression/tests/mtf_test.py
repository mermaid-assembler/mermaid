import unittest

import mtf

# WRONG!

class TestMTF(unittest.TestCase):
    def test_mtf(self):
        out = mtf.mtf("ACGATGTA")
        self.assertEqual(out, [0,1,2,2,3,2,1,2])

    def test_imtf(self):
        input = [0,1,3,3,1,2,2,3]
        out = mtf.imtf(input)
        self.assertEqual(''.join(out), "ACTGTCGA")

    def test_identity(self):
        encoding = mtf.mtf("ACTGTGTA")
        decoding = mtf.imtf(encoding)
        self.assertEqual(''.join(decoding), "ACTGTGTA")

if __name__ == '__main__':
    unittest.main()
