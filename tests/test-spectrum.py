import unittest
import sys
import os

import carmc


class GetMoleculesDataTest(unittest.TestCase):
    # def setUp(self):
    #     sys.path.insert(0, "../carmc")
    #     from spectrum import MassSpectrum

    def test_create_spectrum(self):
        mzMLFile = os.path.join("tests", "testdata", "1AB-1001.mzML")
        spectrum = MassSpectrum(mzMLfilepath=mzMLFile)
        self.assertIsNotNone(spectrum)


if __name__ == "__main__":
    unittest.main()
