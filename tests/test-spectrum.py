import unittest
import sys
import os

from mscheck.spectrum import MassSpectrum


class GetMoleculesDataTest(unittest.TestCase):
    def test_create_spectrum(self):

        mzMLFile = os.path.join("tests", "testdata", "1AB-1001.mzML")
        spectrum = mscheck.MassSpectrum(mzMLfilepath=mzMLFile)
        self.assertIsNotNone(spectrum)


if __name__ == "__main__":
    unittest.main()
