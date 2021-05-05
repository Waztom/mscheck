import unittest
import sys
import os

from mscheck import MassSpectrum


class SpectrumTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mzMLFile = os.path.join(os.getcwd(), "tests", "testdata", "1AB-1001.mzML")
        cls.spectrum_positive = MassSpectrum(mzMLfilepath=mzMLFile, mode="Positive")
        cls.spectrum_negative = MassSpectrum(mzMLfilepath=mzMLFile, mode="Negative")

    def test_create_positive_spectrum(self):
        self.assertIsNotNone(self.spectrum_positive)
        self.assertEqual(len(self.spectrum_positive.MSdata), 3)
        self.assertEqual(sum(self.spectrum_positive.MSdata["TIC"]), 520031)

    def test_create_negative_spectrum(self):
        self.assertIsNotNone(self.spectrum_negative)
        self.assertEqual(len(self.spectrum_negative.MSdata), 3)
        self.assertEqual(sum(self.spectrum_negative.MSdata["TIC"]), 301270)


if __name__ == "__main__":
    unittest.main()
