import unittest
import sys
import os

from mscheck import AnalyseSpectrum
from mscheck import utils


class AnalyseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mzMLFile = os.path.join(os.getcwd(), "tests", "testdata", "1AB-1001.mzML")
        cls.spectrum = AnalyseSpectrum(mzMLfilepath=mzMLFile, mode="Positive")
        cls.spectrum.analyse(
            compoundsmiles="O=C(c1ccco1)N1CCN(C(=O)N2CCN(c3ccccc3)CC2)CC1",
            ionstoadd=["[H]", "[Na]", "[K]", "[NH4+]"],
            tolerance=1,
        )

    def test_create_positive_spectrum(self):
        self.assertIsNotNone(self.spectrum)
        self.assertEqual(len(self.spectrum.MSpeakdata), 4)
        self.assertEqual(sum(self.spectrum.MSpeakdata["TIC"]), 435108)

    def test_get_max_mz(self):
        mz_data = self.spectrum.MSpeakdata["mz_data"][10]
        max_mz = self.spectrum.get_max_mz(mz_data)
        self.assertEqual(max_mz, 169)

    def test_match_data(self):
        ions_matched = self.spectrum.Matchdata["ions"]
        expected_match = [("[H]", 369), ("[Na]", 391)]
        self.assertEqual(ions_matched, expected_match)
        mz_strongest = self.spectrum.Matchdata["mz_strongest"][0][0]
        expected_mz_strongest = [
            101.1,
            106.9,
            112.9,
            181.2,
            268.4,
            369.2,
            370.2,
            371.1,
            390.9,
            391.3,
        ]
        self.assertTrue((mz_strongest == expected_mz_strongest).all())

    def test_get_match_indices(self):
        ions_to_add = ["[H]", "[Na]", "[K]", "[NH4+]"]
        mass_ions = [utils.get_MW(utils.get_mol(smiles)) for smiles in ions_to_add]
        self.assertEqual(mass_ions, [1, 23, 39, 18])
        parent_masses = [self.spectrum.compound_MW + mass for mass in mass_ions]
        self.assertEqual(parent_masses, [369, 391, 407, 386])
        test_matches = [self.spectrum.get_match_indices(mass) for mass in parent_masses]
        known_indices = [(369, [63, 64, 65, 66]), (391, [118]), (407, []), (386, [])]
        for test_match, known in zip(test_matches, known_indices):
            self.assertEqual(test_match, known)


if __name__ == "__main__":
    unittest.main()
