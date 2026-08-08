import unittest
import os
import numpy as np
from mscheck import AnalyseMS
from mscheck import utils


class AnalyseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mzMLFile = os.path.join(os.getcwd(), "tests", "testdata", "1AB-1001.mzML")
        cls.spectrum = AnalyseMS(filepath=mzMLFile, mode="Positive")
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
        ions_matched = self.spectrum.analysedata["ions"]
        expected_match = [("[H]", 369), ("[Na]", 391), ("[NH4+]", 387)]
        self.assertEqual(ions_matched, expected_match)
        mz_strongest = self.spectrum.analysedata["mz_strongest"][0][0]
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

        # analyse() stores matched ion/mass pairs and their retention times in
        # self.analysedata; recover the peak indices those RT values came from.
        matched_indices = {}
        for (_ion, mass), rt_values in zip(
            self.spectrum.analysedata["ions"], self.spectrum.analysedata["RT"]
        ):
            indices = [
                int(np.where(self.spectrum.MSpeakdata["RT"] == rt)[0][0])
                for rt in rt_values
            ]
            matched_indices[mass] = sorted(indices)

        # [K] (407) and the theoretical [NH4+] mass (386) did not match exactly;
        # [NH4+] instead matched at 387 within the search tolerance.
        known_indices = {369: [63, 64, 65, 66], 391: [118], 387: [17, 48, 55]}
        self.assertEqual(matched_indices, known_indices)


if __name__ == "__main__":
    unittest.main()
