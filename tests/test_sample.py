import unittest
import fO2calculate as fc
import numpy as np
import pandas as pd

class TestCreateSample(unittest.TestCase):
    def setUp(self):
        self.majors = pd.Series({'SiO2':    37.46,
                                 'TiO2':    0.34,
                                 'Al2O3':   8.76,
                                 'Cr2O3':   0.01,
                                 'FeO':     0.36,
                                 'MgO':     36.61,
                                 'CaO':     9.97,
                                 'Na2O':    3.24,
                                 'K2O':     0.24,
                                 'MnO':     0.59,
                                 'CoO':     0.004,
                                 'NiO':     0.014,
                                 'Ta2O5':   1.25,
                                 'WO3':     0.006,
                                 'ReO3':    0.024,
                                 'Si':      1.69,
                                 'Cr':      0.15,
                                 'Re':      7.98,
                                 'Fe':      80.55,
                                 'Mn':      0.41,
                                 'Co':      2,
                                 'Ni':      2.01,
                                 'W':       2.84,
                                 'Ta':      0.06,
                                 'S':       1.014,
                    })

        self.majors_normed = pd.Series({'SiO2':    37.885,
                                        'TiO2':    0.344,
                                        'Al2O3':   8.859,
                                        'Cr2O3':   0.010,
                                        'FeO':     0.364,
                                        'MgO':     37.025,
                                        'CaO':     10.083,
                                        'Na2O':    3.277,
                                        'K2O':     0.243,
                                        'MnO':     0.597,
                                        'CoO':     0.004,
                                        'NiO':     0.014,
                                        'Ta2O5':   1.264,
                                        'WO3':     0.006,
                                        'ReO3':    0.024,
                                        'Si':      1.712,
                                        'Cr':      0.152,
                                        'Re':      8.085,
                                        'Fe':      81.608,
                                        'Mn':      0.415,
                                        'Co':      2.026,
                                        'Ni':      2.036,
                                        'W':       2.877,
                                        'Ta':      0.061,
                                        'S':       1.027
                                        })

        # Mol fractions calculated externally
        self.majors_mol = pd.Series({'SiO2':    0.3332,
                                     'TiO2':    0.0023,
                                     'Al2O3':   0.0459,
                                     'Cr2O3':   0.0000,
                                     'FeO':     0.0027,
                                     'MgO':     0.4854,
                                     'CaO':     0.0950,
                                     'Na2O':    0.0279,
                                     'K2O':     0.0014,
                                     'MnO':     0.0044,
                                     'CoO':     0.0000,
                                     'NiO':     0.0001,
                                     'Ta2O5':   0.0015,
                                     'WO3':     0.0000,
                                     'ReO3':    0.0001,
                                     'Si':      0.0352,
                                     'Cr':      0.0017,
                                     'Re':      0.0251,
                                     'Fe':      0.8441,
                                     'Mn':      0.0044,
                                     'Co':      0.0418,
                                     'Ni':      0.0200,
                                     'W':       0.0090,
                                     'Ta':      0.0002,
                                     'S':       0.0185
                                    })

        self.sample = fc.Sample(self.majors)

    def test_createSample(self):
        for ox in self.majors.index:
            self.assertEqual(self.sample._composition[ox],self.majors[ox])

    def test_setdefault_noargs(self):
        self.assertEqual(self.sample.default_normalization,'none')
        self.assertEqual(self.sample.default_units,'wtpt')

    def test_setdefaults_none_wtpt(self):
        sample = fc.Sample(self.majors, default_normalization='none',default_units='wtpt')
        self.assertEqual(sample.default_normalization,'none')
        self.assertEqual(sample.default_units,'wtpt')

    def test_setdefaults_standard_mol(self):
        sample = fc.Sample(self.majors, default_normalization='standard',default_units='mol')
        self.assertEqual(sample.default_normalization,'standard')
        self.assertEqual(sample.default_units,'mol')

    def test_setdefaults_garbageNorm(self):
        with self.assertRaises(fc.core.InputError):
            fc.Sample(composition=self.majors,default_normalization='garbage')

    def test_setdefaults_garbageType(self):
        with self.assertRaises(fc.core.InputError):
            fc.Sample(composition=self.majors,default_units='garbage')

    def test_type_garbage(self):
        with self.assertRaises(fc.core.InputError):
            fc.Sample(composition=self.majors,units='garbage')

    def test_type_wtpt(self):
        sample = fc.Sample(self.majors,units='wtpt')
        for ox in self.majors.index:
            self.assertEqual(self.sample._composition[ox],self.majors[ox])

    def test_type_mol(self):
        sample = fc.Sample(self.majors_mol,units='mol')
        for ox in self.majors.index:
            self.assertEqual(np.round(sample.get_composition(species=ox, units='mol'),4),np.round(self.majors_mol[ox],4))


class TestGetComposition(unittest.TestCase):
    def setUp(self):
        self.majors = pd.Series({'SiO2':    37.46,
                                 'TiO2':    0.34,
                                 'Al2O3':   8.76,
                                 'Cr2O3':   0.01,
                                 'FeO':     0.36,
                                 'MgO':     36.61,
                                 'CaO':     9.97,
                                 'Na2O':    3.24,
                                 'K2O':     0.24,
                                 'MnO':     0.59,
                                 'CoO':     0.004,
                                 'NiO':     0.014,
                                 'Ta2O5':   1.25,
                                 'WO3':     0.006,
                                 'ReO3':    0.024,
                                 'Si':      1.69,
                                 'Cr':      0.15,
                                 'Re':      7.98,
                                 'Fe':      80.55,
                                 'Mn':      0.41,
                                 'Co':      2,
                                 'Ni':      2.01,
                                 'W':       2.84,
                                 'Ta':      0.06,
                                 'S':       1.014,
                    })

        self.majors_normed = pd.Series({'SiO2':    37.885,
                                        'TiO2':    0.344,
                                        'Al2O3':   8.859,
                                        'Cr2O3':   0.010,
                                        'FeO':     0.364,
                                        'MgO':     37.025,
                                        'CaO':     10.083,
                                        'Na2O':    3.277,
                                        'K2O':     0.243,
                                        'MnO':     0.597,
                                        'CoO':     0.004,
                                        'NiO':     0.014,
                                        'Ta2O5':   1.264,
                                        'WO3':     0.006,
                                        'ReO3':    0.024,
                                        'Si':      1.712,
                                        'Cr':      0.152,
                                        'Re':      8.085,
                                        'Fe':      81.608,
                                        'Mn':      0.415,
                                        'Co':      2.026,
                                        'Ni':      2.036,
                                        'W':       2.877,
                                        'Ta':      0.061,
                                        'S':       1.027
                                        })

        self.majors_normed_additionalVolatiles = pd.Series({'SiO2':    37.885,
                                                            'TiO2':    0.344,
                                                            'Al2O3':   8.859,
                                                            'Cr2O3':   0.010,
                                                            'FeO':     0.364,
                                                            'MgO':     37.025,
                                                            'CaO':     10.083,
                                                            'Na2O':    3.277,
                                                            'K2O':     0.243,
                                                            'MnO':     0.597,
                                                            'CoO':     0.004,
                                                            'NiO':     0.014,
                                                            'Ta2O5':   1.264,
                                                            'WO3':     0.006,
                                                            'ReO3':    0.024,
                                                            'Si':      1.730,
                                                            'Cr':      0.154,
                                                            'Re':      8.169,
                                                            'Fe':      82.455,
                                                            'Mn':      0.420,
                                                            'Co':      2.047,
                                                            'Ni':      2.058,
                                                            'W':       2.907,
                                                            'Ta':      0.061,
                                                            'S':       1.014
                                                            })

        self.majors_normed_fixedVolatiles = pd.Series({'SiO2':    37.885,
                                                       'TiO2':    0.344,
                                                       'Al2O3':   8.859,
                                                       'Cr2O3':   0.010,
                                                       'FeO':     0.364,
                                                       'MgO':     37.025,
                                                       'CaO':     10.083,
                                                       'Na2O':    3.277,
                                                       'K2O':     0.243,
                                                       'MnO':     0.597,
                                                       'CoO':     0.004,
                                                       'NiO':     0.014,
                                                       'Ta2O5':   1.264,
                                                       'WO3':     0.006,
                                                       'ReO3':    0.024,
                                                       'Si':      1.712,
                                                       'Cr':      0.152,
                                                       'Re':      8.086,
                                                       'Fe':      81.619,
                                                       'Mn':      0.415,
                                                       'Co':      2.027,
                                                       'Ni':      2.037,
                                                       'W':       2.878,
                                                       'Ta':      0.061,
                                                       'S':       1.014
                                                       })

        # Mol fractions calculated externally
        self.majors_mol = pd.Series({'SiO2':    0.3332,
                                     'TiO2':    0.0023,
                                     'Al2O3':   0.0459,
                                     'Cr2O3':   0.0000,
                                     'FeO':     0.0027,
                                     'MgO':     0.4854,
                                     'CaO':     0.0950,
                                     'Na2O':    0.0279,
                                     'K2O':     0.0014,
                                     'MnO':     0.0044,
                                     'CoO':     0.0000,
                                     'NiO':     0.0001,
                                     'Ta2O5':   0.0015,
                                     'WO3':     0.0000,
                                     'ReO3':    0.0001,
                                     'Si':      0.0360,
                                     'Cr':      0.0017,
                                     'Re':      0.0256,
                                     'Fe':      0.8630,
                                     'Mn':      0.0045,
                                     'Co':      0.0203,
                                     'Ni':      0.0205,
                                     'W':       0.0092,
                                     'Ta':      0.0002,
                                     'S':       0.0189
                                    })

        self.majors_mol_exclV = pd.Series({'SiO2':    0.3332,
                                           'TiO2':    0.0023,
                                           'Al2O3':   0.0459,
                                           'Cr2O3':   0.0000,
                                           'FeO':     0.0027,
                                           'MgO':     0.4854,
                                           'CaO':     0.0950,
                                           'Na2O':    0.0279,
                                           'K2O':     0.0014,
                                           'MnO':     0.0044,
                                           'CoO':     0.0000,
                                           'NiO':     0.0001,
                                           'Ta2O5':   0.0015,
                                           'WO3':     0.0000,
                                           'ReO3':    0.0001,
                                           'Si':      0.0367,
                                           'Cr':      0.0018,
                                           'Re':      0.0261,
                                           'Fe':      0.8796,
                                           'Mn':      0.0046,
                                           'Co':      0.0207,
                                           'Ni':      0.0209,
                                           'W':       0.0094,
                                           'Ta':      0.0002,
                                           'S':       0.0
                                          })

        self.sample = fc.Sample(self.majors)

    def test_default(self):
        composition = self.sample.get_composition()
        for ox in composition.index:
            self.assertEqual(composition[ox],self.majors[ox])

    def test_normnone(self):
        composition = self.sample.get_composition(normalization='none')
        for ox in composition.index:
            self.assertEqual(composition[ox],self.majors[ox])

    def test_normnone_exclV(self):
        composition = self.sample.get_composition(normalization='none',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(composition[ox],self.majors[ox])

    def test_normstd(self):
        composition = self.sample.get_composition(normalization='standard')
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_normed[ox],3))

    def test_normstd_exclV(self):
        composition = self.sample.get_composition(normalization='standard',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_normed[ox],3))

    def test_normfixedVolatiles(self):
        composition = self.sample.get_composition(normalization='fixedvolatiles')
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_normed_fixedVolatiles[ox],3))

    def test_wtptoxides_fixedVolatiles_exclV(self):
        composition = self.sample.get_composition(normalization='fixedvolatiles',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_normed_fixedVolatiles[ox],3))

    def test_wtptoxides_additionalVolatiles(self):
        composition = self.sample.get_composition(normalization='additionalvolatiles')
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_normed_additionalVolatiles[ox],3))

    def test_wtptoxides_additionalVolatiles_exclV(self):
        composition = self.sample.get_composition(normalization='additionalvolatiles',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_normed_additionalVolatiles[ox],3))

    def test_mol(self):
        composition = self.sample.get_composition(units='mol')
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_mol[ox],3))

    def test_moloxides_exclV(self):
        composition = self.sample.get_composition(units='mol',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_mol_exclV[ox],3))
