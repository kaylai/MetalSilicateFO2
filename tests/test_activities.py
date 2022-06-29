import unittest
import fO2calculate as fc
import pandas as pd

class Test_single_sample_gammas(unittest.TestCase):
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

		self.sample = fc.Sample(self.majors)
		self.temperature = 1300.0 # in C
		self.gamma_Fe_metal = 1.0050007836459114
		self.gamma_solutes_metal = {'S': 3.4335928934135906,
									'C': 0,
									'O': 0,
									'Ni': 0.6412192814527147,
									'Cu': 0,
									'Si': 0.0007963668649021159,
									'Mn': 1.1211931441281593,
									'Cr': 0.6117140306124951,
									'Ga': 0,
									'Nb': 0,
									'Ta': 0.04965273646150327}

	def test_calculate_gamma_Fe_metal_single(self):
		result = fc.calculate_gamma_Fe_metal(sample=self.sample,
											 temperature=self.temperature).result
		self.assertAlmostEqual(result, self.gamma_Fe_metal, places=4)

	def test_calculate_gamma_solute_metal_single(self):
		for spec in fc.core.standard_interactions:
			gamma = fc.calculate_gamma_solute_metal(sample=self.sample,
													temperature=self.temperature,
													species=spec,
													print_warnings=False).result
			self.assertAlmostEqual(gamma, self.gamma_solutes_metal[spec], places=4)

class Test_batch_sample_gammas(unittest.TestCase):
	def setUp(self):
		self.myfile = fc.BatchFile('test_data.xlsx')
		self.temperature = 1300.0 # in C

		self.gamma_Fe_metal_batch = {'NSA-36': 1.134278,
									 'MK15':   1.010205}

		self.gamma_solutes_NSA36 = {'S': 12.162458,
									'C': 0.000000,
									'O': 0.000000,
									'Ni': 0.000000,
									'Cu': 0.000000,
									'Si': 0.004036,
									'Mn': 0.896305,
									'Cr': 0.713339,
									'Ga': 0.000000,
									'Nb': 0.000000,
									'Ta': 0.000000}

		self.gamma_solutes_MK15 = {'S': 3.918409,
								   'C': 0.000000,
								   'O': 0.000000,
								   'Ni': 0.645856,
								   'Cu': 0.000000,
								   'Si': 0.000658,
								   'Mn': 1.278680,
								   'Cr': 1.009483,
								   'Ga': 0.000000,
								   'Nb': 0.000000,
								   'Ta': 0.074424}

	def test_calculate_gamma_Fe_metal_batch(self):
		result = self.myfile.calculate_gamma_Fe_metal(temperature=self.temperature,
													  print_warnings=False)
		res_dict = {index: row['gamma_Fe_metal'] for index, row in result.iterrows()}
		for key, value in res_dict.items():
			self.assertAlmostEqual(res_dict[key], self.gamma_Fe_metal_batch[key], places=4)

	def test_calculate_gamma_solute_metal_batch(self):
		result = self.myfile.calculate_gamma_solute_metal(temperature=self.temperature,
														  print_warnings=False)
		res_dict_NSA36 = {spec: result['gamma_' + spec + '_metal']['NSA-36'] for spec in fc.core.standard_interactions}
		res_dict_MK15 = {spec: result['gamma_' + spec + '_metal']['MK15'] for spec in fc.core.standard_interactions}
		for key, value in res_dict_NSA36.items():
			self.assertAlmostEqual(res_dict_NSA36[key], self.gamma_solutes_NSA36[key], places=4)
		for key, value in res_dict_MK15.items():
			self.assertAlmostEqual(res_dict_MK15[key], self.gamma_solutes_MK15[key], places=4)

