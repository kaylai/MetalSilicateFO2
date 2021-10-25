import pandas as pd
import numpy as np
import warnings as w
import math
from scipy.signal import savgol_filter 

from copy import deepcopy, copy

from fO2calculate import core
from fO2calculate import fO2bufferplotter

from fO2calculate import tavern as tv

def calc_dIW_from_fO2(P, T, fO2):
	"""Translates from absolute (not Log) fO2 value into number of log units away from the
	Iron-Wüstite buffer (dIW).

	Parameters
	----------
	P 	float
		Pressure in bars

	T 	float
		Temperature in degrees C

	fO2 	float
		Absolute fO2 value

	Returns
	-------
	float
		fO2 in terms of dIW
	"""
	if fO2 <= 0:
		return np.nan
	else:
		P_GPa = P / 10000
		T_K = T + 273.15
		log_fO2 = math.log10(fO2)

		log_IW_value = fO2bufferplotter.buffers.calc_IW(P_GPa, T_K)
		dIW = log_fO2 - log_IW_value

		return dIW

def calc_dIW_from_fO2_Francis(T, fO2):
	"""Uses Francis McCubbin spreadsheet for equation for IW buffer

	Parameters
	----------
	T 	float
		Temperature in degrees C

	fO2 	float
		Absolute fO2 value

	Returns
	-------
	float
		fO2 in terms of dIW
	"""
	if fO2 <= 0:
		return np.nan
	else:
		T_K = T + 273.15
		log_fO2 = math.log10(fO2)

		log_IW_value = (2 * 
						math.log10((math.exp((-275723 + 
											  126.3154 * T_K - 
											  8.11767 * T_K * 
											  math.log(T_K))/
											 (8.314*T_K)))))
		dIW = log_fO2 - log_IW_value

		return dIW

def calc_gammas(temp, press, species='all'):
		"""Returns fugacity coefficients calculated using the Redlich Kwong Equation of State.
		Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 
		30 October 2003.

		Parameters
		----------
		temp: float
			Temperature in degrees C.

		press: float
			Pressure in bars.

		species: str
			Choose which species to calculate.
			Options are: 'CH4', CO', 'CO2', 'H2', 'H2O', 'H2S', 'O2', 'S2', 'SO2', or all.
			If all is passed, a dictionary of values is returned. Default value is 'all'.

		Returns
		-------
		float or dict
			Fugacity coefficient for passed species.
			If single species is passed, float.
			If "all" is passed, dictionary with keys 'CH4', CO', 'CO2', 'H2', 'H2O', 'H2S', 'O2',
			'S2', 'SO2'
		"""

		tempK = temp + 273.15
		R = 8.3145

		gamma_dict = {}

		for species in core.fluid_species_names:
			#Calculate a and b parameters (depend only on critical parameters)...
			a = 0.42748 * R**2.0 * core.critical_params[species]["cT"]**(2.5) / (core.critical_params[species]["cP"] * 10.0**5)
			b = 0.08664 * R * core.critical_params[species]["cT"] / (core.critical_params[species]["cP"] * 10.0**5)
			kappa = 0.0

			#Calculate coefficients in the cubic equation of state...
			#coeffs: (C0, C1, C2, A, B)
			A = a * press * 10.0**5 / (math.sqrt(tempK) * (R * tempK)**2.0)
			B = b * press * 10.0**5 / (R * tempK)
			C2 = -1.0
			C1 = A - B - B * B
			C0 = -A * B

			#Solve the cubic equation for Z0 - Z2, D...
			Q1 = C2 * C1 / 6.0 - C0 / 2.0 - C2**3.0 / 27.0
			P1 = C2**2.0 / 9.0 - C1 / 3.0
			D = Q1**2.0 - P1**3.0

			if D >= 0:
				kOneThird = 1.0 / 3.0

				absQ1PSqrtD = math.fabs(Q1 + math.sqrt(D))
				temp1 = absQ1PSqrtD**kOneThird
				temp1 *= (Q1 + math.sqrt(D)) / absQ1PSqrtD

				absQ1MSqrtD = math.fabs(Q1 - math.sqrt(D))
				temp2 = absQ1MSqrtD**kOneThird
				temp2 *= (Q1 - math.sqrt(D)) / absQ1MSqrtD

				Z0 = temp1 + temp2 - C2 / 3.0
			else:
				temp1 = Q1**2.0 / (P1**3.0)
				temp2 = math.sqrt(1.0 - temp1) / math.sqrt(temp1)
				temp2 *= Q1 / math.fabs(Q1)

				gamma = math.atan(temp2)

				if gamma < 0:
					gamma = gamma + math.pi 

				Z0 = 2.0 * math.sqrt(P1) * math.cos(gamma/3.0) - C2 / 3.0
				Z1 = 2.0 * math.sqrt(P1) * math.cos((gamma + 2.0 * math.pi) / 3.0) - C2/3.0
				Z2 = 2.0 * math.sqrt(P1) * math.cos((gamma + 4.0 * math.pi) / 3.0) - C2/3.0

				if Z0 < Z1:
					temp0 = Z0
					Z0 = Z1
					Z1 = temp0

				if Z1 < Z2:
					temp0 = Z1
					Z1 = Z2
					Z2 = temp0

				if Z0 < Z1:
					temp0 = Z0
					Z0 = Z1
					Z1 = temp0

			#Determine the fugacity coefficient of first root and departure functions...
			#calcdepfns(coeffs[3], 	coeffs[4], 	paramsab[0], 	Z[0])
			#calcdepfns(A, 			B, 			kappa, 			Z)

			#Calculate Departure Functions
			gamma = math.exp(Z0 - 1.0 - math.log(Z0-B) - A * math.log(1.0+B/Z0)/B)
			Hdep = R * tempK * (Z0 - 1.0 - 1.5*A*math.log(1.0+B/Z0)/B)
			Sdep = R * (math.log(Z0-B) - 0.5*A*math.log(1.0+B/Z0)/B)

			gamma_dict[species] = gamma

		if species is 'all':
			return gamma_dict
		elif isinstance(species, str):
			try:
				return gamma_dict[species]
			except:
				raise core.InputError('Passed species not recognized.')
		else:
			raise core.InputError('species must be type string.')
		

class Fluid(object):
	"""
	Creates a Fluid object with parameters defined here.
	"""

	def __init__(self, folder_name, DSC_filename='DSC.txt', masses_filename='masses.txt',
				 TG_filename='TG.txt', highT=False):
		""" Initiates the Fluid class

		Parameters
		----------
		folder_name	str
			Path to folder containing three text files with DSC, mass, and TG data. This is the
			standard data format for Netchze gas files.

		DSC_filename	str
			OPTIONAL. Default is 'DSC.txt'. Name of file with DSC data.

		masses_filename	str
			OPTIONAL. Default is 'masses.txt'. Name of file with masses data.

		TG_filename	str
			OPTIONAL. Default is 'TG.txt'. Name of file with TG data.

		highT	float
			OPTIONAL. Default is False. If float, this value will be used as the minimum
			temperature at which to keep gas data.
		"""
		# Read in files
		dsc_data = pd.read_csv(folder_name + '/' + DSC_filename, delimiter="\t") # DSC file
		masses_data = pd.read_csv(folder_name + '/' + masses_filename, delimiter="\t") # masses file
		tg_data = pd.read_csv(folder_name + '/' + TG_filename, delimiter="\t") # TG file

		# drop last column of masses (should be empty)
		masses_data = masses_data.iloc[:, :-1]

		# remove white space in column names
		dsc_data.columns = dsc_data.columns.str.replace(' ', '')
		masses_data.columns = masses_data.columns.str.replace(' ', '')
		tg_data.columns = tg_data.columns.str.replace(' ', '')

		# coerce all data to numeric
		dsc_data = dsc_data.apply(pd.to_numeric, errors='coerce')
		masses_data = masses_data.apply(pd.to_numeric, errors='coerce')
		tg_data = tg_data.apply(pd.to_numeric, errors='coerce')

		# drop all rows with any NaN
		dsc_data = dsc_data.dropna()
		masses_data = masses_data.dropna()
		tg_data = tg_data.dropna()

		# delete trailing data from cool-down phase
		max_temp_val = masses_data['temp_2'].max()
		max_temp_val = max_temp_val//1 # get floor (round down)
		masses_data = masses_data.loc[(masses_data['temp_2']<max_temp_val)]
		tg_data = tg_data.loc[(tg_data['temp_TG']<max_temp_val)]
		dsc_data = dsc_data.loc[(dsc_data['temp_DSC']<max_temp_val)]

		# if requested, drop all low temp data
		if isinstance(highT, float) or isinstance(highT, int):
			minT = highT
			masses_data = masses_data.loc[(masses_data['temp_2']>=minT)]
			tg_data = tg_data.loc[(tg_data['temp_TG']>=minT)]
			dsc_data = dsc_data.loc[(dsc_data['temp_DSC']>=minT)]

		self.max_temp_round = int(5 * round(max_temp_val/5)) #save max temp val rounded to nearest 5

		# normalize small numbers to large
		masses_data['ic_2'] = masses_data['ic_2'].apply(lambda x: x/10**-10)
		masses_data['ic_16'] = masses_data['ic_16'].apply(lambda x: x/10**-10)
		masses_data['ic_18'] = masses_data['ic_18'].apply(lambda x: x/10**-10)
		masses_data['ic_28'] = masses_data['ic_28'].apply(lambda x: x/10**-10)
		masses_data['ic_32'] = masses_data['ic_32'].apply(lambda x: x/10**-10)
		masses_data['ic_34'] = masses_data['ic_34'].apply(lambda x: x/10**-10)
		masses_data['ic_44'] = masses_data['ic_44'].apply(lambda x: x/10**-10)
		masses_data['ic_64'] = masses_data['ic_64'].apply(lambda x: x/10**-10)

		self.dsc_data = dsc_data
		self.masses_data = masses_data
		self.tg_data = tg_data

		self.highT = highT

		# subtract background from all masses we care about
		self.masses_data = self._subtract_gas_background()

	def get_smoothed_masses(self, window_length=11, polyorder=2):
		""" Applies scipy.signal.savgol_filter to all gas data

		Parameters
		----------
		window_length:	int
			After scipy.signal.savgol_filter. The length of the filter window (i.e., the number of
			coefficients). window_length must be a positive odd integer. If mode is ‘interp’,
			window_length must be less than or equal to the size of x.

		polyorder:	int
			After scipy.signal.savgol_filter. The order of the polynomial used to fit the samples.
			polyorder must be less than window_length.

		Returns
		-------
		pandas DataFrame
		"""
		_masses_data = self.masses_data.copy()
		for i, col_name in enumerate(core.fluid_col_names):
			x = _masses_data['temp_2'].tolist()
			gas = _masses_data[col_name].tolist()
			gas_smooth = savgol_filter(gas, window_length, polyorder)

			# drop non_smoothed column values
			_masses_data.drop(col_name, axis=1, inplace=True)

			# add smoothed values in place of dropped ones
			_masses_data[col_name] = gas_smooth

		return _masses_data

	def _subtract_gas_background(self):
		""" Finds minimum cps value for all species and subtracts this from all values in
		that species curve.

		Returns
		-------
		pandas DataFrame
		"""
		_masses_data = self.masses_data.copy()
		for i, col_name in enumerate(core.fluid_col_names):
			minvalue = _masses_data[col_name].min()
			_masses_data[col_name] -= minvalue

		return _masses_data

	def calc_fO2_from_gas_ratios(self, oxidized_species):
		"""
		Returns fO2 value of gas calculated using molar ratios of redox couples.

		Parameters
		----------
		oxidized_species	str
			Name of the oxidized species in the redox couple. Currently can take one of 'CO2'
			or 'H2O'. 'CO2' will perform calculation on reaction: CO + 1/2O2 = CO2. 'H2O' will
			perform calculation on reaction: H2 + 1/2O2 = H2O.

		Returns
		-------
		pandas DataFrame
			With temperature and computed fO2 values
		"""
		if oxidized_species is 'CO2':
			reduced_species = 'CO'
			temp_col = 'temp_28'
			ox_mass_col = 'ic_44'
			red_mass_col = 'ic_28'
			ox_MW = 44.01
			red_MW = 28
		elif oxidized_species is 'H2O':
			reduced_species = 'H2'
			temp_col = 'temp_2'
			ox_mass_col = 'ic_18'
			red_mass_col = 'ic_2'
			ox_MW = 18.02
			red_MW = 2.02
		else:
			raise core.InputError("oxidized_species must be one of 'CO2' or 'H2O'.")

		# Calculate absolute fO2 values
		fO2_vals = []
		temps = []
		for index, row in self.masses_data.iterrows():
			K_F = tv.calc_Ks(row[temp_col], species=oxidized_species)
			gamma_oxidized = calc_gammas(row[temp_col], press=1, species=oxidized_species)
			gamma_reduced = calc_gammas(row[temp_col], press=1, species=reduced_species)
			ox_moles = row[ox_mass_col] #* ox_MW
			red_moles = row[red_mass_col] #* red_MW
			molar_ratio = ox_moles / red_moles
			fO2_vals.append((molar_ratio * gamma_oxidized/gamma_reduced * 1/K_F)**2)
			temps.append(row[temp_col])

		abs_data = pd.DataFrame({'temp': temps, 'fO2': fO2_vals})

		# Calculate dIW values as well
		dIW_vals = []
		for index, row in abs_data.iterrows():
			dIW_vals.append(calc_dIW_from_fO2_Francis(row['temp'], row['fO2']))

		return_data = abs_data.copy()
		return_data['dIW'] = dIW_vals

		return return_data

	def calc_fCO(self):
		""" Calculate fCO as:
		fCO = KF * (XCH4/XH2)*(XH2O/XH2**2) * gammaCH4 * gammaH2O / gammaH2**3

		Derivation:
		-----------
		CH4 + H2O = CO + 3H2

		KF = XCO*gammaCO * (XH2*gammaH2)**3 / XCH4*gammaCH4 * XH2O*gammaH2O
		KF = fCO * (XH2*gammaH2)**3 / XCH4*gammaCH4 * XH2O*gammaH2O
		KF = fCO * (XH2/XCH4) * (gammaH2/gammaCH4) * (XH2**2/XH2O) * (gammaH2**2/gammaH2O)
		KF = fCO * (XH2/XCH4) * (XH2**2/XH2O) * (gammaH2**3/gammaCH4*gammaH2O)
		fCO = KF * (XCH4/XH2) * (XH2O/XH2**2) * gammaCH4 * gammaH2O / gammaH2**3

		Returns
		-------
		pandas DataFrame
			With temperature and computed fCO values
		"""
		fCO_vals = []
		temp_vals = []
		for index, row in self.masses_data.iterrows():
			temp = row['temp_16']
			K_F = tv.calc_Ks(temp, species='CH4')
			gammaH2 = calc_gammas(temp, press=1, species='H2')
			gammaCH4 = calc_gammas(temp, press=1, species='CH4')
			gammaH2O = calc_gammas(temp, press=1, species='H2O')
			temp_vals.append(temp)
			fCO = (K_F * (row['ic_16']/row['ic_2']) *
						 (row['ic_18']/(row['ic_2']**2)) *
						 gammaCH4 * gammaH2O / (gammaH2**3))
			fCO_vals.append(fCO)

		return_data = pd.DataFrame({'temp': temp_vals, 'fCO': fCO_vals})

		return return_data

	def calc_fCO_alt(self, fO2_how='CO2'):
		""" Calculate fCO as:
		fCO = KF * sqrt(fO2) * (XCH4/XH2**2) * (gammaCH4/gammaH2**2)

		Derivation:
		-----------
		CH4 + 1/2O2 = CO + 2H2

		KF = (fCO * XH2**2*gammaH2**2)/(XCH4*gammaCH4 * sqrt(fO2))
		KF = fCO * (1/sqrt(fO2)) * (XH2**2/XCH4) * (gammaH2**2)/(gammaCH4)
		fCO = KF * sqrt(fO2) * (XCH4/XH2**2) * (gammaCH4/gammaH2**2)

		Parameters:
		-----------
		fO2_how:	str
			Method to calculate fO2. Default is 'CO2', which uses calc_fO2_from_gas_ratios('CO2').
			Can also be 'H2O', which uses calc_fO2_from_gas_ratios('H2O').

		Returns
		-------
		pandas DataFrame
			With temperature and computed fCO values
		"""
		# get fO2 values
		fO2_data = self.calc_fO2_from_gas_ratios(fO2_how)

		# copy over masses_data
		_masses_data = self.masses_data.copy()

		# add fO2 values to masses_data
		data = pd.concat([_masses_data, fO2_data], axis=1)

		fCO_vals = []
		temp_vals = []
		for index, row in data.iterrows():
			temp = row['temp_16']
			K_F = tv.calc_Ks(temp, species='CH4_alt')
			gammaH2 = calc_gammas(temp, press=1, species='H2')
			gammaCH4 = calc_gammas(temp, press=1, species='CH4')
			fO2 = row['fO2']
			temp_vals.append(temp)
			fCO = (K_F * math.sqrt(fO2) *
				   (row['ic_16']/row['ic_2']**2) *
				   (gammaCH4/gammaH2**2))
			fCO_vals.append(fCO)

		return_data = pd.DataFrame({'temp': temp_vals, 'fCO': fCO_vals})

		return return_data

