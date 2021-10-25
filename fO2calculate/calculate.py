from mendeleev import element
import numpy as np
from abc import abstractmethod

from fO2calculate import core
from fO2calculate import batchfile
from fO2calculate import sample_class
from fO2calculate import interactionparameters as ip
from fO2calculate import activitycoefficients as ac


class Calculate(object):
	""" The Calculate object is a template for implementing user-friendly
	methods for running calculations using the various models implemented
	here. Results of the calculation are always returned by accessing
	self.result.
	"""
	def __init__(self, sample, **kwargs):
		"""
		Initializes the calculation.

		Parameters
		----------
		sample:	Sample class
			Composition as a sample object.
		"""
		self.sample = sample
		self.result = self.calculate(sample=self.sample, **kwargs)

	@abstractmethod
	def calculate(self):
		""" """

class calculate_dIW(Calculate):
	""" Calculates the fO2 of a sample in terms of number of log units from the
	iron-wüstite buffer (dIW).

	Parameters
	----------
	sample:	Sample class
		Composition of silicate melt and metal phases as a sample object.

	temperature:	float
		Temperature in degrees C.

	interactions:	list
		OPTIONAL. List of strings of element names. Elements are solutes in a
		metal Fe liquid alloy for which interaction parameters are known and
		for which the user wishes to calculate the effects of interaction within
		the alloy. Elements need not be infinitely dilute. Compositional and
		temperature ranges for which interaction parameters are known are given
		in the interaction_parameters script within this library.
	"""
	def calculate(self, sample, temperature,
				  interactions=core.standard_interactions, **kwargs):
		_silicate_comp = sample.get_composition(how='silicate')
		_metal_comp = sample.get_composition(how='metal')

		dIW = calc_dIW(_silicate_comp, _metal_comp, temperature=temperature,
					   interactions=interactions, **kwargs)

		return dIW

def calc_dIW(silicate_comp, metal_comp, temperature=None, gammaFe="calculate", 
			 gammaFeO="calculate", interactions=core.standard_interactions,
			 print_warnings=False):
	"""
	Returns fO2 in terms of delta Iron-Wustite. Calculation is performed using
	mole fractions and activity coefficients of Fe in the metal and FeO in the
	silicate.

	Parameters
	----------
	silicate_comp: dict
		Dictionary of the composition of the silicate in wt% oxides.

	metal_comp: dict
		Dictionary of compositional information only for a metal, in terms of
		wt% elements

	temperature: float
		Temperature in degrees C.

	gammaFe and gammaFeO: float
		OPTIONAL. Default is "calculate" in which case the gammaFe or gammaFeO value will be
		calculated here. If the gammaFe or gammaFeO value is already known, it can be passed here
		to avoid having to calculate it again.

	interactions: list
		OPTIONAL. List of strings of element names. Elements are solutes in a metal Fe liquid
		alloy for which interaction parameters are known and for which the user wishes to calculate
		the effects of interaction within the alloy. Elements need not be infinitely dilute.
		Compositional and temperature ranges for which interaction parameters are known are given
		in the interaction_parameters script within this library.

	Returns
	-------
	float
		logfO2 in terms of delta Iron-Wustite
	"""
	silicate_sample = sample_class.Sample(silicate_comp)
	metal_sample = sample_class.Sample(metal_comp)
	molfrac_silicate = silicate_sample.get_composition(units='mol')
	molfrac_metal = metal_sample.get_composition(units='mol')

	X_FeO_silicate = molfrac_silicate['FeO']
	X_Fe_metal = molfrac_metal['Fe']

	if gammaFeO == "calculate":
		gamma_FeO_silicate = calc_gamma_FeO_silicate(silicate_comp)
	else:
		gamma_FeO_silicate = gammaFeO

	if gammaFe == "calculate":
		if temperature == None:
			raise InputError("Error: no temperature given. Cannot calculate " +
							 "gammaFe without temperature.")
		gamma_Fe_metal = calc_gamma_Fe_metal(metal_comp=metal_comp, temperature=temperature, 
											 interactions=interactions, 
											 print_warnings=print_warnings)
	else:
		gamma_Fe_metal = gammaFe

	return 2*np.log10(X_FeO_silicate/X_Fe_metal) + 2*np.log10(gamma_FeO_silicate/gamma_Fe_metal)

def calc_ln_gamma_naught_at_temperature(species, temperature):
	"""
	Calculates the reference value for the activity coefficient gamma at the given temperature
	based on a known value for gamma at a reference temperature. NOTA BENE: if there is no
	tabulated value for the reference gamma, ideality will be assumed (the reference gamma will be
	set equal to 1).

	Parameters
	----------
	species: str
		String of the name of the element for which to calculate gamma_naught

	temperature: float
		Temperature at which to calculate the reference gamma, in degrees C.
	"""
	user_T_K = temperature + 273.15

	try:
		ref_gamma = ac.activity_coeffs.loc[species]['gamma_naught_i']
		ref_T_K = ac.activity_coeffs.loc[species]['Temp_K']

		if isinstance(ref_gamma, float) == False and isinstance(ref_gamma, int) == False:
			ln_gamma_naught = 1
		else:
			ln_gamma_naught = (ref_T_K/user_T_K)*np.log(ref_gamma)
	except:
		ln_gamma_naught = 1

	return ln_gamma_naught

def calc_epsilon_at_temperature(i, j, temperature):
	"""
	Calculates the value for the interaction parameter epsilon between elements i and j at the
	given temperature based on a known value for e at a reference temperature.

	Parameters
	----------
	i: str
		String of the name of the first of two elements for which to calculate epsilon

	j: str
		String of the name of the second of two elements for which to calculate epsilon

	temperature: float
		Temperature at which to calculate the reference gamma, in degrees C.
	"""
	i_j = i + ' ' + j
	ref_e = ip.interaction_params.loc[i_j]['ei(j,k)']

	#convert e_i_j to epsilon_i_j using Equation 20 in Steelmaking Data Sourcebook
	epsilon_i_j_ref = (ref_e/0.00434 - 1)*(element(j).atomic_weight/55.85)+1

	#Convert epsilon to value at sample temp using equation 14 from Corgne et al. (2008)
	ref_T_K = ip.interaction_params.loc[i_j]['Temp_K']
	user_T_K = temperature + 273.15
	epsilon_i_j = (ref_T_K/user_T_K)*epsilon_i_j_ref

	return epsilon_i_j

def calc_gamma_Fe_metal(metal_comp, temperature,
						interactions=core.standard_interactions, 
						print_warnings=True):
	"""
	Calculates the activity coefficient, gamma, for iron. Interaction parameters epsilon are 
	computed for all elements passed to interactions, so long as interaction parameter values 
	are known.

	Parameters
	----------
	metal_comp: Sample object
		Dictionary of compositional information only for a metal, in terms of wt% elements

	temperature: float
		Temperature at which to perform the calculation, in degrees C.

	interactions: list
		OPTIONAL. List of strings of element names. Elements are solutes in a metal Fe liquid alloy
		for which interaction parameters are known and for which the user wishes to calculate the
		effects of interaction within the alloy. Elements need not be infinitely dilute.
		Compositional and temperature ranges for which interaction parameters are known are given
		in the interaction_parameters script within this library.
	
	print_warnings: bool
			OPTIONAL. Default is True. If set to True, any warnings related to the lack of
			compositional data or interaction parameters will be printed.
	"""
	# Note: Exclude Fe from this calculation by default

	metal_sample = sample_class.Sample(metal_comp)
	molfrac = metal_sample.get_composition(units='mol')

	def calc_indiv_term1(i, molfrac, epsilon_i_i):
		return epsilon_i_i*(molfrac[i] + np.log(1-molfrac[i]))

	def calc_indiv_term2(j, k, molfrac, epsilon_j_k):
		return epsilon_j_k*molfrac[j]*molfrac[k]*(1+((np.log(1-molfrac[j]))/molfrac[j]) +
												  (np.log(1-molfrac[k])/molfrac[k]))

	def calc_indiv_term3(i, k, molfrac, epsilon_i_k):
		return epsilon_i_k*molfrac[i]*molfrac[k]*(1+((np.log(1-molfrac[k]))/molfrac[k]) -
												  (1/(1-molfrac[i])))

	def calc_indiv_term4(j, k, molfrac, epsilon_j_k):
		return epsilon_j_k*molfrac[j]**2*molfrac[k]**2*((1/(1-molfrac[j])) + (1/(1-molfrac[k])) - 1)

	def calc_indiv_term5(i, k, molfrac, epsilon_i_k):
		return epsilon_i_k*molfrac[i]**2*molfrac[k]**2*((1/(1-molfrac[i])) + (1/(1-molfrac[k])) +
														(molfrac[i]/(2*(1-molfrac[i])**2)) - 1)

	user_interactions = interactions.copy()
	interactions = []
	skipped_elements = []
	for element in user_interactions:
		if element in metal_comp and metal_comp[element]>0:
			interactions.append(element)
		else:
			skipped_elements.append(element)

	interactions_NminusOne = interactions[:-1]
	interactions_JplusOne = interactions[1:]

	term1 = 0
	term2 = 0
	term3 = 0
	term4 = 0
	term5 = 0

	no_interaction_param = []
	for i in interactions:
		skip_this = False
		try:
			epsilon_i_i = calc_epsilon_at_temperature(i, i, temperature)
		except:
			skip_this = True
			param_string = str(i) + " and " + str(i)
			no_interaction_param.append(param_string)
		if skip_this == False:
			term1 += calc_indiv_term1(i, molfrac, epsilon_i_i)

	for j in interactions_NminusOne:
		for k in interactions_JplusOne:
			skip_this = False
			try:
				epsilon_j_k = calc_epsilon_at_temperature(j, k, temperature)
			except:
				skip_this = True
				param_string = str(j) + " and " + str(k)
				no_interaction_param.append(param_string)
			if skip_this == False:
				term2 += calc_indiv_term2(j, k, molfrac, epsilon_j_k)
				term4 += calc_indiv_term4(j, k, molfrac, epsilon_j_k)
			

	for i in interactions:
		for k in interactions:
			skip_this = False
			if i == k:
				pass
			else:
				try:
					epsilon_i_k = calc_epsilon_at_temperature(i, k, temperature)
				except:
					skip_this = True
					param_string = str(i) + " and " + str(k)
					no_interaction_param.append(param_string)
				if skip_this == False:
					term3 += calc_indiv_term3(i, k, molfrac, epsilon_i_k)
					term5 += calc_indiv_term5(i, k, molfrac, epsilon_i_k)

	ln_gamma_Fe = term1 - term2 + term3 + 0.5*term4 - term5
	gamma_Fe = np.exp(ln_gamma_Fe)

	if print_warnings == True:
		if len(skipped_elements) > 0:
			w.warn("No compositional information given for: " +
				', '.join(skipped_elements) +
				". These elements were skipped.",RuntimeWarning,stacklevel=2)
		if len(no_interaction_param) > 0:
			w.warn("No interaction parameters known for: " +
				', '.join(no_interaction_param),RuntimeWarning,stacklevel=2)

	return gamma_Fe

def calc_gamma_solute_metal(metal_comp, species, temperature, gammaFe="calculate",
							interactions=core.standard_interactions,
							print_warnings=True, **kwargs):
	"""
	Calculates the activity coefficient, gamma, for any solutes in an Fe-rich metal alloy.
	Interaction parameters epsilon are computed for all elements passed to interactions, so long as
	interaction parameter values are known.

	Parameters
	----------
	metal_comp: dict
		Dictionary of compositional information only for a metal, in terms of wt% elements

	species: str
		String of the name of the element for which to calculate gamma

	temperature: float
		Temperature at which to perform the calculation, in degrees C.

	gammaFe: float
		OPTIONAL. Default is "calculate" in which case the gammaFe value will be calculated here.
		If the gammaFe value is already known, it can be passed here to avoid having to calculate
		it again.

	elements: list
		OPTIONAL. List of elements for which to calculate gamma values, if that list is different
		than the list of interactions. Default value is None, in which case the elements
		list = interactions.

	interactions: list
		OPTIONAL.List of strings of element names. Elements are solutes in a metal Fe liquid alloy
		for which interaction parameters are known and for which the user wishes to calculate the
		effects of interaction within the alloy. Elements need not be infinitely dilute.
		Compositional and temperature ranges for which interaction parameters are known are given
		in the interaction_parameters script within this library.
	
	print_warnings: bool
			OPTIONAL. Default is True. If set to True, any warnings related to the lack of
			compositional data or interaction parameters will be printed.
	"""
	# Note: Should Fe be included in the interactions for this calculation?
	metal_sample = sample_class.Sample(metal_comp)

	user_interactions = interactions.copy()
	interactions = []
	no_comp_interact = []
	no_interaction_param = []
	for element in user_interactions:
		if element in metal_comp and metal_comp[element]>0:
			interactions.append(element)
		else:
			no_comp_interact.append(element)

	if gammaFe == "calculate":
		gammaFe = calc_gamma_Fe_metal(metal_comp=metal_comp, temperature=temperature,
									  interactions=interactions, print_warnings=print_warnings)
	ln_gamma_Fe_metal = np.log(gammaFe)
	
	i = species
	if metal_comp[i]>0:	
		ln_gamma_species_naught_ref = calc_ln_gamma_naught_at_temperature(species=species,
																		  temperature=temperature)
		molfrac = metal_sample.get_composition(units='mol')
		try:
			epsilon_i_i = calc_epsilon_at_temperature(i=species, j=species, temperature=temperature)
		except:
			epsilon_i_i = 0 #TODO should this be 1?

		term1 = 0
		term2 = 0

		for j in interactions:
			skip_this = False
			if j == i:
				pass
			else:
				try:
					epsilon_i_j = calc_epsilon_at_temperature(i=i, j=j, temperature=temperature)
				except:
					skip_this = True
					param_string = str(i) + " and " + str(j)
					no_interaction_param.append(param_string)

				if skip_this == False:
					term1 += epsilon_i_j*molfrac[j]*(1 + ((np.log(1-molfrac[j]))/molfrac[j]) -
														  (1/(1-molfrac[i])))
					term2 += epsilon_i_j*molfrac[j]**2*molfrac[i]*((1/(1-molfrac[i])) +
																   (1/(1-molfrac[j])) +
																   (molfrac[i]/(2*(1-
																   	molfrac[i])**2))-1)

		ln_gamma_i = (ln_gamma_Fe_metal + ln_gamma_species_naught_ref -
					  epsilon_i_i*np.log(1-molfrac[species]) - term1 + term2)

		gamma_i = np.exp(ln_gamma_i)

		if print_warnings == True:
			if len(no_comp_interact) > 0:
				w.warn("No compositional information given for: " +
					', '.join(no_comp_interact) +
					". Interaction with this element not calculated.",RuntimeWarning,stacklevel=2)
			if len(no_interaction_param) > 0:
				w.warn("No interaction parameter known for: " +
					', '.join(no_interaction_param) + ".")

		return gamma_i
	else:
		if print_warnings == True:
			w.warn("No compositional information given for species " + species +
				". Cannot calculate gamma without concentration.",RuntimeWarning,stacklevel=2)
		return 0

def calc_gamma_FeO_silicate(silicate_comp):
	"""
	Returns a value for gammaFeO in the silicate. Parameterization is based on Holdzheid,
	where gammaFeO is taken as a constant value from 1.7-3, dependent only upon MgO content. 

	Parameters
	----------
	silicate_comp: dict
		Dictionary of the composition of the silicate in wt% oxides.

	Returns
	-------
	float
		gammaFeO in the silicate melt
	"""

	if silicate_comp['MgO'] <= 20:
		return 1.7
	else:
		gamma_value =  1.7 + 0.1*(silicate_comp['MgO']-20)
		if gamma_value <= 3:
			return gamma_value
		else:
			return 3

def calc_activity(X, gamma):
	"""
	Returns the value of the activity of any species given the concentration of that species in
	mol fraction, X, and the activity coefficient gamma.

	Parameters
	----------
	X: float
		Concentration of the given species in mol fraction

	gamma: float
		Activity coefficient of the given species

	Returns
	-------
	float
		Activity of the given species
	"""

	return X * gamma

def metal_activity_from_composition(metal_comp, species, temperature, 
									interactions=core.standard_interactions):
	"""
	Returns the activity of the given species in an Fe-rich metal, calculated as X times gamma.

	Parameters
	----------
	metal_comp: dict
		Dictionary of compositional information only for a metal, in terms of wt% elements

	species: string
		Name of desired species for which to calculate the activity. Must match form of elements
		used in MetalSilicate (e.g., 'Fe', 'W', 'Ti')

	temperature: float
		Temperature at which to perform the calculation, in degrees C.

	interactions: list
		OPTIONAL. List of strings of element names. Elements are solutes in a metal Fe liquid
		alloy for which interaction parameters are known and for which the user wishes to calculate
		the effects of interaction within the alloy. Elements need not be infinitely dilute.
		Compositional and temperature ranges for which interaction parameters are known are given
		in the interaction_parameters script within this library.

	Returns
	-------
	float
		Activity of the given species in an Fe-rich metal
	"""
	metal_sample = sample_class.Sample(metal_comp)

	if isinstance(species, str):
		pass
	else:
		raise InputError("Species must be passed as string. You have passed " + str(type(species)))

	molfrac_metal = metal_sample.get_composition(units='mol')
	try:
		X = molfrac_metal[species]
	except:
		raise InputError("No compositional information given for " + species)

	if X == 0:
		raise InputError("No compositional information given for " + species)

	if species == "Fe":
		gamma = calc_gamma_Fe_metal(metal_comp=metal_comp, temperature=temperature, 
									interactions=interactions, print_warnings=False)
	elif species in elements:
		gamma = calc_gamma_solute_metal(metal_comp=metal_comp, species=species, 
										temperature=temperature, gammaFe="calculate",
										interactions=core.standard_interactions,
										print_warnings=False)
	else:
		raise InputError("Species must be one of: " + str(elements))

	if gamma == 0:
		raise GeneralError("Something went wrong. Unable to calculate an" +
						   " activity coefficient gamma for " + species)

	return calc_activity(X=X, gamma=gamma)

def calc_IW(pressure, temperature):
	"""
	Fe-FeO (Iron-Wustite)
	=====================
	Define IW buffer value at P
	
	References
	----------
	Campbell et al. (2009) High-pressure effects on the iron-iron oxide and nickel-nickel oxide
	oxygen fugacity buffers

	Parameters
	----------
	pressure: float
		Pressure in GPa

	temperature: float or numpy array
		Temperature in degrees C

	Returns
	-------
	float or numpy array
		log(fO2)

	Polynomial coefficients
	-----------------------
	log fO2  =  (a0+a1*P) + (b0+b1*P+b2*P^2+b3*P^3)/T
	a0: 6.54106 
	a1: 0.0012324
	b0: -28163.6
	b1: 546.32
	b2: -1.13412
	b3: 0.0019274               
	"""
	#convert T from C to K
	T_K = temperature+273.15

	log_fO2 = ((6.54106+0.0012324*pressure) +
			   (-28163.6+546.32*pressure-1.13412*pressure**2+0.0019274*pressure**3)/T_K)

	return log_fO2

def calc_logfO2_from_IW(pressure, temperature, dIW="calculate"):
	"""
	Calculates the absolute fO2 value (as log(fO2)) based on deltaIW value, pressure, and
	temperature.

	Parameters
	----------
	pressure: float
		Pressure in GPa

	temperature: float
		Temperature in degrees C

	dIW: float
		fO2 in terms of deltaIW

	Returns
	-------
	float
		log(fO2) absolute value
	"""

	IW_at_PT = calc_IW(pressure,temperature)
	log_fO2 = IW_at_PT + dIW

	return log_fO2