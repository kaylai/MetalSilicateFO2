from mendeleev import element
from mendeleev import get_table
import sys

# ----------SOME UNIVERSAL DEFINITIONS---------- #
OxygenNum = {'Al2O3': 3,'CaO': 1, 'Cl': 0, 'CO2': 2, 'CoO': 1, 'Cr2O3': 3,'CuO': 1, 'F': 0,
             'Fe2O3': 3,'FeO': 1, 'Ga2O3': 3,'H2O': 1, 'K2O': 1, 'MgO': 1, 'MnO': 1, 'MoO2': 2,
             'Na2O': 1, 'Nb2O5': 5,'NiO': 1, 'O': 1, 'P2O5': 5, 'PbO': 1, 'ReO3': 3, 'S': 0,
             'SiO2': 2, 'Ta2O5': 5,'TeO2': 2, 'ThO': 1, 'TiO2': 2, 'UO2': 2, 'V2O3': 3, 'WO3': 3,
             'ZnO': 1}

CationNum = {'Al2O3': 2,'CaO': 1, 'Cl': 1, 'CO2': 1, 'CoO': 1, 'Cr2O3': 2,'CuO': 1, 'F': 1,
             'Fe2O3': 2,'FeO': 1, 'Ga2O3': 2,'H2O': 2, 'K2O': 2, 'MgO': 1, 'MnO': 1, 'MoO2': 1,
             'Na2O': 2, 'Nb2O5': 2,'NiO': 1, 'O': 0, 'P2O5': 2, 'PbO': 1, 'ReO3': 1, 'S': 1,
             'SiO2': 1, 'Ta2O5': 2,'TeO2': 1, 'ThO': 1, 'TiO2': 1, 'UO2': 1, 'V2O3': 2, 'WO3': 1,
             'ZnO': 1}

#Critical parameters cP, cT, o for relevant species
#dict of dicts
critical_params = {'CH4':{  "cT":   304.2,
                            "cP":   73.82,
                            "o":    0.228
                        },
                    'CO':{  "cT":   133.15,
                            "cP":   34.9571,
                            "o":    0.049
                        },
                   'CO2':{  "cT":   304.15,
                            "cP":   73.8659,
                            "o":    0.225
                        },
                   'H2':{   "cT":   33.25,
                            "cP":   12.9696,
                            "o":    -0.218
                        },
                   'H2O':{  "cT":   647.25,
                            "cP":   221.1925,
                            "o":    0.334
                        },
                   'H2S':{  "cT":   373.55,
                            "cP":   90.0779,
                            "o":    0.081
                        },
                   'O2':{   "cT":   154.75,
                            "cP":   50.7638,
                            "o":    0.021
                        },
                   'S2':{   "cT":   208.15,
                            "cP":   72.954,
                            "o":    0.0 #need omega value for S2
                        },
                   'SO2':{  "cT":   430.95,
                            "cP":   77.87295,
                            "o":    0.256
                        }
                    }

"""Names"""
silicate_oxides = ['Al2O3','CaO', 'CO2', 'CoO', 'Cr2O3','CuO', 'Fe2O3','FeO', 'Ga2O3','H2O', 
                   'K2O', 'MgO', 'MnO', 'MoO2','Na2O', 'Nb2O5','NiO', 'P2O5', 'PbO', 'ReO3',
                   'SiO2', 'Ta2O5','TeO2', 'ThO', 'TiO2', 'UO2', 'V2O3', 'WO3', 'ZnO']
metal_elements = ['Al','Ca', 'Cl', 'C', 'Co', 'Cr','Cu', 'F', 'Fe3','Fe', 'Ga','H', 
                  'K', 'Mg', 'Mn', 'Mo','Na', 'Nb','Ni', 'O', 'P', 'Pb', 'Re', 'S', 'Si', 
                  'Ta','Te', 'Th', 'Ti', 'U', 'V', 'W', 'Zn']
oxides = ['Al2O3','CaO', 'Cl', 'CO2', 'CoO', 'Cr2O3','CuO', 'F', 'Fe2O3','FeO', 'Ga2O3','H2O', 
          'K2O', 'MgO', 'MnO', 'MoO2','Na2O', 'Nb2O5','NiO', 'O', 'P2O5', 'PbO', 'ReO3','S',
          'SiO2', 'Ta2O5','TeO2', 'ThO', 'TiO2', 'UO2', 'V2O3', 'WO3', 'ZnO']
elements = ['Al','Ca', 'Cl', 'C', 'Co', 'Cr','Cu', 'F', 'Fe3','Fe', 'Ga','H', 
          'K', 'Mg', 'Mn', 'Mo','Na', 'Nb','Ni', 'O', 'P', 'Pb', 'Re', 'S', 'Si', 
          'Ta','Te', 'Th', 'Ti', 'U', 'V', 'W', 'Zn']
oxides_and_elements = oxides + elements
major_oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'FeO', 'MgO', 'CaO', 'Na2O', 'K2O', 'MnO', 'P2O5']
volatiles = ['CO2', 'H2O', 'F', 'Cl', 'S']
all_volatiles = ['C', 'CO', 'CO2', 'CH4', 'H', 'H2O', 'H2', 'OH', 'F', 'Cl', 'S', 'SO2', 'SO3']
fluid_species_names = ['H2', 'CH4', 'H2O', 'CO', 'O2', 'H2S', 'CO2', 'S2', 'SO2']
fluidMass = [2, 16, 18, 28, 32, 34, 44, 64, 64]
fluid_col_names = ['ic_2', 'ic_16', 'ic_18', 'ic_28', 'ic_32', 'ic_34', 'ic_44', 'ic_64', 'ic_64']

"""Transformations"""
oxides_to_cations = dict(zip(oxides, elements))
cations_to_oxides = dict(zip(elements, oxides))

"""Masses"""
trunc_elements = elements.copy()
trunc_elements.remove('Fe3')
elementMass = {i : element(i).atomic_weight for i in trunc_elements}
elementMass.update({'Fe3': element('Fe').atomic_weight})
oxideMass = {k: CationNum[k]*elementMass[v]+OxygenNum[k]*element('O').atomic_weight for k, v in oxides_to_cations.items()}

"""Other standard values"""
corgne_species = ['Ni', 'Cu', 'Si', 'Mn', 'Cr', 'Ga', 'Nb', 'Ta']
standard_interactions = ['S', 'C', 'O', 'Ni', 'Cu', 'Si', 'Mn', 'Cr', 'Ga', 'Nb', 'Ta']
norris_interactions = ['Fe', 'C', 'O', 'Si', 'P', 'S', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Ga',
                       'Ge', 'Zr', 'Nb', 'Mo', 'Hf', 'Ta', 'W']

def status_bar(percent, sample_name, barLen=20):
    """
    Prints a status bar to the terminal.

    percent: float
        Percent value of progress from 0 to 1

    barLen: int
        Length of bar to print
    """ 
    sys.stdout.write("\r")
    sys.stdout.write("[{:<{}}] {:.0f}%".format("=" * int(barLen * percent), barLen, percent * 100))
    sys.stdout.write("  Working on sample " + str(sample_name))
    if percent == 1.0:
        sys.stdout.write("\n")
    sys.stdout.flush()

# ---------ERROR HANDLING-------- #
class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class GeneralError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
