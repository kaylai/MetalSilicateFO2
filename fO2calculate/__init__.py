"""
fO2calculate

A tool for calculating non-idea activity coefficients for Fe in 
Fe-rich metal alloy and FeO in co-existing silicate melt. These
values can then be used to calculate fO2 relative to the IW
buffer.
"""

__version__ = "0.1.0"
__author__ = 'Kayla Iacovino'

# ----------------- IMPORTS ----------------- #

import pandas as pd

import fO2calculate.core
import fO2calculate.activitycoefficients
import fO2calculate.interactionparameters

import fO2calculate.batchfile

import fO2calculate.sample_class
from fO2calculate.sample_class import Sample

import fO2calculate.calculate 
from fO2calculate.calculate import calculate_dIW
from fO2calculate.calculate import calculate_gamma_Fe_metal
from fO2calculate.calculate import calculate_gamma_solute_metal
from fO2calculate.calculate import calculate_dSiSiO2

import fO2calculate.batch_calculate
from fO2calculate.batch_calculate import BatchFile
