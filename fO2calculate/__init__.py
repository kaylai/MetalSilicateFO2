"""
fO2calculate

A tool for calculating non-idea activity coefficients for Fe in 
Fe-rich metal alloy and FeO in co-existing silicate melt. These
values can then be used to calculate fO2 relative to the IW
buffer.
"""

__version__ = "0.0.1"
__author__ = 'Kayla Iacovino'

# ----------------- IMPORTS ----------------- #

import pandas as pd

import fO2calculate.core
import fO2calculate.activitycoefficients
import fO2calculate.interactionparameters

import fO2calculate.batchfile

import fO2calculate.sample_class
from fO2calculate.sample_class import Sample

import fO2calculate.fluid_class
from fO2calculate.fluid_class import Fluid

import fO2calculate.gasplots
from fO2calculate.gasplots import make_plots
from fO2calculate.gasplots import make_comparison_plots

import fO2calculate.calculate 
from fO2calculate.calculate import calculate_dIW

import fO2calculate.batch_calculate
from fO2calculate.batch_calculate import BatchFile
