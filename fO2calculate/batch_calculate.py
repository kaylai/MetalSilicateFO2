import pandas as pd
import warnings as w
import numpy as np

from fO2calculate import core
from fO2calculate import sample_class
from fO2calculate import batchfile
from fO2calculate import calculate

class BatchFile(batchfile.BatchFile):
    """Performs model functions on a batchfile.BatchFile object
    """
    pass

    def calculate_dIW(self, temperature,
                      interactions=core.standard_interactions, **kwargs):
        """ Calculates the fO2 of all samples in a BatchFile in terms of number
        of log units from the iron-wüstite buffer (dIW).

        Parameters
        ----------
        sample: Sample class
            Composition of silicate melt and metal phases as a sample object.

        temperature:    float, int, or str
            Temperature in degrees C. Can be passed as float, in which case
            the passed value is used as the temperature for all samples.
            Alternatively, temperature information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the BatchFile
            object.

        interactions:   list
            OPTIONAL. List of strings of element names. Elements are solutes in
            a metal Fe liquid alloy for which interaction parameters are known
            and for which the user wishes to calculate the effects of
            interaction within the alloy. Elements need not be infinitely
            dilute. Compositional and temperature ranges for which interaction
            parameters are known are given in the interaction_parameters script
            within this library.

        Returns
        -------
        pandas DataFrame
            Original data passed plus newly calculated dIW values are returned.
        """
        batch_data = self.get_data().copy()

        if isinstance(temperature, str):
            file_has_temp = True
            temp_name = temperature
        elif isinstance(temperature, float) or isinstance(temperature, int):
            file_has_temp = False
        else:
            raise core.InputError("temperature must be type str, int, or float")

        dIW_vals = []
        warnings = []
        for index, row in batch_data.iterrows():
            try:
                if file_has_temp:
                    temperature = row[temp_name]

                # Get sample comps as Sample class with defaults
                sample_comp = self.get_sample_composition(index,
                                                          units='wtpt',
                                                          asSampleClass=True,
                                                          how='combined')
                sample_comp.set_default_units(self.default_units)
                sample_comp.set_default_normalization(
                                                    self.default_normalization)

                # Run the calculation
                calc = calculate.calculate_dIW(sample=sample_comp,
                                              temperature=temperature,
                                              interactions=interactions,
                                              **kwargs)
                dIW_vals.append(calc.result)
            except Exception:
                dIW_vals.append(np.nan)
                warnings.append("Calculation Failed.")
        batch_data["dIW_Calculated"] = dIW_vals
        if warnings:
            batch_data["Warnings"] = warnings
        else:
            pass

        return batch_data

