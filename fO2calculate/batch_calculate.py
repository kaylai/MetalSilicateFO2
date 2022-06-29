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

    def calculate_dIW(self, temperature, pressure=None,
                      returnlogfO2="automatic",
                      interactions=core.standard_interactions, **kwargs):
        """ Calculates the fO2 of all samples in a BatchFile in terms of number
        of log units from the iron-wüstite buffer (dIW).

        Parameters
        ----------
        temperature:    float, int, or str
            Temperature in degrees C. Can be passed as float, in which case
            the passed value is used as the temperature for all samples.
            Alternatively, temperature information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the BatchFile
            object.

        pressure: float, int, or str
            Pressure in bars. Only required if it is desired that the fO2 be
            returned both as dIW and as logfO2. Can be passed as float, in which
            case the passed value is used as the pressure for all samples.
            Alternatively, pressure information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the BatchFile
            object.

        returnlogfO2: str or bool
            If set to True, function will return fO2 values both in terms of
            dIW and as logfO2. If True, a pressure must also be passed. If
            set to "automatic" (the default), will be set to True if a pressure
            is passed or set to False if no pressure is passed.

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

        if isinstance(pressure, str):
            file_has_pressure = True
            pressure_name = pressure
            if returnlogfO2 == "automatic":
                returnlogfO2 = True
        elif isinstance(pressure, float) or isinstance(pressure, int):
            file_has_pressure = False
            if returnlogfO2 == "automatic":
                returnlogfO2 = True
        elif pressure is None:
            file_has_pressure = False
            if returnlogfO2 == "automatic":
                returnlogfO2 = False
            else:
                raise core.InputError("If returnlogfO2 is True, you must " +
                                      "also pass a pressure.")
        else:
            raise core.InputError("pressure must be type None, str, int, or" +
                                  " float")

        dIW_vals = []
        logfO2_vals = []
        warnings = []
        iterno = 0
        for index, row in batch_data.iterrows():
            iterno += 1
            percent = iterno/len(batch_data.index)
            batchfile.status_bar.status_bar(percent, index)

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

            try: # Also return logfO2 if called for
                if returnlogfO2 is True:
                    if file_has_pressure:
                        pressure = row[pressure_name]

                    dIW = calc.result
                    logfO2 = calculate.calc_logfO2_from_IW(pressure=pressure,
                                                           temperature=temperature,
                                                           dIW=dIW)
                    logfO2_vals.append(logfO2)
            except Exception:
                logfO2_vals.append(np.nan)

        batch_data["dIW_Calculated"] = dIW_vals
        if returnlogfO2 is True:
            batch_data["logfO2_dIW"] = logfO2_vals
        if warnings:
            batch_data["Warnings"] = warnings
        else:
            pass

        if file_has_temp is False:
            batch_data["Temperature_C_Modeled"] = [temperature]*len(batch_data.index)
        if file_has_pressure is False and pressure is not None:
            batch_data["Pressure_Modeled"] = [pressure]*len(batch_data.index)

        return batch_data

    def calculate_dSiSiO2(self, temperature, aSiO2, pressure=None,
                          returnlogfO2="automatic",
                          returndIW="automatic",
                          interactions=core.standard_interactions, **kwargs):
        """ Calculates the fO2 of all samples in a BatchFile in terms of number
        of log units from the Si-SiO2 buffer.

        Parameters
        ----------
        temperature:    float, int, or str
            Temperature in degrees C. Can be passed as float, in which case
            the passed value is used as the temperature for all samples.
            Alternatively, temperature information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the BatchFile
            object.

        aSiO2: float, int, or str
            Activity of SiO2 in the melt. Can be passed as float, in which case
            the passed value is used as the aSiO2 for all samples.
            Alternatively, aSiO2 information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the BatchFile
            object.

        pressure: float, int, or str
            Pressure in bars. Only required if it is desired that the fO2 be
            returned both as dSiSiO2 and as logfO2. Can be passed as float, in
            which case the passed value is used as the pressure for all samples.
            Alternatively, pressure information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the BatchFile
            object.

        returnlogfO2: str or bool
            If set to True, function will return fO2 values both in terms of
            dSiSiO2 and as logfO2. If True, a pressure must also be passed. If
            set to "automatic" (the default), will be set to True if a pressure
            is passed or set to False if no pressure is passed.

        returndIW: str or bool
            If set to True, function will return fO2 values both in terms of
            dSiSiO2 and as dIW. If True, a pressure must also be passed. If
            set to "automatic" (the default), will be set to True if a pressure
            is passed or set to False if no pressure is passed.

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
            Original data passed plus calculated dSiSiO2 values are returned.
        """
        batch_data = self.get_data().copy()

        if isinstance(temperature, str):
            file_has_temp = True
            temp_name = temperature
        elif isinstance(temperature, float) or isinstance(temperature, int):
            file_has_temp = False
        else:
            raise core.InputError("temperature must be type str, int, or float")

        if isinstance(aSiO2, str):
            file_has_aSiO2 = True
            aSiO2_name = aSiO2
        elif isinstance(aSiO2, float) or isinstance(aSiO2, int):
            file_has_aSiO2 = False
        else:
            raise core.InputError("aSiO2 must be type str, int, or float")

        if isinstance(pressure, str):
            file_has_pressure = True
            pressure_name = pressure
            if returnlogfO2 == "automatic":
                returnlogfO2 = True
            if returndIW == "automatic":
                returndIW = True
        elif isinstance(pressure, float) or isinstance(pressure, int):
            file_has_pressure = False
            if returnlogfO2 == "automatic":
                returnlogfO2 = True
            if returndIW == "automatic":
                returndIW = True
        elif pressure is None:
            file_has_pressure = False
            if returnlogfO2 == "automatic":
                returnlogfO2 = False
            if returndIW == "automatic":
                returndIW = False
            else:
                raise core.InputError("If returnlogfO2 or returndIW is True," +
                                      " you must also pass a pressure.")
        else:
            raise core.InputError("pressure must be type None, str, int, or" +
                                  " float")

        dSiSiO2_vals = []
        logfO2_vals = []
        dIW_vals = []
        warnings = []
        iterno = 0
        for index, row in batch_data.iterrows():
            iterno += 1
            percent = iterno/len(batch_data.index)
            batchfile.status_bar.status_bar(percent, index)

            try:
                if file_has_temp:
                    temperature = row[temp_name]

                if file_has_aSiO2:
                    aSiO2 = row[aSiO2_name]

                # Get sample comps as Sample class with defaults
                sample_comp = self.get_sample_composition(index,
                                                          units='wtpt',
                                                          asSampleClass=True,
                                                          how='combined')
                sample_comp.set_default_units(self.default_units)
                sample_comp.set_default_normalization(
                                                    self.default_normalization)

                # Run the calculation
                calc = calculate.calculate_dSiSiO2(sample=sample_comp,
                                                   temperature=temperature,
                                                   aSiO2=aSiO2,
                                                   interactions=interactions,
                                                   **kwargs)
                dSiSiO2_vals.append(calc.result)
            except Exception:
                dSiSiO2_vals.append(np.nan)
                warnings.append("Calculation Failed.")

            # Also return logfO2/dIW if called for
            if returndIW is True or returnlogfO2 is True:
                try: 
                    if file_has_pressure:
                        pressure = row[pressure_name]

                    dSiSiO2 = calc.result
                    logfO2 = calculate.calc_logfO2_from_SiSiO2(pressure=pressure,
                                                        temperature=temperature,
                                                        dSiSiO2=dSiSiO2)
                    logfO2_vals.append(logfO2)
                    
                    if returndIW is True:
                        dIW = calculate.calc_dIW_from_logfO2(pressure=pressure,
                                                        temperature=temperature,
                                                        logfO2=logfO2)
                        dIW_vals.append(dIW)
                except Exception:
                    logfO2_vals.append(np.nan)
                    dIW_vals.append(np.nan)


        batch_data["dSiSiO2_Calculated"] = dSiSiO2_vals
        if returnlogfO2 is True:
            batch_data["logfO2_dSiSiO2"] = logfO2_vals
        if returndIW is True:
            batch_data["dIW_SiSiO2"] = dIW_vals
        if warnings:
            batch_data["Warnings"] = warnings
        else:
            pass

        if file_has_temp is False:
            batch_data["Temperature_C_Modeled"] = [temperature]*len(batch_data.index)
        if file_has_aSiO2 is False:
            batch_data["aSiO2_Modeled"] = [aSiO2]*len(batch_data.index)
        if file_has_pressure is False and pressure is not None:
            batch_data["Pressure_Modeled"] = [pressure]*len(batch_data.index)

        return batch_data

    def calculate_gamma_Fe_metal(self, temperature,
                                 interactions=core.standard_interactions,
                                 print_warnings=False, **kwargs):
        """ Calculates the activity coefficient, gamma, for iron. Interaction
        parameters epsilon are computed for all elements passed to interactions,
        so long as interaction parameter values are known.

        Parameters
        ----------
        temperature:    float, int, or str
            Temperature in degrees C. Can be passed as float, in which case
            the passed value is used as the temperature for all samples.
            Alternatively, temperature information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the BatchFile
            object.

        interactions: list
            OPTIONAL. List of strings of element names. Elements are solutes in
            a metal Fe liquid alloy for which interaction parameters are known
            and for which the user wishes to calculate the effects of
            interaction within the alloy. Elements need not be infinitely
            dilute. Compositional and temperature ranges for which interaction
            parameters are known are given in the interaction_parameters script
            within this library.
        
        print_warnings: bool
                OPTIONAL. Default is True. If set to True, any warnings related
                to the lack of compositional data or interaction parameters will
                be printed.

        Returns
        -------
        pandas DataFrame
            Original data passed plus newly calculated values are returned.
        """
        batch_data = self.get_data().copy()

        if isinstance(temperature, str):
            file_has_temp = True
            temp_name = temperature
        elif isinstance(temperature, float) or isinstance(temperature, int):
            file_has_temp = False
        else:
            raise core.InputError("temperature must be type str, int, or float")

        gammaFe_vals = []
        warnings = []
        iterno = 0
        for index, row in batch_data.iterrows():
            iterno += 1
            percent = iterno/len(batch_data.index)
            batchfile.status_bar.status_bar(percent, index)

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
                calc = calculate.calculate_gamma_Fe_metal(sample=sample_comp,
                                                    temperature=temperature,
                                                    interactions=interactions,
                                                    **kwargs)
                gammaFe_vals.append(calc.result)
            except Exception:
                gammaFe_vals.append(np.nan)
                warnings.append("Calculation Failed.")
        batch_data["gamma_Fe_metal"] = gammaFe_vals
        if warnings:
            batch_data["Warnings"] = warnings
        else:
            pass

        if file_has_temp is False:
            batch_data["Temperature_C_Modeled"] = [temperature]*len(batch_data.index)

        return batch_data

    def calculate_gamma_solute_metal(self, temperature,
                                     species=core.standard_interactions,
                                     interactions=core.standard_interactions,
                                     gammaFe="calculate",
                                     print_warnings=False,
                                     **kwargs):
        """
        Calculates the activity coefficient, gamma, for multiple solutes in an
        Fe-rich metal alloy. Interaction parameters epsilon are computed for all
        elements passed to interactions, so long as interaction parameter values
        are known.

        Parameters
        ----------
        sample: Sample class
            Composition of silicate melt and metal phases as a sample object.

        species: str
            String of the name of the element for which to calculate gamma

        temperature: float, int, or str
            Temperature in degrees C. Can be passed as float, in which case
            the passed value is used as the temperature for all samples.
            Alternatively, temperature information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the BatchFile
            object.

        gammaFe: float, int, or str
            OPTIONAL. Default is "calculate" in which case the gammaFe value
            will be calculated here. If the gammaFe value is already known, it
            can be passed here to avoid having to calculate it again. Can be
            passed as float or int, in which case the passed value is used as
            the gammaFe for all samples. Alternatively, gammaFe values for each
            individual sample may already be present in the BatchFile object.
            If so, pass the str value corresponding to the column title in the
            BatchFile object.

        elements: list
            OPTIONAL. List of elements for which to calculate gamma values, if
            that list is different than the list of interactions. Default value
            is None, in which case the elements list = interactions.

        interactions: list
            OPTIONAL. List of strings of element names. Elements are solutes in
            a metal Fe liquid alloy for which interaction parameters are known
            and for which the user wishes to calculate the effects of
            interaction within the alloy. Elements need not be infinitely
            dilute. Compositional and temperature ranges for which interaction
            parameters are known are given in the interaction_parameters script
            within this library.
        
        print_warnings: bool
                OPTIONAL. Default is True. If set to True, any warnings related
                to the lack of compositional data or interaction parameters will
                be printed.

        Returns
        -------
        pandas DataFrame
            Original data passed plus newly calculated values are returned.
        """
        batch_data = self.get_data().copy()

        if isinstance(temperature, str):
            file_has_temp = True
            temp_name = temperature
        elif isinstance(temperature, float) or isinstance(temperature, int):
            file_has_temp = False
        else:
            raise core.InputError("temperature must be type str, int, or float")

        if isinstance(gammaFe, str):
            if gammaFe == "calculate":
                file_has_gammaFe = False
            else:
                file_has_gammaFe = True
                gammaFe_name = gammaFe
        elif isinstance(gammaFe, float) or isinstance(gammaFe, int):
            file_has_gammaFe = False
        else:
            raise core.InputError("If passed, gammaFe argument must be" +
                                  " 'calculate' or type float, int, or str." +
                                  " Not sure what you want here? Don't pass" +
                                  " this argument.")

        if isinstance(species, str):
            species = [species]

        specno = 0
        iterno = 0
        for spec in species:
            specno += 1
            gamma_solute_vals = []
            for index, row in batch_data.iterrows():
                iterno += 1
                percno = iterno*specno
                length = len(batch_data.index) * len(species)
                percent = iterno/length
                batchfile.status_bar.status_bar(percent, index)

                if file_has_temp is True:
                    temperature = row[temp_name]

                if file_has_gammaFe is True:
                    gammaFe = row[gammaFe_name]

                sample_comp = self.get_sample_composition(index,
                                                          units='wtpt',
                                                          asSampleClass=True,
                                                          how='combined',
                                                          silence_warnings=True)
                # run the calculation
                calc = calculate.calculate_gamma_solute_metal(sample=sample_comp, 
                                                species=spec, 
                                                temperature=temperature,
                                                gammaFe=gammaFe,
                                                interactions=interactions,
                                                print_warnings=print_warnings,
                                                **kwargs)
                gamma_solute_vals.append(calc.result)
            column_head_title = "gamma_" + spec + "_metal"
            batch_data[column_head_title] = gamma_solute_vals

        if file_has_temp is False:
            batch_data["Temperature_C_Modeled"] = [temperature]*len(batch_data.index)

        return batch_data


