import pandas as pd
import numpy as np
import warnings as w

from fO2calculate import core

from copy import deepcopy, copy

class Sample(object):
    """Based on the sample_class module of VESIcal.

    The sample class stores compositional information for samples, and contains methods for
    normalization and other compositional calculations. Designed to understand both silicate
    melt data and metal alloy data in a single sample. Silicate melt data must be in terms
    of oxides, and metal alloy data must be in terms of elements.
    """

    def __init__(self, composition, units='wtpt', default_normalization='none',
                 default_units='wtpt', silence_warnings=False):
        """ Initialises the sample class.

        The composition is stored as wtpt. If the composition is provided as wtpt, no
        normalization will be applied. If the composition is supplied as mols, the composition
        will be normalized to 100 wt%.

        Parameters
        ----------
        composition     dict or pandas.Series
            The composition of the sample in the format specified by the composition_type
            parameter. Default is oxides in wtpt.

        units     str
            Specifies the units and type of compositional information passed in the composition
            parameter. Choose from 'wtpt', 'mol'.

        default_normalization:     None or str
            The type of normalization to apply to the data by default. One of:
                - None (no normalization)
                - 'standard' (default): Normalizes an input composition to 100%.
                - 'fixedvolatiles': Normalizes major element oxides to 100 wt%, including
                  volatiles. The volatile wt% will remain fixed, whilst the other major element
                  oxides are reduced proportionally so that the total is 100 wt%.
                - 'additionalvolatiles': Normalises major element oxide wt% to 100%, assuming it
                  is volatile-free. If H2O or CO2 are passed to the function, their un-normalized
                  values will be retained in addition to the normalized non-volatile oxides,
                  summing to >100%.

        default_units     str
            The type of composition to return by default, one of:
            - wtpt (default)
            - mol
        """

        self._silence_warnings = silence_warnings
        composition = deepcopy(composition)

        if isinstance(composition, dict):
            composition = pd.Series(composition, dtype='float64')
        elif isinstance(composition, pd.Series) is False:
            raise core.InputError("The composition must be given as either a dictionary or a "
                                  "pandas Series.")

        if units == 'wtpt':
            self._composition = composition
        elif units == 'mol':
            self._composition = self._mol_to_wtpercent(composition)
        else:
            raise core.InputError("Units must be one of 'wtpt' or 'mol'.")

        # check if ambigous volatiles are passed in, warn user
        ambiguous_voltiles = ['Cl', 'F', 'S']
        if silence_warnings is False:
            for vol in ambiguous_voltiles:
                if vol in composition.index:
                    w.warn("Element " + str(vol) + " was passed as part of sample composition. This" +
                            " will be assumed to be part of the metal composition.")

        self.set_default_normalization(default_normalization)
        self.set_default_units(default_units)

    def set_default_normalization(self, default_normalization):
        """ Set the default type of normalization to use with the get_composition() method.

        Parameters
        ----------
        default_normalization:    str
            The type of normalization to apply to the data. One of:
            - 'none' (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            - 'fixedvolatiles': Normalizes major element oxides to 100 wt%, including volatiles.
              The volatile wt% will remain fixed, whilst the other major element oxides are
              reduced proportionally so that the total is 100 wt%.
            - 'additionalvolatiles': Normalises major element oxide wt% to 100%, assuming it is
              volatile-free. If H2O or CO2 are passed to the function, their un-normalized values
              will be retained in addition to the normalized non-volatile oxides, summing to >100%.
        """
        if default_normalization in ['none', 'standard', 'fixedvolatiles', 'additionalvolatiles']:
            self.default_normalization = default_normalization
        else:
            raise core.InputError("The normalization method must be one of 'none', 'standard', "
                                  "'fixedvolatiles', or 'additionalvolatiles'.")

    def set_default_units(self, default_units):
        """ Set the default units of composition to return when using the get_composition() method.

        Parameters
        ----------
        default_units     str
            The type of composition to return, one of:
            - wtpt (default)
            - mol
        """
        if default_units in ['wtpt', 'mol']:
            self.default_units = default_units
        else:
            raise core.InputError("The units must be one of 'wtpt' or 'mol'.")

    def get_composition(self, species=None, normalization=None, units=None,
                        exclude_volatiles=False, asSampleClass=False, oxide_masses={},
                        how='combined', **kwargs):
        """ Returns the silicate and metal composition in the format requested, normalized as
        requested.

        Parameters
        ----------
        species:    NoneType or str
            The name of the oxide or cation to return the concentration of. If NoneType (default)
            the whole composition will be returned as a pandas.Series. If an oxide is passed, the
            value in wtpt will be returned unless units is set to 'mol', even if the
            default units for the sample object are mol. If an element is passed, the
            concentration will be returned as mol_cations, unless 'mol_singleO' is specified as
            units, even if the default units for the sample object are mol_singleO. Unless
            normalization is specified in the method call, none will be applied.

        normalization:     NoneType or str
            The type of normalization to apply to the data. One of:
            - 'none' (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            - 'fixedvolatiles': Normalizes major element oxides to 100 wt%, including volatiles.
              The volatile wt% will remain fixed, whilst the other major element oxides are reduced
              proportionally so that the total is 100 wt%.
            - 'additionalvolatiles': Normalises major element oxide wt% to 100%, assuming it is
              volatile-free. If H2O or CO2 are passed to the function, their un-normalized values
              will be retained in addition to the normalized non-volatile oxides, summing to >100%.
            If NoneType is passed the default normalization option will be used
            (self.default_normalization).

        units:     NoneType or str
            The units of composition to return, one of:
            - wtpt (default)
            - mol
            If NoneType is passed the default units option will be used (self.default_type).

        exclude_volatiles   bool
            If True, volatiles will be excluded from the returned composition, prior to
            normalization and conversion.

        asSampleClass:  bool
            If True, the sample composition will be returned as a sample class, with default
            options. In this case any normalization instructions will be ignored.

        oxide_masses:  dict
            Specify here any oxide masses that should be changed from the VESIcal default. This
            might be useful for recreating other implementations of models that use slightly
            different molecular masses. The default values in VESIcal are given to 3 dp.

        how:    str
            Specify which composition to return. Either: 'combined' for both metal and silicate
            composition (default); 'metal' for only the metal composition; 'silicate' for only
            the silicate composition. Intended to be used by get_metal_composition() and 
            get_silicate_composition() functions.

        Returns
        -------
        pandas.Series, float, or Sample class
            The sample composition, as specified.
        """
        # Process the oxide_masses variable, if necessary:
        oxideMass = copy(core.oxideMass)
        for ox in oxide_masses:
            if ox in oxideMass:
                oxideMass[ox] = oxide_masses[ox]
            else:
                raise core.InputError("The oxide name provided in oxide_masses is not recognised.")

        # Fetch the default return types if not specified in function call
        if normalization is None and species is None:
            normalization = self.default_normalization
        if units is None and species is None:
            units = self.default_units

        # Check whether to exclude volatiles
        # note that here composition is gotten as wtpt
        if exclude_volatiles:
            composition = self._composition.copy()
            for vol in core.all_volatiles:
                if vol in composition.index:
                    composition = composition.drop(index=vol)
        else:
            composition = self._composition.copy()

        # Check for a species being provided, if so, work out which units to return.
        if isinstance(species, str):
            if species in composition.index:  # if the requested species has a value, proceed
                if species in core.elements or species in core.oxides:
                    pass
                else:
                    raise core.InputError(species + " was not recognised, check spelling, " +
                                          "capitalization and stoichiometry.")
                if normalization is None:
                    normalization = 'none'
            else:
                return 0.0  # if the requested species has no set value, return a float of 0.0
        elif species is not None:
            raise core.InputError("Species must be either a string or a NoneType.")

        # Get the requested type of composition
        if units == 'wtpt':
            converted = composition
        elif units == 'mol':
            converted = self._wtpercent_to_mol(composition)
        else:
            raise core.InputError("The units must be one of 'wtpt' or 'mol'.")

        # Do requested normalization
        if normalization == 'none':
            final = converted
        elif normalization == 'standard':
            final = self._normalize_Standard(converted, units=units)
        elif normalization == 'fixedvolatiles':
            final = self._normalize_FixedVolatiles(converted, units=units)
        elif normalization == 'additionalvolatiles':
            final = self._normalize_AdditionalVolatiles(converted, units=units)
        else:
            raise core.InputError("The normalization method must be one of 'none', 'standard', "
                                  "'fixedvolatiles', or 'additionalvolatiles'.")

        if how is "combined":
            if species is None:
                if asSampleClass is False:
                    return final
                else:
                    return Sample(final, **kwargs)
            elif isinstance(species, str):
                if asSampleClass:
                    w.warn("Cannot return single species as Sample class. Returning as float.",
                           RuntimeWarning, stacklevel=2)
                return final[species]

        elif how is "silicate":
            for i in final.index:
                if i not in core.silicate_oxides:
                    final = final.drop(index=i)
            if asSampleClass is False:
                return final
            else:
                return Sample(final, **kwargs)

        elif how is "metal":
            for i in final.index:
                if i not in core.metal_elements:
                    final = final.drop(index=i)
            if asSampleClass is False:
                return final
            else:
                return Sample(final, **kwargs)

    def get_silicate_composition(self, **kwargs):
        """
        Returns only the silicate composition. Inherits all arguments from get_composition()
        """

        return self.get_composition(how='silicate', **kwargs)

    def get_metal_composition(self, **kwargs):
        """
        Returns only the metal composition. Inherits all arguments from get_composition()
        """

        return self.get_composition(how='metal', **kwargs)

    def change_composition(self, new_composition, units='wtpt', inplace=True):
        """
        Change the concentration of some component of the composition.

        If the units are moles, they are read as moles relative to the present composition,
        i.e. if you wish to double the moles of MgO, if the present content is 0.1 moles,
        you should provide {'MgO':0.2}. The composition will then be re-normalized. If the
        original composition was provided in un-normalized wt%, the unnormalized total will
        be lost.

        Parameters
        ----------
        new_composition:    dict or pandas.Series
            The components to be updated.
        units:      str
            The units of new_composition. Should be one of:
            - wtpt (default)
            - mol
        inplace:    bool
            If True the object will be modified in place. If False, a copy of the Sample
            object will be created, modified, and then returned.

        Returns
        -------
        Sample class
            Modified Sample class.
        """

        # if new_composition is pandas.Series, convert to dict
        if isinstance(new_composition, pd.Series):
            new_composition = dict(new_composition)

        if inplace is False:
            newsample = deepcopy(self)
            return newsample.change_composition(new_composition, units=units)

        if units == 'wtpt':
            for ox in new_composition:
                self._composition[ox] = new_composition[ox]

        elif units == 'mol':
            _comp = self.get_composition(units='mol')
            for ox in new_composition:
                _comp[ox] = new_composition[ox]
            self._composition = self._mol_to_wtpercent(_comp)

        else:
            raise core.InputError("Units must be one of 'wtpt' or 'mol'.")

        return self

    def check_oxide(self, oxide):
        """
        Check whether the sample composition contains the given oxide.

        Parameters
        ----------
        oxide:  str
            Oxide name to check composition for.

        Returns
        -------
        bool
            Whether the composition contains the given oxide, or not.
        """

        if oxide not in core.oxides:
            w.warn("Oxide name not recognised. If it is in your sample, unexpected behaviour "
                   "might occur!",
                   RuntimeWarning, stacklevel=2)
        return oxide in self._composition

    def check_cation(self, cation):
        """
        Check whether the sample composition contains the given cation.

        Parameters
        ----------
        cation:    str
            The element name to check the composition for.

        Returns
        -------
        bool
            Whether the composition contains the given element, or not.
        """

        if cation not in core.elements:
            w.warn("Cation name not recognised. If it is in your sample, unexpected behaviour "
                   "might occur!",
                   RuntimeWarning, stacklevel=2)
        return cation in self.get_composition(units='mol')

    def _normalize_Standard(self, composition, units='wtpt'):
        """
        Normalizes the given composition to 100 wt%, including volatiles. This method
        is intended only to be called by the get_composition() method.

        Parameters
        ----------
        composition:     pandas.Series
            A rock composition with oxide names as keys and concentrations as values.

        units:      str
            The units of composition. Should be one of:
            - wtpt (default)
            - mol

        Returns
        -------
        pandas.Series
            Normalized oxides in wt%.
        """
        ox_comp = self.get_silicate_composition()
        ox_comp = dict(ox_comp)
        elem_comp = self.get_metal_composition()
        elem_comp = dict(elem_comp)

        if units == 'wtpt':
            normed_ox = {k: 100.0 * v / sum(ox_comp.values()) for k, v in ox_comp.items()}
            normed_elem = {k: 100.0 * v / sum(elem_comp.values()) for k, v in elem_comp.items()}
            normed = pd.Series({**normed_ox, **normed_elem})
        elif units == 'mol':
            normed_ox = {k: v / sum(ox_comp.values()) for k, v in ox_comp.items()}
            normed_elem = {k: v / sum(elem_comp.values()) for k, v in elem_comp.items()}
            normed = pd.Series({**normed_ox, **normed_elem})
        else:
            raise core.InputError("Units must be one of 'wtpt' or 'mol'.")

        return normed

    def _normalize_FixedVolatiles(self, composition, units='wtpt'):
        """
        Normalizes major element oxides to 100 wt%, including volatiles. The volatile wt% will
        remain fixed, whilst the other major element oxides are reduced proportionally so that the
        total is 100 wt%.

        Intended to be called only by the get_composition() method.

        Parameters
        ----------
        composition:     pandas Series
            Major element composition

        units:      str
            The units of composition. Should be one of:
            - wtpt (default)
            - mol

        Returns
        -------
        pandas Series
            Normalized major element oxides.
        """
        ox_comp = self.get_silicate_composition()
        ox_comp = dict(ox_comp)
        elem_comp = self.get_metal_composition()
        elem_comp = dict(elem_comp)
        ox_norm = {}
        elem_norm = {}

        ox_vols = 0
        elem_vols = 0
        for vol in core.all_volatiles:
            if vol in ox_comp.keys():
                ox_vols += ox_comp[vol]
            if vol in elem_comp.keys():
                elem_vols += elem_comp[vol]

            for k, v in ox_comp.items():
                if k != vol:
                    ox_norm.update({k: v})

            for k, v in elem_comp.items():
                if k != vol:
                    elem_norm.update({k: v})

        if units is 'wtpt':
            ox_norm = {k: (100-ox_vols) * v / (sum(ox_norm.values())-ox_vols) for k, v in ox_norm.items()}
            elem_norm = {k: (100-
                             elem_vols) * v / (sum(elem_norm.values())-elem_vols) for k, v in elem_norm.items()}
        elif units is 'mol':
            ox_norm = {k: (1-ox_vols) * v / (sum(ox_norm.values())-ox_vols) for k, v in ox_norm.items()}
            elem_norm = {k: (1-
                             elem_vols) * v / (sum(elem_norm.values())-elem_vols) for k, v in elem_norm.items()}
        else:
            raise core.InputError("Units must be one of 'wtpt' or 'mol'.")

        for vol in core.all_volatiles:
            if vol in ox_norm.keys():
                ox_norm.update({vol: ox_comp[vol]})
            if vol in elem_norm.keys():
                elem_norm.update({vol: elem_comp[vol]})

        final = pd.Series({**ox_norm, **elem_norm})

        return final

    def _normalize_AdditionalVolatiles(self, composition, units='wtpt'):
        """
        Normalises major element oxide wt% to 100%, assuming it is volatile-free. If volatiles
        are passed to the function, their un-normalized values will be retained in addition to the
        normalized non-volatile oxides, summing to >100%.

        Intended to be called only by the get_composition() method.

        Parameters
        ----------
        composition:     pandas.Series
            Major element composition

        units:      str
            The units of composition. Should be one of:
            - wtpt (default)
            - mol

        Returns
        -------
        pandas.Series
            Normalized major element oxides.
        """
        ox_comp = self.get_silicate_composition()
        ox_comp = dict(ox_comp)
        elem_comp = self.get_metal_composition()
        elem_comp = dict(elem_comp)
        ox_norm = {}
        elem_norm = {}

        ox_vols = 0
        elem_vols = 0
        for vol in core.all_volatiles:
            if vol in ox_comp.keys():
                ox_vols += ox_comp[vol]
            if vol in elem_comp.keys():
                elem_vols += elem_comp[vol]

            for k, v in ox_comp.items():
                if k != vol:
                    ox_norm.update({k: v})

            for k, v in elem_comp.items():
                if k != vol:
                    elem_norm.update({k: v})

        if units is 'wtpt':
            ox_norm = {k: 100 * v / (sum(ox_norm.values()) - ox_vols) for k, v in ox_norm.items()}
            elem_norm = {k: 100 * v / (sum(elem_norm.values()) - 
                                       elem_vols) for k, v in elem_norm.items()}
        elif units is 'mol':
            ox_norm = {k: 100 * v / (sum(ox_norm.values()) - ox_vols) for k, v in ox_norm.items()}
            elem_norm = {k: 100 * v / (sum(elem_norm.values()) -
                                       elem_vols) for k, v in elem_norm.items()}
        else:
            raise core.InputError("Units must be one of 'wtpt' or 'mol'.")

        for vol in core.all_volatiles:
            if vol in ox_norm.keys():
                ox_norm.update({vol: ox_comp[vol]})
            if vol in elem_norm.keys():
                elem_norm.update({vol: elem_comp[vol]})

        final = pd.Series({**ox_norm, **elem_norm})

        return final

    def _wtpercent_to_mol(self, composition, oxideMass=core.oxideMass, 
                          elementMass=core.elementMass):
        """
        Converts a wt% oxide composition to mol oxides, normalised to 1 mol.

        Intended to be called only by the get_composition() method.

        Parameters
        ----------
        composition:    pandas.Series
            Major element composition in wt%

        oxideMass:  dict
            The molar mass of the oxides. Default is the VESIcal default
            molar masses.

        elementMass: dict
            The molar mass of the elements. Default is the VESIcal default
            element oxides.

        Returns
        -------
        pandas.Series
            Molar proportions of major element oxides, normalised to 1.
        """
        mols_ox = {}
        mols_elem = {}
        comp = composition.copy()
        specieslist = list(comp.index)

        for spec in specieslist:
            if spec in core.silicate_oxides:
                mols_ox[spec] = comp[spec]/oxideMass[spec]
            if spec in core.metal_elements:
                mols_elem[spec] = comp[spec]/elementMass[spec]

        # do normalization
        ox_sum = sum(mols_ox.values())
        ox_final = {k: v / ox_sum for k, v in mols_ox.items()}
        elem_sum = sum(mols_elem.values())
        elem_final = {k: v / elem_sum for k, v in mols_elem.items()}

        combined = {**ox_final, **elem_final}
        combined = pd.Series(combined)

        return combined

    def _mol_to_wtpercent(self, composition, oxideMass=core.oxideMass,
                          elementMass=core.elementMass):
        """
        Converts mol to wt%. Returned composition is normalized to 100 wt%.

        Parameters
        ----------
        composition:     pandas.Series
            mol fraction 

        oxideMass:  dict
            The molar mass of the oxides. Default is the VESIcal default 
            molar masses.

        elementMass: dict
            The molar mass of the elements. Default is the VESIcal default
            element oxides.

        Returns
        -------
        pandas.Series
            wt% oxides normalized to 100 wt%.
        """
        wtpt_ox = {}
        wtpt_elem = {}
        comp = composition.copy()
        specieslist = list(comp.index)

        for spec in specieslist:
            if spec in core.silicate_oxides:
                wtpt_ox[spec] = comp[spec]*oxideMass[spec]
            if spec in core.metal_elements:
                wtpt_elem[spec] = comp[spec]*elementMass[spec]

        # do normalization
        ox_sum = sum(wtpt_ox.values())
        ox_final = {k: 100*v / ox_sum for k, v in wtpt_ox.items()}
        elem_sum = sum(wtpt_elem.values())
        elem_final = {k: 100*v / elem_sum for k, v in wtpt_elem.items()}

        combined = {**ox_final, **elem_final}
        combined = pd.Series(combined)

        return combined