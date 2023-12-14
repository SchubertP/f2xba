"""Implementation of RbaMedium class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import pandas as pd

MAX_MEDIUM_CONC = 10.0


class RbaMedium:

    def __init__(self):
        self.concentrations = {}

    def import_xml(self, model_dir):
        # actually an import from tsv file
        file_name = os.path.join(model_dir, 'medium.tsv')
        if os.path.exists(model_dir) is True:
            df = pd.read_table(file_name, usecols=['Metabolite', 'Concentration'], index_col='Metabolite')
            self.concentrations = df.to_dict()['Concentration']
        else:
            print(f'{file_name} not found!')

    def from_df(self, m_dict):
        if 'medium' in m_dict:
            self.concentrations = m_dict['medium'].to_dict()['Concentration']
        else:
            print(f'medium not imported!')

    def from_xba(self, general_params, xba_model):
        """Configure Medium based xba_model reactions.

        Medium defines the nutrient environment for the model simulation.
        In FBA based simulations, exchange reactions import metabolites into the model.
        Import of a specific metabolite is allowed when the corresponding
        exchange reaction is unblocked in the updake direction (usually a negative lower flux bound).

        RBA modelling does not use exchange reactions, instead if defines positive
        metabolite concentrations for the selected environment.
        These concentrations are configured in the RBA medium component.

        Therefore, potential nutrients are defined based on available exchange reactions
        in the metabolic model (FBA model).
        An initial medium is configured based on the flux bounds of exchange reactions.

        Metabolites in the RBA do not contain compartment postfixes.
            i.e. metabolite M_glc__D instead of species id M_glc__D_e

        In RBA, species in the external compartment (those created via exchange reactions in FBA)
        are configured with boundaryCondition set to True, i.e. their concentration is not
        affected during simulation
        Usually ,using Outer Membrane Transporters, external species are transported
        from external compartment to periplasm. These OM transporters will have an associated
        protein cost linear with the total import rate.

        In E. coli metabolite concentration in external medium and periplasm are close to identical.
        RBA then implement Michaelis Menten kinetics for importing species from periplasm to cytosol,
        using as species concentration the respective concentration of the medium metabolite.

        Maximum medium concentration from rba_params['general']['medium_max_conc']['value']

        :param general_params: general RBA parameters loaded from file
        :type general_params: dict
        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: Class XbaModel
        """
        medium_max_conc = general_params.get('medium_max_conc', MAX_MEDIUM_CONC)

        nutrients = 0
        for rid, r in xba_model.reactions.items():
            if r.kind == 'exchange':
                if len(r.reactants) == 1:
                    sid = list(r.reactants)[0]
                    flux_bnd = xba_model.parameters[r.fbc_lower_bound].value
                else:
                    sid = list(r.products)[0]
                    flux_bnd = xba_model.parameters[r.fbc_upper_bound].value
                mid, cid = sid.rsplit('_', 1)
                assert (cid == xba_model.external_compartment)
                if abs(flux_bnd) > 0.0:
                    nutrients += 1
                self.concentrations[mid] = min(abs(flux_bnd), medium_max_conc)
        print(f'{len(self.concentrations):4d} medium metabolites ({nutrients} > 0.0 mmol/l)')

    def export_xml(self, model_dir):
        # actually an export to tsv file
        file_name = os.path.join(model_dir, 'medium.tsv')
        df = self.to_df()['medium']
        df.to_csv(file_name, sep='\t')

    def to_df(self):
        df = pd.DataFrame(self.concentrations.values(), index=list(self.concentrations),
                          columns=['Concentration'])
        df.index.name = 'Metabolite'
        return {'medium': df}
