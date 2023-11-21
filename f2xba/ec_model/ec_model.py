"""Implementation of EcModel class.


Peter Schubert, HHU Duesseldorf, October 2023
"""

import re
import numpy as np
import pandas as pd


class EcModel:
    """Class EcModel

    Build upon XbaModel (genome-scale metabolic model), create
    an enzyme constraint model (GECKO, ccFBA, MOMENT, MOMENTmr),
    by modifying the XbaModel.

    .. code-block:: python

        xba_model = XbaModel('iJO1366.xml')
        xba_model.configure('iJO1366_xba_GECKO_parameters.xlsx')

        ec_model = EcModel(xba_model)
        ec_model.configure(('iJO1366_ec_GECKO_parameters.xlsx')
        if ec_model.validate():
            ec_model.export('iJO1366_GECKO_batch.xml')

    """

    def __init__(self, xba_model):
        """Instantiate EcModel

        XbaModel created from a genome-scale metabolic model and
        adequately configured

        :param xba_model: XbaModel
        :type xba_model: Class f2xba.XbaModel
        """
        self.model = xba_model

    def configure(self, ecm_params_fname):
        """Convert the xba model to an enzyme constraint model using parameters in fname.

        to be executed after xba model has been configured
        This method adds EC-model specific parametrizations to the xba model

        Excel spreadsheet document requires sheet 'general' with configuration data for
        'ecm_type', 'p_total', 'avg_enz_sat', ('p_abundance_fname' or 'pm2totpm')

        Notes:
        - GECKO considers isoenzymes, promiscuous enzymes and stoichiometry
          of enzymes/complexes.
        - MOMENTmr is similar to GECKO, except that stoichiometry of enzymes and enzyme
          complexes is fixed at 1 for participating proteins
        - ccFBA is similar to GECKO, with the modification that for reactions catalyzed
          by isoenzymes only one isoreaction with the lease expensive MW/kcat isoenzymes
          is retained in the model, the other isoreactions are removed
          (i.e. loss of isoenzyme redundancy)
        - MOMENT model is similar to MOMENTmr, with the modification the for promiscuous
          enzymes only the most expensive reaction is relevant, all other reactions
          catalyzed by same promiscuous enzyme have no added protein cost.
          Protein species have to be added, before making reactions irreversible

        :param ecm_params_fname: file name for configuration parameters
        :type ecm_params_fname: str
        :return: success/failure of model configuration
        :rtype: bool
        """
        sheet = 'general'
        with pd.ExcelFile(ecm_params_fname) as xlsx:
            ecm_params = pd.read_excel(xlsx, sheet_name=sheet, index_col=0)['value'].to_dict()
            print(f'EC model configuration parameters loaded from {ecm_params_fname}')

        ecm_type = ecm_params['ecm_type']
        p_total = ecm_params['p_total']
        avg_enz_sat = ecm_params['avg_enz_sat']

        if 'p_abundance_fname' in ecm_params:
            pax_db_fname = ecm_params['p_abundance_fname']
            pm2totpm = self.model.get_modelled_protein_mass_fraction(pax_db_fname)
            print(f'{pm2totpm:.3f} g_protein_in_model / g_total_protein, based on {pax_db_fname}')
        else:
            assert('pm2totpm' in ecm_params)
            pm2totpm = ecm_params['pm2totpm']

        protein_pool = p_total * pm2totpm * avg_enz_sat
        print(f'{protein_pool:.3f} g_protein/gDW, assuming avg enzyme saturation of {avg_enz_sat:.2f}')

        # create enzyme constraint model
        if ecm_type in ['GECKO', 'MOMENTmr']:
            self.model.add_isoenzyme_reactions(create_arm=True)
            self.model.make_irreversible()
            self._add_gecko_protein_species()
        elif ecm_type == 'MOMENT':
            self.model.add_isoenzyme_reactions(create_arm=True)
            self._add_moment_protein_species()
            self.model.make_irreversible()
        elif ecm_type == 'ccFBA':
            self._select_ccfba_isoenzyme()
            self.model.make_irreversible()
            self.model.remove_unused_gps()
            self._add_gecko_protein_species()
        else:
            print(f'model not constructed, invalid ecm_type {ecm_type}')
            return False

        self._add_total_protein_constraint(protein_pool)
        self.model.reaction_enzyme_coupling()

        # remove model components that are not used in the model
        self.model.clean()

        # add some parameter values for reference (after model.clean() so they are not removed)
        self.model.add_parameter('frac_prot_totprot', pm2totpm, 'dimensionless')
        self.model.add_parameter('frac_totprot_cdw', p_total, 'dimensionless')
        self.model.add_parameter('frac_enzyme_sat', avg_enz_sat, 'dimensionless')

        # modify some model attributs and create L3V2 SBML model
        self.model.model_attrs['id'] += f'_{ecm_type}'
        self.model.model_attrs['name'] = f'{ecm_type} model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        self.model.print_size()
        return True

    def _add_gecko_protein_species(self):
        """add protein species to the GECKO/ccFBA/MOMENTmr model.

        Note: Protein compartments in XBA model are derived from
        reaction compartments, which in turn is a concatenation of
        species compartments being substrates to the reaction.

        Here we only need one of the species compartments.

        add protein species to the model
        add protein drain reactions to the model
        """
        protein_sids = {}
        for uid, p in self.model.proteins.items():
            prot_sid = f'M_prot_{uid}_{p.cid}'
            p.link_sid(prot_sid)
            prot_metaid = f'meta_prot_{uid}'
            prot_annot = f'bqbiol:is, uniprot/{uid}'
            protein_sids[prot_sid] = [p.name, p.cid, False, False, False, prot_metaid, prot_annot]

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition',
                'constant', 'metaid', 'miriamAnnotation']
        df_add_species = pd.DataFrame(protein_sids.values(), index=list(protein_sids), columns=cols)
        print(f'{len(df_add_species):4d} protein constraints to add')
        self.model.add_species(df_add_species)

    def _add_moment_protein_species(self):
        """Add protein species to MOMENT model.

        Execute before makeing model irreversilbe.

        In MOMENT promiscuous enzymes can catalyze several reactions in parallel
        To implement this constraint, we add for each catalzyed (iso)reaction
        a specific protein species

        add protein species to the model
        """
        # add dedicated protein species for each (iso)reaction
        protein_sids = {}
        for rid, r in self.model.reactions.items():
            assert(len(r.enzymes) <= 1)
            if len(r.enzymes) == 1:
                enz = self.model.enzymes[r.enzymes[0]]
                for locus in enz.composition:
                    uid = self.model.locus2uid[locus]
                    p = self.model.proteins[uid]
                    prot_sid = f'M_prot_{uid}_{rid}_{p.cid}'
                    p.link_sid(prot_sid)
                    prot_metaid = f'meta_prot_{uid}_{rid}'
                    prot_annot = f'bqbiol:is, uniprot/{uid}'
                    protein_sids[prot_sid] = [p.name, p.cid, False, False, False, prot_metaid, prot_annot]

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition',
                'constant', 'metaid', 'miriamAnnotation']
        df_add_species = pd.DataFrame(protein_sids.values(), index=list(protein_sids), columns=cols)
        print(f'{len(df_add_species):4d} moment protein constraints to add')
        self.model.add_species(df_add_species)

    def _select_ccfba_isoenzyme(self):
        """Select least cost enzyme for reactions catalyzed by several isoenzymes.

        Other isoenzymes get deleted and gpa updated.

        Selection is based on molecular weight and turnover number.
        Select enzyme with minimal MW/kcat.

        Different isoenzymes might be selected for forward/reverse direction.
        Therefore, reversible reactions with isoenzymes get split in fwd/bwd reaction.
        """
        # Note model.reactions get updated within the loop, therefore fix rids initially
        rids = list(self.model.reactions)
        for rid in rids:
            r = self.model.reactions[rid]
            if len(r.enzymes) >= 2:
                directional_rs = [r]
                if r.reversible:
                    # split reaction, i.e. create a new irreversible reverse reaction
                    rev_r = self.model.split_reversible_reaction(r)
                    directional_rs.append(rev_r)

                for dir_r in directional_rs:
                    idx = np.argmin([self.model.enzymes[e].mw / kcat for e, kcat in zip(dir_r.enzymes, dir_r.kcatf)])
                    eid = dir_r.enzymes[idx]
                    gpa = ' and '.join(sorted([self.model.uid2gp[self.model.locus2uid[locus]]
                                               for locus in self.model.enzymes[eid].composition]))
                    if ' and ' in gpa:
                        gpa = '(' + gpa + ')'
                    dir_r.enzymes = [dir_r.enzymes[idx]]
                    dir_r.kcatf = [dir_r.kcatf[idx]]
                    dir_r.gene_product_assoc = gpa

    def _add_total_protein_constraint(self, total_protein):
        """add total protein constraint to the enzyme constraint model.

        add protein pool species to the model
        add protein drain reactions to the model
        add protein pool exchange reation to the model

        :param total_protein: total active protein mass balance in g/gDW
        :type total_protein: float
        """
        # add total protein pool species in default compartment
        pool_cid = sorted([(count, cid) for cid, count in self.model.get_used_cids().items()])[-1][1]
        pool_sid = f'M_prot_pool_{pool_cid}'
        pool_name = 'total protein pool'
        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame([[pool_name, pool_cid, False, False, False]],
                                      index=[pool_sid], columns=cols)
        self.model.add_species(df_add_species)

        # add total protein pool exchange reaction
        protein_vars = {}
        draw_rid = f'R_conc_active_protein'
        draw_name = f'total concentration active protein'
        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'mmol_per_gDW', 'conc_0_bound')
        ub_pid = self.model.get_fbc_bnd_pid(total_protein, 'mmol_per_gDW', 'conc_active_protein')
        protein_vars[draw_rid] = [draw_name, f' => {pool_sid}', lb_pid, ub_pid, None]

        # add protein drain reactions - supporting MOMENT with specific protein split in protein species per reaction
        ub_pid = self.model.get_fbc_bnd_pid(1000.0, 'mmol_per_gDW', 'conc_default_ub')
        for uid, p in self.model.proteins.items():
            draw_rid = f'R_conc_prot_{uid}'
            draw_name = f'conc_prot_{uid}'
            products = ' + '.join([sid for sid in p.linked_sids])
            reaction_string = f'{p.mw / 1000.0} {pool_sid} => {products}'
            protein_vars[draw_rid] = [draw_name, reaction_string, lb_pid, ub_pid, self.model.uid2gp[uid]]

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'fbcGeneProdAssoc']
        df_add_rids = pd.DataFrame(protein_vars.values(), index=list(protein_vars), columns=cols)
        print(f'{len(df_add_rids):4d} protein variables to add')
        self.model.add_reactions(df_add_rids)

    def query_kcats(self):
        """Query enzyme-reaction kcat values.

        Generate a table that can be used to configure kcats

        :return: table with kcat values
        :rtype: pandas DataFrame
        """
        notes = f'kcat export from {self.model.model_attrs.id}'
        rid_kcats = []
        for sub_rid, r in self.model.reactions.items():
            tmp_rid = re.sub('_REV$', '', sub_rid)
            rid = re.sub(r'_iso\d$', '', tmp_rid)
            if len(r.enzymes) == 1:
                dirxn = -1 if re.search('_REV$', sub_rid) else 1
                enzyme = ', '.join([locus.strip() for locus in r.enzymes[0].split('_')[1:]])
                rid_kcats.append([rid, sub_rid, dirxn, enzyme, r.kcatf, notes])
        df_rid_kcats = pd.DataFrame(rid_kcats, columns=['rid', 'sub_rid', 'dirxn', 'enzyme', 'kcat', 'notes'])
        df_rid_kcats.sort_values(by=['sub_rid', 'dirxn'], inplace=True)
        df_rid_kcats.set_index('rid', inplace=True)
        return df_rid_kcats

    def validate(self):
        """Validate compliance to SBML standards

        :return: flag of compliance
        :rtype: bool
        """
        return self.model.validate()

    def export(self, fname):
        """Export ec model to either SBML or Excel

        :param fname: filename with extension .xml or .xlsx
        :type fname: str
        :return: success flag
        :rtype: bool
        """
        return self.model.export(fname)
