"""Implementation of EcModel class.


Peter Schubert, HHU Duesseldorf, October 2023
"""

import re
import numpy as np
import pandas as pd
import f2xba.prefixes as pf

MAX_CONC_PROT = 100  # mg/gDW maximum protein concentration in kpmf (kilo protein mass fraction)


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

    def configure(self, fname):
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

        :param fname: file name for ECM configuration parameters
        :type fname: str
        :return: success/failure of model configuration
        :rtype: bool
        """
        sheets = ['general', 'rescale_reactions']
        ecm_params = {}
        with pd.ExcelFile(fname) as xlsx:
            for sheet in sheets:
                if sheet in xlsx.sheet_names:
                    ecm_params[sheet] = pd.read_excel(xlsx, sheet_name=sheet, index_col=0)
            print(f'{len(ecm_params)} tables with EC model configuration parameters loaded from {fname}')

        general_params = ecm_params['general']['value'].to_dict()
        ecm_type = general_params.get('ecm_type', 'GECKO')
        p_total = general_params.get('p_total', 0.5)
        avg_enz_sat = general_params.get('avg_enz_sat', 0.5)

        if 'rescale_reactions' in ecm_params:
            self._rescale_reactions(ecm_params['rescale_reactions'])

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

        if 'p_abundance_fname' in general_params:
            pax_db_fname = general_params['p_abundance_fname']
            pm2totpm = self.get_modelled_protein_mass_fraction(pax_db_fname)
            print(f'modeled protein fraction {pm2totpm:.4f} g/g, based on {pax_db_fname}')
        else:
            assert('pm2totpm' in general_params)
            pm2totpm = general_params['pm2totpm']

        protein_pool = p_total * pm2totpm * avg_enz_sat
        print(f'protein constraint: {protein_pool*1000.0:.2f} mg/gDW, '
              f'assuming avg enzyme saturation of {avg_enz_sat:.2f}')

        self._add_total_protein_constraint(protein_pool)
        self._reaction_enzyme_coupling()

        # remove model components that are not used in the model
        self.model.clean()

        # add some parameter values for reference (after model.clean() so they are not removed)
        add_params = {'frac_prot_totprot': {'value': pm2totpm, 'name': 'protein mass fraction modelled'},
                      'frac_totprot_cdw': {'value': p_total, 'name': 'protein mass fraction of total dry mass'},
                      'frac_enzyme_sat': {'value': avg_enz_sat, 'name': 'average enzyme saturation level'}}
        for pid, p_data in add_params.items():
            self.model.add_parameter(pid, p_data)

        # modify some model attributs and create L3V2 SBML model
        self.model.model_attrs['id'] += f'_{ecm_type}'
        if 'name' in self.model.model_attrs:
            self.model.model_attrs['name'] = f'{ecm_type} model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        self.model.print_size()
        return True

    def _get_enzyme_gene_products(self):
        """Determine used genes products for enzymes.

        Note: Model can contain additional gene products for
        gene product associations where enzyme kcat have not
        been provided. Useful for gene deletion studies.

        :return: gene product ids used by the enzymes in the model
        :rtype: set of str
        """
        used_enzymes = set()
        for rid, r in self.model.reactions.items():
            for eid in r.enzymes:
                used_enzymes.add(eid)
        used_enz_gps = set()
        for eid in used_enzymes:
            e = self.model.enzymes[eid]
            for gene_id in e.composition:
                used_enz_gps.add(self.model.locus2gp[gene_id])
        return used_enz_gps

    def _add_gecko_protein_species(self):
        """Add protein species to the GECKO/ccFBA/MOMENTmr model.

        Note: Protein compartments in XBA model are derived from
        reaction compartments, which in turn is a concatenation of
        species compartments being substrates to the reaction.

        Here we only need one of the species compartments.

        We only create proteins when gene product is used by any of the
        enzymes

        add protein species to the model (also ensure that protein species name is valid)
        add protein drain reactions to the model
        """
        used_enz_gps = self._get_enzyme_gene_products()

        protein_sids = {}
        for gpid in used_enz_gps:
            uid = self.model.gps[gpid].uid
            checked_uid = re.sub(r'\W', '', uid)
            p = self.model.proteins[uid]
            prot_sid = f'{pf.M_prot}_{checked_uid}_{p.cid}'
            p.link_sid(prot_sid)
            prot_metaid = f'meta_prot_{checked_uid}'
            prot_annot = f'bqbiol:is, uniprot/{uid}'
            protein_sids[prot_sid] = [p.name, p.cid, False, False, False, prot_metaid, prot_annot]

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition',
                'constant', 'metaid', 'miriamAnnotation']
        df_add_species = pd.DataFrame(protein_sids.values(), index=list(protein_sids), columns=cols)
        print(f'{len(df_add_species):4d} protein constraints to add')
        self.model.add_species(df_add_species)

    def _add_moment_protein_species(self):
        """Add protein species to MOMENT model.

        Execute before making model irreversilbe.

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
                    checked_uid = re.sub(r'\W', '', uid)
                    p = self.model.proteins[uid]
                    prot_sid = f'{pf.M_prot}_{checked_uid}_{rid}_{p.cid}'
                    p.link_sid(prot_sid)
                    prot_metaid = f'meta_prot_{checked_uid}_{rid}'
                    prot_annot = f'bqbiol:is, uniprot/{uid}'
                    protein_sids[prot_sid] = [p.name, p.cid, False, False, False, prot_metaid, prot_annot]

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition',
                'constant', 'metaid', 'miriamAnnotation']
        df_add_species = pd.DataFrame(protein_sids.values(), index=list(protein_sids), columns=cols)
        print(f'{len(df_add_species):4d} moment protein constraints to add')
        self.model.add_species(df_add_species)

    def _select_ccfba_isoenzyme(self):
        """Select least cost enzyme for reactions catalyzed by several isoenzymes.

        Consider number of active sites when determining enzyme kcat value
        Other isoenzymes get deleted and gpa updated.

        Integration with Thermodynamic FBA (TFA):
            TFA splits TD related reations in fwd/rev (_REV)
            TFA also split reactions that in base model are not reversible
            i.e. we may have irreversible '<rid>_REV' with no kcat

        Selection is based on molecular weight and turnover number.
        Select enzyme with minimal MW/kcat.

        Different isoenzymes might be selected for forward/reverse direction.
        Therefore, reversible reactions with isoenzymes get split in fwd/bwd reaction.
        """
        # Note model.reactions get updated within the loop, therefore fix rids initially

        rids = list(self.model.reactions)
        for rid in rids:
            r = self.model.reactions[rid]
            if len(r.enzymes) >= 2 and r.kcatsf is not None and np.isfinite(r.kcatsf[0]):
                directional_rs = [r]
                if r.reversible:
                    # split reaction, i.e. create a new irreversible reverse reaction
                    rev_r = self.model.split_reversible_reaction(r)
                    directional_rs.append(rev_r)

                for dir_r in directional_rs:
                    valid_kcatsf = (sum([1 for kcat in dir_r.kcatsf if np.isfinite(kcat)])
                                    if dir_r.kcatsf is not None else 0)
                    # consider reactions with undefined kcatf (due to TD integration)
                    if len(r.enzymes) == valid_kcatsf:
                        idx = np.argmin([self.model.enzymes[eid].mw / (kcat * self.model.enzymes[eid].active_sites)
                                         for eid, kcat in zip(dir_r.enzymes, dir_r.kcatsf)])
                        eid = dir_r.enzymes[idx]
                        gpa = ' and '.join(sorted([self.model.uid2gp[self.model.locus2uid[locus]]
                                                   for locus in self.model.enzymes[eid].composition]))
                        if ' and ' in gpa:
                            gpa = '(' + gpa + ')'
                        dir_r.enzymes = [dir_r.enzymes[idx]]
                        dir_r.kcatsf = [dir_r.kcatsf[idx]]
                        dir_r.gene_product_assoc = gpa
                    else:
                        dir_r.enzymes = []
                        dir_r.kcatsf = None
                        dir_r.kcatsr = None

    def _add_total_protein_constraint(self, total_protein):
        """add total protein constraint to the enzyme constraint model.

        add protein pool species to the model
        add protein drain reactions to the model
        add protein pool exchange reation to the model

        protein concentration in model has units of mg/gDW (i.e. 1000.0 * pmf)

        :param total_protein: total active protein mass balance in g/gDW
        :type total_protein: float
        """
        # add total protein pool species in default compartment
        pool_cid = sorted([(count, cid) for cid, count in self.model.get_used_cids().items()])[-1][1]
        pool_sid = f'{pf.M_prot_pool}_{pool_cid}'
        pool_name = 'total protein pool'
        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame([[pool_name, pool_cid, False, False, False]],
                                      index=[pool_sid], columns=cols)
        self.model.add_species(df_add_species)

        # add total protein pool exchange reaction
        protein_vars = {}
        draw_rid = pf.V_PC_total_active
        draw_name = f'total concentration active protein'
        zero_mg_pid = self.model.get_fbc_bnd_pid(0.0, 'mg_per_gDW', 'zero_mass_conc_mg')
        total_prot_mg_pid = self.model.get_fbc_bnd_pid(total_protein * 1000.0, 'mg_per_gDW', 'conc_active_protein_mg')
        protein_vars[draw_rid] = [draw_name, f' => {pool_sid}', zero_mg_pid, total_prot_mg_pid, None, 'protein']

        # add protein drain reactions - supporting MOMENT with specific protein split in protein species per reaction
        max_conc_mg_pid = self.model.get_fbc_bnd_pid(MAX_CONC_PROT, 'mg_per_gDW', 'max_conc_prot_mg')

        # create only for used proteins
        used_enz_gps = self._get_enzyme_gene_products()
        for gpid in used_enz_gps:
            uid = self.model.gps[gpid].uid
            p = self.model.proteins[uid]
            checked_uid = re.sub(r'\W', '', uid)
            draw_rid = f'{pf.V_PC}_{checked_uid}'
            draw_name = f'conc_prot_{uid}'
            # supporting MOMENT model with protein split into protein species per reaction
            products = ' + '.join([sid for sid in p.linked_sids])
            # reaction_string = f'{p.mw / 1000.0} {pool_sid} => {products}'
            reaction_string = f'{pool_sid} => {products}'
            protein_vars[draw_rid] = [draw_name, reaction_string, zero_mg_pid, max_conc_mg_pid,
                                      self.model.uid2gp[uid], 'protein']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'fbcGeneProdAssoc', 'kind']
        df_add_rids = pd.DataFrame(protein_vars.values(), index=list(protein_vars), columns=cols)
        print(f'{len(df_add_rids):4d} protein variables to add')
        self.model.add_reactions(df_add_rids)

    def get_modelled_protein_mass_fraction(self, fname):
        """Determine protein mass fraction based on Pax-db.org downlaod.

        Determine which part of organism protein mass fraction is modelled.
        Protein abundance file needs first to be downloaded from Pax-db.org

        :param fname: file path/name of protein abundance data collected from Pax-db.org
        :type fname: str
        :return: relative pmf of model based
        :rtype: float
        """
        # parse file
        df_ppm = pd.read_table(fname, comment='#', header=None, usecols=[1, 2], index_col=0)
        df_ppm.index = [locus.split('.')[1] for locus in df_ppm.index]
        ppm_abundance = df_ppm.iloc[:, 0].to_dict()

        used_enz_gps = self._get_enzyme_gene_products()
        used_uids = {self.model.gps[gpid].uid for gpid in used_enz_gps}

        p_rel_total = 0.0
        p_rel_model = 0.0
        for locus, ppm in ppm_abundance.items():
            if locus in self.model.uniprot_data.locus2uid:
                uid = self.model.uniprot_data.locus2uid[locus]
                rel_mass = ppm / 1.0e6 * self.model.uniprot_data.proteins[uid].mass
                p_rel_total += rel_mass
                if uid in used_uids:
                    p_rel_model += rel_mass
        return p_rel_model / p_rel_total

    def _reaction_enzyme_coupling(self):
        """Couple reactions with enzyme/protein requirement via kcats.

        applicable to enzyme catalyzed reactions.
        from reaction get enzyme and corresponding kcat
        convert kcat from s-1 to h-1
        scaled total enzyme kcat by number of active sites

        protein constraints are in units of mg/gDW
        1/kcat [s] * 1/3600 [h/s]  * 1/n_AS * enz_stoic * flux [mmol/gDWh] * MW [mg/mmol] = 1 * V_PC_xxxx [mg/gDW]

        from enzyme get proteins and their stoichiometry in the enzyme complex
        add to reaction reactants the proteins with 1/kcat_per_h scaled by stoic.
        """
        # scale = 1e3
        for r in self.model.reactions.values():
            assert (len(r.enzymes) <= 1)
            if len(r.enzymes) == 1 and r.kcat is not None:
                e = self.model.enzymes[r.enzymes[0]]
                enz_kcat_per_h = r.kcat * 3600.0 * e.active_sites

                for locus, stoic in e.composition.items():
                    uid = self.model.locus2uid[locus]
                    p = self.model.proteins[uid]
                    prot_sid = None
                    # in MOMENT model, a promiscuous enzyme/protein has specific proteins per reactions
                    # e.g. b2215, ompC, P06996: this protein is linked to all (reversible) reactions it catalyzes
                    #  reactions R_23CCMPtex_iso4 and R_23CCMPtex_iso4_REV both coupled to constraint
                    #  M_prot_P06996_R_23CAMPtex_iso4_e
                    if len(p.linked_sids) == 1:
                        prot_sid = list(p.linked_sids)[0]
                    else:
                        rev_rid = re.sub('_REV$', '', r.id)
                        for linked_sid in p.linked_sids:
                            if rev_rid in linked_sid:
                                prot_sid = linked_sid
                                break
                    # r.reactants[prot_sid] = stoic / enz_kcat_per_h * scale
                    r.reactants[prot_sid] = stoic / enz_kcat_per_h * p.mw

    def _rescale_reactions(self, df_rescale):
        """rescale reactions (e.g. Biomass) as per ecYeast7 (Sanchez, 2017)

        biomass dict contains
        - biomass reaction id 'rid'
        - rescale tasks 'rescale' a list of dict, each with
            - either rescale factor 'factor' or an absolute value 'value', flaat
            - 'rectants' and/or 'products', consisting of list of species ids

        amino acids get scaled by parameter f_p
        carbohydrates get scaled by parameter f_c
        GAM related species get set to gam level

        :param df_rescale: biomass rescaling parameters
        :type df_rescale: pandas DataFrame
        :return:
        """
        modify_reactants = {}
        modify_products = {}
        for rid, row in df_rescale.iterrows():
            r = self.model.reactions[rid]
            if type(row['reactants']) is str:
                scale_reacs = {sid.strip() for sid in row['reactants'].split(',')}
                if rid not in modify_reactants:
                    modify_reactants[rid] = r.reactants.copy()
                new_srefs = {}
                if np.isfinite(row['factor']):
                    factor = row['factor']
                    for sid, stoic in modify_reactants[rid].items():
                        new_srefs[sid] = factor * stoic if sid in scale_reacs else stoic
                elif np.isfinite(row['value']):
                    value = row['value']
                    for sid, stoic in modify_reactants[rid].items():
                        new_srefs[sid] = value if sid in scale_reacs else stoic
                modify_reactants[rid] = new_srefs

            if type(row['products']) is str:
                scale_prods = {sid.strip() for sid in row['products'].split(',')}
                if rid not in modify_products:
                    modify_products[rid] = r.products.copy()
                new_srefs = {}
                if np.isfinite(row['factor']):
                    factor = row['factor']
                    for sid, stoic in modify_products[rid].items():
                        new_srefs[sid] = factor * stoic if sid in scale_prods else stoic
                elif np.isfinite(row['value']):
                    value = row['value']
                    for sid, stoic in modify_products[rid].items():
                        new_srefs[sid] = value if sid in scale_prods else stoic
                modify_products[rid] = new_srefs

        modify_attrs = []
        for rid, srefs in modify_reactants.items():
            reactants = '; '.join([f'species={sid}, stoic={stoic}' for sid, stoic in srefs.items()])
            modify_attrs.append([rid, 'reaction', 'reactants', reactants])
        for rid, srefs in modify_products.items():
            products = '; '.join([f'species={sid}, stoic={stoic}' for sid, stoic in srefs.items()])
            modify_attrs.append([rid, 'reaction', 'products', products])
        cols = ['id', 'component', 'attribute', 'value']
        df_modify_attrs = pd.DataFrame(modify_attrs, columns=cols)
        df_modify_attrs.set_index('id', inplace=True)
        print(f'reaction rescaling')
        self.model.modify_attributes(df_modify_attrs, 'reaction')

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
                eid = r.enzymes[0]
                genes = ', '.join([locus.strip() for locus in eid.split('_')[1:]])
                active_sites = self.model.enzymes[eid].active_sites,
                rid_kcats.append([rid, sub_rid, dirxn, genes, r.kcatsf, active_sites, notes])
        cols = ['rid', 'sub_rid', 'dirxn', 'genes', 'kcat_per_s', 'active_sites', 'notes']
        df_rid_kcats = pd.DataFrame(rid_kcats, columns=cols)
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
