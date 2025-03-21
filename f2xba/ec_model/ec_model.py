"""Implementation of EcModel class.

Extend the XbaModel to become an enzyme constraint model.
The implementation of GECKO is based on the GECKO 1.0 matlab package (Sánchez et al., 2017).
The implementation of ccFBA and MOMENTmr is based on the R-package sybilccFBA (Desouki, 2015)
The implementation of MOMENT is based on literature (Adadi et al., 2012)

Peter Schubert, HHU Duesseldorf, October 2023
"""

import re
import numpy as np
import pandas as pd
import f2xba.prefixes as pf
from ..utils.mapping_utils import valid_sbml_sid, load_parameter_file

MAX_CONC_PROT = 100  # mg/gDW maximum protein concentration in kpmf (kilo protein mass fraction)


class EcModel:
    """Extend a XbaModel to an enzyme constraint model.

    The implementation of GECKO is based on the GECKO 1.0 matlab package (Sánchez et al., 2017),
    ccFBA and MOMENTmr on the R-package sybilccFBA (Desouki, 2015) and MOMENT on literature
    (Adadi et al., 2012).

    .. code-block:: python

        xba_model = XbaModel('iJO1366.xml')
        xba_model.configure('iJO1366_GECKO_xba_parameters.xlsx')

        ec_model = EcModel(xba_model)
        ec_model.configure('iJO1366_GECKO_ecm_parameters.xlsx')
        if ec_model.validate():
            ec_model.export('iJO1366_GECKO.xml')
    """

    def __init__(self, xba_model):
        """Instantiate the EcModel instance.

        :param xba_model: a reference to a XbaModel instance
        :type xba_model: :class:`f2xba.XbaModel`
        """
        self.model = xba_model
        self.arm_flag = False

    def configure(self, fname):
        """Configuration with ECM configuration data.

        Accepted table: 'general'

        :param str fname: filename of ECM configuration file (.xlsx)
        :return: success
        :rtype: bool
        """
        sheet_names = ['general', 'rescale_reactions']
        ecm_params = load_parameter_file(fname, sheet_names)
        if 'general' not in ecm_params.keys():
            print(f'mandatory table "general" not found in the document')
            raise ValueError

        general_params = ecm_params['general']['value'].to_dict()
        ecm_type = general_params.get('ecm_type', 'GECKO')
        p_total = general_params.get('p_total', 0.5)
        avg_enz_sat = general_params.get('avg_enz_sat', 0.5)
        self.arm_flag = general_params.get('arm_flag', False)

        if 'pm2totpm' in general_params:
            pm2totpm = general_params['pm2totpm']
        else:
            assert 'pm2totpm_val_or_paxdb' in general_params, 'parameter "pm2totpm_val_or_paxdb" is required'
            if type(general_params['pm2totpm_val_or_paxdb']) is float:
                pm2totpm = general_params['pm2totpm_val_or_paxdb']
            else:
                pm2totpm = self.get_modelled_protein_mass_fraction(general_params['pm2totpm_val_or_paxdb'])
        print(f'modeled protein fraction of total protein mass {pm2totpm:.4f} g/g')

        if 'rescale_reactions' in ecm_params:
            self._rescale_reactions(ecm_params['rescale_reactions'])

        # create enzyme constraint model
        if ecm_type in ['GECKO', 'MOMENTmr']:
            self.model.add_isoenzyme_reactions()
            self.model.make_irreversible()
            self._add_gecko_protein_species()
        elif ecm_type == 'MOMENT':
            self.model.add_isoenzyme_reactions()
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

        if self.arm_flag is True:
            self.model.add_arm_reactions()

        protein_pool = p_total * pm2totpm
        print(f'protein constraint: {protein_pool*1000.0:.2f} mg/gDW - modeled protein')
        self._add_total_protein_constraint(protein_pool)
        self._reaction_enzyme_coupling(avg_enz_sat)

        # remove model components that are not used in the model
        self.model.clean()

        # add some parameter values for reference (after model.clean() so they are not removed)
        add_params = {'gmodeledP_per_gP': {'value': pm2totpm, 'name': 'protein mass fraction modeled'},
                      'gP_per_gDW': {'value': p_total, 'name': 'protein mass fraction of total dry mass'},
                      'avg_enz_sat': {'value': avg_enz_sat, 'name': 'average enzyme saturation level'}}
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
            p = self.model.proteins[uid]
            prot_sid = pf.C_prot_ + valid_sbml_sid(uid)
            p.link_sid(prot_sid)
            prot_metaid = f'meta_prot_{valid_sbml_sid(uid)}'
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
                    p = self.model.proteins[uid]
                    prot_sid = pf.C_prot_ + valid_sbml_sid(uid) + f'_{rid}'
                    p.link_sid(prot_sid)
                    prot_metaid = f'meta_prot_{valid_sbml_sid(uid)}_{rid}'
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
        """add total modeled protein constraint to the enzyme constraint model.

        add protein pool species to the model
        add protein drain reactions to the model
        add protein pool exchange reation to the model

        protein concentration in model has units of mg/gDW (i.e. 1000.0 * pmf)

        :param float total_protein: total modeled protein mass in g/gDW
        """
        # add total protein pool species in default compartment
        pool_sid = f'{pf.C_prot_pool}'
        pool_name = 'total protein pool'
        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame([[pool_name, self.model.main_cid, False, False, False]],
                                      index=[pool_sid], columns=cols)
        self.model.add_species(df_add_species)

        # add total protein pool exchange reaction
        protein_vars = {}
        draw_rid = pf.V_PC_total
        draw_name = f'total concentration modeled protein'
        zero_mg_pid = self.model.get_fbc_bnd_pid(0.0, 'mg_per_gDW', 'zero_mass_conc_mg_per_gDW')
        total_prot_mg_pid = self.model.get_fbc_bnd_pid(total_protein * 1000.0, 'mg_per_gDW', 'conc_protein_mg_per_gDW')
        protein_vars[draw_rid] = [draw_name, f' => {pool_sid}', zero_mg_pid, total_prot_mg_pid, None, 'protein']

        # add protein drain reactions - supporting MOMENT with specific protein split in protein species per reaction
        max_conc_mg_pid = self.model.get_fbc_bnd_pid(MAX_CONC_PROT, 'mg_per_gDW', 'max_conc_prot_mg_per_gDW')

        # create only for used proteins
        used_enz_gps = self._get_enzyme_gene_products()
        for gpid in used_enz_gps:
            uid = self.model.gps[gpid].uid
            p = self.model.proteins[uid]
            draw_rid = pf.V_PC_ + valid_sbml_sid(uid)
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

        :param str fname: file path/name of protein abundance data collected from Pax-db.org
        :return: relative pmf of model based
        :rtype: float
        """
        # parse PacDb file
        name = ''
        publication = ''
        last_comment_line = ''
        gene_idx = None
        ppm_abundance = {}
        with open(fname) as file:
            for line in file:
                line = line.rstrip()
                if re.match('#', line):
                    if re.match('#name: ', line):
                        name = line.split(': ')[1]
                    elif re.match('#publication_year: ', line):
                        publication = line.split(': ')[1]
                    last_comment_line = line
                else:
                    if gene_idx is None:
                        cols = re.sub('^#', '', last_comment_line).split('\t')
                        gene_idx = cols.index('string_external_id')
                        ppm_idx = cols.index('abundance')
                    record = line.split('\t')
                    gene = record[gene_idx].split('.')[1]
                    ppm_abundance[gene] = float(record[ppm_idx])
        print(f'PaxDb file {name}, year {publication}, covering {len(ppm_abundance)} proteins '
              f'(total {sum(ppm_abundance.values()):.1f} ppm) loaded from {fname}')

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

    def _reaction_enzyme_coupling(self, avg_enz_sat):
        """Couple reactions with enzyme/protein requirement via kcats.

        applicable to enzyme catalyzed reactions.
        from reaction get enzyme and corresponding kcat
        convert kcat from s-1 to h-1
        scaled total enzyme kcat by number of active sites and average enzyme saturation

        protein constraint for a specific protein:
          reaction flux [mmol/gDWh] ≤ kcat [1/s] * 3600 [s/h] * n_AS * 1/stoic * avg_enz_sat * P_conc [mg/gDW]
                                      * 1/MW [mmol/mg]
          - stoic: number of protein copies in the enzyme
          - n_AS: number of active sites
          - MW: protein molecular weight in Da = g/mol = mg/mmol
        Except for MOMENT we assume equality to optimize protein allocation
          constraint per protein: C_prot_uid couples reaction flux R_xxx with protein concentration V_PC_uid
          C_prot_uid: (-1/(kcat * 3600 * n_AS * avg_enz_sat) * stoic * MW) * R_xxx + 1 * V_PC_uid = (≤) 0

        :param float avg_enz_sat: average enzyme saturation applicable to all reactions
        """
        for r in self.model.reactions.values():
            assert len(r.enzymes) <= 1, ('Something went wrong, at this stage of isoreactions, we expect a maximum'
                                         'of one enzyme catalyzing this reaction')
            if len(r.enzymes) == 1 and r.kcat is not None:
                e = self.model.enzymes[r.enzymes[0]]
                enz_kcat_per_h = r.kcat * 3600.0 * e.active_sites * avg_enz_sat

                for locus, stoic in e.composition.items():
                    uid = self.model.locus2uid[locus]
                    p = self.model.proteins[uid]
                    prot_sid = None
                    # in MOMENT model, a promiscuous enzyme/protein has specific proteins per reactions
                    # e.g. b2215, ompC, P06996: this protein is linked to all (reversible) reactions it catalyzes
                    #  reactions R_23CCMPtex_iso4 and R_23CCMPtex_iso4_REV both coupled to constraint
                    #  C_prot_P06996_R_23CAMPtex_iso4
                    if len(p.linked_sids) == 1:
                        prot_sid = list(p.linked_sids)[0]
                    else:
                        rev_rid = re.sub('_REV$', '', r.id)
                        for linked_sid in p.linked_sids:
                            if rev_rid in linked_sid:
                                prot_sid = linked_sid
                                break
                    r.reactants[prot_sid] = (stoic / enz_kcat_per_h) * p.mw

    def _rescale_reactions(self, df_rescale):
        """Rescale reactions (e.g. Biomass) as per ecYeast7 (Sanchez, 2017)

        biomass dict contains
        - biomass reaction id 'rid'
        - rescale tasks 'rescale' a list of dict, each with
            - either rescale factor 'factor' or an absolute value 'value', flaat
            - 'rectants' and/or 'products', consisting of list of species ids

        amino acids get scaled by parameter f_p
        carbohydrates get scaled by parameter f_c
        GAM related species get set to gam level

        :param df_rescale: biomass rescaling parameters
        :type df_rescale: pandas.DataFrame
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
        :rtype: pandas.DataFrame
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
        """Validate compliance to SBML standards.

        :return: success
        :rtype: bool
        """
        return self.model.validate()

    def export(self, fname):
        """Export EcModel to SBML coded file or in tabular format.

        :param str fname: filename (with extension '.xml' or '.xlsx')
        :return: success
        :rtype: bool
        """
        return self.model.export(fname)
