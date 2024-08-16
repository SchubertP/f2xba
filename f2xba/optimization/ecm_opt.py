"""Implementation of EcmOptimization class.

Support functions for CobraPy and GurobiPy optimization of EC (enzyme constraint) models
that have been created using the f2xba modelling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import re
import numpy as np
import pandas as pd
from collections import defaultdict
# Gurobipy should not be a hard requirement, unless used in this context
try:
    import gurobipy as gp
except ImportError:
    gp = None
    pass

import f2xba.prefixes as pf
from .optimize import Optimize, status2text


class EcmOptimization(Optimize):

    def __init__(self, fname, cobra_model=None):
        """Instantiate EcmOptimization class

        :param fname: path name of SBML coded metabolic model
        :type fname: str
        :param cobra_model: (optional) the corresponding cobra model
        :type cobra_model: cobra.core.model.Model if supplied
        """
        super().__init__(fname, cobra_model)
        self.orig_coeffs = {}

        self.ecm_type = self.m_dict['modelAttrs'].get('id', '_GECKO').rsplit('_', 1)[1]
        if self.ecm_type.endswith('MOMENT'):
            self.configure_moment_model_constraints()

    def configure_moment_model_constraints(self):
        """configure constraints related to MOMENT modelling
            in MOMENT formalism, promiscuous enszymes can catalyze alternative reactions
            without additional cost. Only most costly reacions flux need to be supported
            and all alternative reactions come for free
        """
        if self.is_gpm:
            for constr in self.gpm.getConstrs():
                if re.match(pf.C_prot_, constr.ConstrName) and re.match(pf.C_prot_pool, constr.ConstrName) is None:
                    constr.sense = '>'
        else:
            for constr in self.model.constraints:
                if re.match(pf.C_prot_, constr.name) and re.match(pf.C_prot_pool, constr.name) is None:
                    constr.ub = 1000.0
        print(f'MOMENT protein constraints configured ≥ 0')

    # Support functions to tune model kcat parameters by rescaling kcat values
    def gp_tune_scale_kcats(self, scale_kcats):
        """Scale selected kcats of the model as per provided scale_kcats.

        Used to manually tune kcat parameters.
        - call gp_tune_reset_kcats() after optimizatin to restore original kcat values
        - example:
          scale_kcats= {'BPNT': 0.2, 'TMPK':0.2, 'TMPK_REV':0.2, 'CLPNS':0.25}
          geo.gp_tune_scale_kcats(scale_kcats)
          solution = geo.optimize()
          geo.gp_tune_reset_kcats()

        Steps:
        - based on reaction id, access model variable and get column entries of S-matrix
        - identify matrix entries related to protein constraints
        - retrieve current coefficient value, so it can be restored after optimization
        - set new coefficient using the scaling value

        :param scale_kcats: selected reaction ids (without 'R_' prefix) with kcat scaling factor
        :type scale_kcats: dict (key: reaction id without 'R_' prefix, val: float)
        """
        assert self.is_gpm is True, "kcat tunining supported on Gurobipy model only"
        orig_coeffs = defaultdict(dict)
        self.gpm.update()
        for ridx, scale in scale_kcats.items():
            var = self.gpm.getVarByName(f'{pf.R_}{ridx}')
            col = self.gpm.getCol(var)
            for idx in range(col.size()):
                constr = col.getConstr(idx)
                if re.match(f'{pf.C_prot_}', constr.getAttr('ConstrName')):
                    coeff = col.getCoeff(idx)
                    orig_coeffs[var][constr] = coeff
                    self.gpm.chgCoeff(constr, var, coeff / scale)
        self.gpm.update()
        self.orig_coeffs = dict(orig_coeffs)

    def gp_tune_reset_kcats(self):
        """reset kcat values to original values.

        Used in manually tuning model kcats
        - call gp_tune_scale_kcats() prior to optmization
        - call gp_tune_reset_kcats() after optimizatin to reset old kcat values
        """
        assert self.is_gpm is True, "kcat tunining supported on Gurobipy model only"
        for var, data in self.orig_coeffs.items():
            for constr, orig_coeff in data.items():
                self.gpm.chgCoeff(constr, var, orig_coeff)
        self.gpm.update()
        self.orig_coeffs = {}

    def gp_fit_kcats(self, proteins_mg_per_gdw, max_scale_factor=1000.0):

        # create mapping of gene to uniprot id based on V_PC_xxxx variables
        gene2uid = {}
        for rid, row in self.m_dict['reactions'].iterrows():
            if re.match(pf.V_PC_, rid) and type(row['fbcGeneProdAssoc']) is str:
                uid = re.sub(f'^{pf.V_PC_}', '', rid)
                gpid = re.sub('^assoc=', '', row['fbcGeneProdAssoc'])
                if gpid in self.m_dict['fbcGeneProducts'].index:
                    gene = self.m_dict['fbcGeneProducts'].at[gpid, 'label']
                    gene2uid[gene] = uid

        # create slack model from original model
        self.gpm.update()
        slack_gpm = self.gpm.copy()

        # add slack variables (positive and negative)
        slack_vars = []
        for gene, mg_per_gdw in proteins_mg_per_gdw.items():
            uid = gene2uid[gene]

            # fix protein concentrations to proteomics values (selected proteins only)
            var = slack_gpm.getVarByName(f'{pf.V_PC_}{uid}')
            var.lb = mg_per_gdw
            var.ub = mg_per_gdw

            # create positive and negative slack variables for protein concentrations
            pslack_var = slack_gpm.addVar(lb=0.0, ub=100, vtype='C', name=f'{pf.V_PSLACK_}{uid}')
            nslack_var = slack_gpm.addVar(lb=0.0, ub=100, vtype='C', name=f'{pf.V_NSLACK_}{uid}')
            slack_vars.append(pslack_var)
            slack_vars.append(nslack_var)

            # couple pos/neg slack variable with protein concentration constraint
            constr = slack_gpm.getConstrByName(f'{pf.C_prot_}{uid}')
            slack_gpm.chgCoeff(constr, pslack_var, 1.0)
            slack_gpm.chgCoeff(constr, nslack_var, -1.0)

        # 1. optimize for maximum growth rate under protein constraint conditions
        slack_gpm.optimize()
        status = status2text[slack_gpm.status]
        if status.lower() != 'optimal':
            print(f'slack model optimization (max growth rate) returned status: {status}')
            return {}

        max_gr = slack_gpm.objval
        print(f'successful slack optimization with growth rate of {max_gr:.3f} h-1')

        # 2. minimize sums of slack variables to avoid loops
        # fix maximum growth rate by adding a new constraint
        model_objective = slack_gpm.getObjective()
        slack_gpm.addLConstr(model_objective, '=', max_gr, 'Max growth rate objective')

        # New Objective: minimize sum of all non-negative reactions
        slack_gpm.update()

        slack_objective = gp.LinExpr(np.ones(len(slack_vars)), slack_vars)
        slack_gpm.setObjective(slack_objective, gp.GRB.MINIMIZE)

        # optimize and collect results (without error handling)
        slack_gpm.optimize()
        status = status2text[slack_gpm.status]
        if status.lower() != 'optimal':
            print(f'slack model optimization (minimize sum of slack) returned status: {status}')
            return {}

        # process slack model results to determine kcat scaling factors
        gene2scale = {}
        for gene, exp_mg_per_gdw in proteins_mg_per_gdw.items():
            uid = gene2uid[gene]
            pos_slack = slack_gpm.getVarByName(f'{pf.V_PSLACK_}{uid}').x
            neg_slack = slack_gpm.getVarByName(f'{pf.V_NSLACK_}{uid}').x
            pred_mg_per_gdw = exp_mg_per_gdw + pos_slack - neg_slack

            scale = pred_mg_per_gdw / exp_mg_per_gdw
            scale = min(scale, max_scale_factor)
            # only capture postive >0 values for scale, i.e. predicted concentration > 0
            if pred_mg_per_gdw > 1e-8:
                gene2scale[gene] = scale

        min_scale = min(gene2scale.values())
        max_scale = max(gene2scale.values())
        print(f'{len(gene2scale)} genes to be scaled by factors of {min_scale:.4f} ... {max_scale:.3f}')
        return gene2scale

    @staticmethod
    def create_fitted_kcats_file(gene2scale, orig_kcats_fname, fitted_kcats_fname):
        """Scale original kcats for selected genes and create fitted kcats Excel document

        - genes products (proteins) can be part of one or several enzymes
        - enzyme may catalyze several reactions, including forward/reverse split isoreactions
        - for each selected gene with a kcat scaling factor, retrieve respective kcat records from original
          kcats Excel Document (one record per isoenzyme)
        - calculated fitted kcats for each kcat record where gene product is used.
        - in case of an enzyme complex consisting of several gene products used in an isoreaction, calculate a single
          kcat value, here we used geometric mean across all individual kcats for the kcat record
        - we also limit kcats to be between 1000 and 0.0001 s-1
        - configure scaled kcat records and store resulting fitted kcat Excel document

        :param gene2scale: gene labels (proteins) with kcat scaling factors
        :type gene2scale: dict (key: gene label, val: kcat scaling factor)
        :param orig_kcats_fname: path/file name of Excel spreadsheet with originao model kcats
        :type orig_kcats_fname: str
        :param fitted_kcats_fname: path/file name of Excel spreadsheet with scaled kcats
        :type fitted_kcats_fname: str
        """

        # load model kcats from which to extract selected transporter enzymes
        with pd.ExcelFile(orig_kcats_fname) as xlsx:
            df_orig_kcats = pd.read_excel(xlsx, sheet_name='kcats', index_col=0)
            print(f'{len(df_orig_kcats):4d} original kcat records loaded from {orig_kcats_fname}')

        gene2records = defaultdict(list)
        for record_id, row in df_orig_kcats.iterrows():
            genes = [item.strip() for item in row['info_genes'].split(',')]
            for gene in genes:
                gene2records[gene].append(record_id)
        gene2records = dict(gene2records)
        # (46) gene2records = {'JCVISYN3A_0641': ['R_5FTHFabc', 'R_COAabc', 'R_NACabc', 'R_P5Pabc', 'R_RIBFLVabc'], ...}

        record2fittedkcats = defaultdict(list)
        for gene, scale in gene2scale.items():
            for record_id in gene2records[gene]:
                orig_kcat = df_orig_kcats.at[record_id, 'kcat_per_s']
                record2fittedkcats[record_id].append(orig_kcat * scale)

        # get a unique fitted kcat value (applicable to enzyme complexes): min, max, mean or geometric mean
        # print(f'{" ":42}  min    max   mean  geo_mean')
        record2fittedkcat = {}
        for record_id, kcats in record2fittedkcats.items():
            if re.match(r'R_\w{3}t2r_', record_id) is None:
                # kcat = min(kcats)
                # kcat = max(kcats)
                # kcat = np.mean(kcats)
                kcat = np.exp(sum([np.log(kcat) for kcat in kcats])/len(kcats))
                kcat_rounded = min(1000.0, max(0.0001, round(kcat, 4)))
                record2fittedkcat[record_id] = kcat_rounded

        # create a new kcats file with updated values
        df_fitted_kcats = df_orig_kcats.copy()
        for record_id, kcat in record2fittedkcat.items():
            df_fitted_kcats.at[record_id, 'kcat_per_s'] = kcat
            df_fitted_kcats.at[record_id, 'notes'] = 'automatically fitted to proteomics'
        print(f'{len(record2fittedkcat):4d} kcats fitted')

        # export fitted kcats file
        with pd.ExcelWriter(fitted_kcats_fname) as writer:
            df_fitted_kcats.to_excel(writer, sheet_name='kcats')
            print(f'fitted kcats exported to {fitted_kcats_fname}')
