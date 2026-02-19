Quick User Guide
================


Extended model creation
-----------------------

This concise user guide offers Python code for the development of different types of extended metabolic models. However, it should be noted that configuration data must be created in advance. The reader is advised to proceed through the tutorial sections, where the process of generating the necessary configuration data will be elucidated. This will be exemplified by the genome-scale metabolic model of *Escherichia coli*, designated "iML1515". Once an initial set of configuration data has been obtained, it can be readily modified using a spreadsheet editor or programmatically. The subsequent code facilitates the creation of different parametrized model versions.

The creation of an extended genome-scale model adheres to a straightforward workflow. A ``XbaModel`` instance serves as an in-memory representation of the model during the extension process. The genome-scale metabolic model is initially loaded from a SBML file, and the in-memory model is prepared for the extension using the ``xba_model.configure()`` function. Subsequently, the in-memory model is provided as a reference to extended model type-specific classes, such as ``EcModel``, which extend the in-memory model based on supplied configuration data.


Create GECKO model
^^^^^^^^^^^^^^^^^^

``GECKO`` is one of four enzyme constraint model types that have been implemented. The creation of the ``ccFBA``, ``MOMENT``, and ``MOMENTmr`` models is achievable by setting the parameter ``ecm_type`` in the ECM configuration file accordingly.

.. code-block:: python

   from f2xba import XbaModel, EcModel

   xba_model = XbaModel('iML1515.xml')
   xba_model.configure('iML1515_GECKO_xba_parameters.xlsx')

   ec_model = EcModel(xba_model)
   ec_model.configure('iML1515_GECKO_ecm_parameters.xlsx')
   ec_model.export('iML1515_GECKO.xml')


Create RBA model
^^^^^^^^^^^^^^^^

RBA models are generated in a similar manner; however, a ``RbaModel`` instance is required.

.. code-block:: python

   from f2xba import XbaModel, RbaModel

   xba_model = XbaModel('iML1515.xml')
   xba_model.configure('iML1515_RBA_xba_parameters.xlsx')

   rba_model = RbaModel(xba_model)
   rba_model.configure('iML1515_RBA_rba_parameters.xlsx')
   rba_model.export('iML1515_RBA.xml')


Create TFA model
^^^^^^^^^^^^^^^^

TFA models are thermodynamics constraint FBA models. They employ an instance of ``TfaModel`` and may not require specific XBA configuration data.

.. code-block:: python

   from f2xba import XbaModel, TfaModel

   xba_model = XbaModel('iML1515.xml')
   xba_model.configure('iML1515_TFA_xba_parameters.xlsx')

   tfa_model = TfaModel(xba_model)
   tfa_model.configure('iML1515_TFA_tfa_parameters.xlsx')
   tfa_model.export('iML1515_TFA.xml')
  

Create TGECKO model
^^^^^^^^^^^^^^^^^^^

The development of a thermodynamics constraint GECKO model, designated TGECKO, comprises two distinct stages. The initial stage entails the implementation of thermodynamics constraints, followed by the subsequent establishment of GECKO-specific constraints. The configuration data created for the TFA and GECKO models is reused. A similar methodology is employed in the creation of other thermodynamics and enzyme constraint models, such as TMOMTENT. It is noteworthy that the XbaModel is configured with GECKO-specific XBA configuration data.

.. code-block:: python

   from f2xba import XbaModel, TfaModel, EcModel

   xba_model = XbaModel('iML1515.xml')
   xba_model.configure('iML1515_GECKO_xba_parameters.xlsx')

   tfa_model = TfaModel(xba_model)
   tfa_model.configure('iML1515_TFA_tfa_parameters.xlsx')

   ec_model = EcModel(xba_model)
   ec_model.configure('iML1515_GECKO_ecm_parameters.xlsx')
   ec_model.export('iML1515_TGECKO.xml')


Create TRBA model
^^^^^^^^^^^^^^^^^

The generation of a thermodynamics constraint RBA model, designated TRBA, follows the same approach as the construction of a TGECKO model.

.. code-block:: python

   from f2xba import XbaModel, TfaModel, RbaModel

   xba_model = XbaModel('iML1515.xml')
   xba_model.configure('iML1515_RBA_xba_parameters.xlsx')

   tfa_model = TfaModel(xba_model)
   tfa_model.configure('iML1515_TFA_tfa_parameters.xlsx')

   rba_model = RbaModel(xba_model)
   rba_model.configure('iML1515_RBA_rba_parameters.xlsx')
   rba_model.export('iML1515_TRBA.xml')


Extended model optimization
---------------------------

This section presents Python code snippets for loading and optimizing extended metabolic models, which have been generated using the f2xba modeling framework. Two interfaces for model optimization are demonstrated: the well-known `cobrapy <https://cobrapy.readthedocs.io/en/latest/>`_ interface and the `gurobipy <https://www.gurobi.com>`_ interface, which directly accesses the Gurobi numerical solvers. It is imperative that users have already installed cobrapy and/or gurobipy, along with numerical solvers that support linear programming (LP) and mixed integer linear programming (MILP). It is assumed that models are available in SBML encoded files. The objective of optimization is to maximize the growth rate, based on the default media conditions configured in the model.

While the processes of model loading and optimization are specific to the interface used, the extraction of optimization solutions and the presentation of results are independent of the interface. The capacity to process results for multiple optimization solutions in parallel is demonstrated in the tutorials section. Depending on the extended model type, different results can be accessed. It is assumed that the proteomics data is available in ``df_mpmf`` to produce correlation reports and plots for predicted to measured protein concentrations.

Depending on the extended model type that is selected, the ``<model_type>Optimization`` and ``<model_type>Results`` classes will be used to support optimization and results processing, respectively.


Optimize GECKO model using cobrapy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The model is loaded initially through the cobrapy interface. The model is also loaded using by ``EcmOptimization``. While this step is not strictly necessary for the purposes of optimization, the EcmOptimization class retrieves data from the SBML model, which can subsequently be incorporated into the optimization results. The subsequent model optimization will be conducted using the ``ecm.optimize()`` function.

.. code-block:: python

   import cobra
   from f2xba import EcmOptimization

   # load GECKO model
   ecm = cobra.io.read_sbml_model('iML1515_GECKO.xml')
   eo = EcmOptimization('iML1515_GECKO.xml', ecm)

   # optimize GECKO model
   solution = ecm.optimize()
   print(f'predicted growth rate: {solution.objective_value:.3f} h-1 ({solution.status})')


Optimize GECKO model using gurobipy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The gurobipy interface only requires the SBML model to be loaded by ``EcmOptimization``. The model is then optimized by ``eo.optimize()``, a method that returns a solution compatible with cobrapy.

.. code-block:: python

   from f2xba import EcmOptimization

   # load GECKO model using gurobipy
   eo = EcmOptimization('iML1515_GECKO.xml')

   # optimize GECKO model
   solution = eo.optimize()
   print(f'predicted growth rate: {solution.objective_value:.3f} h-1 ({solution.status})')


GECKO optimization results analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Regardless of the optimization interface employed, the optimization results are processed by ``EcmResults`` and a range of results can be retrieved. It is required to provide a reference to the EcmOptimization class and a dictionary of optimization results. The provision of measured protein concentrations, denoted by ``df_mpmf``, is optional.

The retrieval of tables can be performed, and these tables can contain variable values for reaction fluxes, net fluxes, and proteins. When measured protein concentrations are supplied, reports on proteomics correlation and protein levels can be printed, and protein correlation can be plotted. It is noteworthy that plots can be exported to PDF files for future reference. Furthermore, results can be exported to Escher-compliant formats, facilitating subsequent visual representation on metabolic maps.

.. code-block:: python

   from f2xba import EcmResults

   er = EcmResults(eo, {'Glucose': solution}, df_mpmf)

   df_fluxes = er.collect_fluxes()
   df_net_fluxes = er.collect_fluxes(net=True)
   df_proteins = er.collect_protein_results()
   er.report_proteomics_correlation(scale='lin')
   er.report_protein_levels('Glucose')
   er.plot_proteins('Glucose')  
   er.save_to_escher(df_net_fluxes['Glucose'], 'iML1515_GECKO')


Optimize RBA model
^^^^^^^^^^^^^^^^^^

RBA optimization employing the cobrapy interface necessitates support from ``RbaOptimization``, which implements the RBAs bisection algorithm and reconfigures the parameters of the optimization problem during the optimization run. RBA employs external metabolite concentrations as opposed to constraints on the uptake fluxes; consequently, medium conditions are established through the utilization of ``ro.set_medium_conc()``. The bisection optimization is activated using ro.solve(). The utilization of the gurobipy interface is subsequently shown for the optimization of the TRBA model.

The ``RbaResults`` module furnishes access to supplementary optimization outcomes, including enzyme concentrations, RNA concentrations, and compartment occupancy levels.

.. code-block:: python

   import re
   import cobra
   from f2xba import RbaOptimization, RbaResults

   # load RBA model
   rbam = cobra.io.read_sbml_model('iML1515_RBA.xml')
   ro = RbaOptimization('iML1515_RBA.xml', rbam)

   # convert uptake fluxes to relative external medium concentrations
   sigma = ro.avg_enz_saturation
   importer_km = ro.importer_km
   rel_mmol_per_l =  (sigma / (1.0-min(sigma, .99))) * importer_km
   ex_sidx2mmol_per_l = {ex_ridx[3:]: rel_mmol_per_l for ex_ridx in rbam.medium}
   ex_sidx2mmol_per_l['h_e'] = 100.0 * sigma   # H-symport reactions not to be constraint by proton concentration
   ro.set_medium_conc(ex_sidx2mmol_per_l)

   # optimize RBA model
   solution = ro.solve(gr_min=0.01, gr_max=1.2, bisection_tol=1e-3)
   print(f'predicted growth rate: {solution.objective_value:.3f} h-1 ({solution.status})')    

   # RBA optimization results analysis
   rr = RbaResults(ro, {'Glucose': solution}, df_mpmf)
   df_fluxes = rr.collect_fluxes()
   df_all_net_fluxes = rr.collect_fluxes(net=True)    
   metabolic_rids = [rid for rid in list(df_all_net_fluxes.index) if re.match('(PROD)|(DEGR)', rid) is None]
   synthesis_rids = [rid for rid in list(df_all_net_fluxes.index) if re.match('(PROD)|(DEGR)', rid)]
   df_net_fluxes = df_all_net_fluxes.loc[metabolic_rids]
   df_synt_fluxes = df_all_net_fluxes.loc[synthesis_rids]

   # additional results with RBA
   df_proteins = rr.collect_protein_results('mg_per_gP')
   df_enzyme_conc = rr.collect_enzyme_results('µmol_per_gDW')
   df_rna_conc = rr.collect_rna_results('µmol_per_gDW')
   df_occupancy = rr.collect_density_results()
   df_capacity = rr.collect_density_results(capacity=True)

   rr.report_proteomics_correlation(scale='lin')
   rr.report_proteomics_correlation(scale='log')
   rr.report_protein_levels('Glucose')
   rr.plot_proteins('Glucose')

Optimize TFA model
^^^^^^^^^^^^^^^^^^

TFA models can be optimized similarly to FBA models. Optionally, concentration ranges of selected metabolites could be further constraint prior to the optimization. Reactant concentrations can be retrieved from the optimization results.

.. code-block:: python

   import numpy as np
   import cobra
   from f2xba import FbaOptimization, FbaResults

   # Load TFA model
   gem = cobra.io.read_sbml_model('iML1515_TFA.xml')

   # Reconfigure variables and constraints
   fo = FbaOptimization('iML1515_TFA.xml', gem)

   # optional - fix metabolite concentrations to metabolomics
   metabolomics = {'atp_c': 0.00356, 'adp_c': 0.000116, 'gtp_c': 0.001660, 'gdp_c': 0.000203}
   molar_bounds = {sid: (mmol_per_l/1.5, mmol_per_l*1.5) for sid, mmol_per_l in  metabolomics.items()}
   for sid, (lb, ub) in molar_bounds.items():
       fo.reactions.get_by_id('V_LC_' + sid).bounds = (np.log10(lb), np.log10(ub))

   # optimize TFA model
   solution = fo.optimize()
   print(f'predicted growth rate: {solution.objective_value:.3f} h-1 ({solution.status})')

   # TFA optimization results analysis
   fr = FbaResults(fo, {'Glucose': solution})
   df_fluxes = fr.collect_fluxes()
   df_net_fluxes = fr.collect_fluxes(net=True)
   df_species_conc = fr.collect_species_conc()

Optimize TGECKO model
^^^^^^^^^^^^^^^^^^^^^

Optimization of a TD constraint enzyme constraint model is similar to optimizing enzyme constraint models. Optionally, concentration ranges of selected metabolites could be further constraint prior to the optimization. Reactant concentrations can be retrieved from the optimization results.

.. code-block:: python

   import numpy as np
   import cobra
   from f2xba import EcmOptimization, EcmResults

   # Load TGECKO model
   ecm = cobra.io.read_sbml_model('iML1515_TGECKO.xml')

   # Reconfigure variables and constraints
   eo = EcmOptimization('iML1515_TGECKO.xml', ecm)

   # optimize TGECKO model
   solution = eo.optimize()
   print(f'predicted growth rate: {solution.objective_value:.3f} h-1 ({solution.status})')

   # TGECKO optimization results analysis
   er = EcmResults(eo, {'Glucose': solution}, df_mpmf)

   df_fluxes = er.collect_fluxes()
   df_net_fluxes = er.collect_fluxes(net=True)
   df_proteins = er.collect_protein_results()
   df_species_conc = er.collect_species_conc()
   er.report_proteomics_correlation(scale='lin')
   er.report_protein_levels('Glucose')


Optimize TRBA model
^^^^^^^^^^^^^^^^^^^

The optimization of the TRBA is demonstrated using the gurobipy interface, which is significantly faster than using cobrapy. Metabolite concentrations are integrated with the function ``ro.set_tfa_metab_concentrations()``. The optimization results obtained from this method represent a combination of the RBA and TFA optimization results.

.. code-block:: python

   import re
   from f2xba import RbaOptimization, RbaResults

   # load TRBA model - using gurobipy interface
   ro = RbaOptimization('iML1515_TRBA.xml')

   # convert uptake fluxes to relative external medium concentrations
   sigma = ro.avg_enz_saturation
   importer_km = ro.importer_km
   uptake_rids = ro.m_dict['reactions'].loc[ro.uptake_rids]['fbcLb']
   rel_mmol_per_l =  (sigma / (1.0-min(sigma, .99))) * importer_km
   ex_sidx2mmol_per_l = {ex_ridx[5:]: rel_mmol_per_l for ex_ridx in uptake_rids[uptake_rids < 0.0].index}
   ex_sidx2mmol_per_l['h_e'] = 100.0 * sigma   # H-symport reactions not to be constraint by proton concentration
   ro.set_medium_conc(ex_sidx2mmol_per_l)

   # optimize TRBA model
   solution = ro.solve(gr_min=0.01, gr_max=1.2, bisection_tol=1e-3)
   print(f'predicted growth rate: {solution.objective_value:.3f} h-1 ({solution.status})')    

   # TRBA optimization results analysis
   rr = RbaResults(ro, {'Glucose': solution}, df_mpmf)
   df_fluxes = rr.collect_fluxes()
   df_all_net_fluxes = rr.collect_fluxes(net=True)    
   metabolic_rids = [rid for rid in list(df_all_net_fluxes.index) if re.match('(PROD)|(DEGR)', rid) is None]
   synthesis_rids = [rid for rid in list(df_all_net_fluxes.index) if re.match('(PROD)|(DEGR)', rid)]
   df_net_fluxes = df_all_net_fluxes.loc[metabolic_rids]
   df_synt_fluxes = df_all_net_fluxes.loc[synthesis_rids]
   df_proteins = rr.collect_protein_results('mg_per_gP')
   df_enzyme_conc = rr.collect_enzyme_results('µmol_per_gDW')
   df_rna_conc = rr.collect_rna_results('µmol_per_gDW')
   df_occupancy = rr.collect_density_results()
   df_capacity = rr.collect_density_results(capacity=True)
   df_species_conc = rr.collect_species_conc()

   rr.report_proteomics_correlation(scale='lin')
   rr.report_proteomics_correlation(scale='log')
   rr.report_protein_levels('Glucose')
   rr.plot_proteins('Glucose')

