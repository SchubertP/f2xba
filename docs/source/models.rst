Supported Models
================

The f2xba modeling framework provides support for a variety of model types, which are grouped into the following categories: enzyme constraint models (GECKO, ccFBA, MOMENT, and MOMENTmr); resource balanced constraint models (RBA); thermodynamics constraint FBA models (TFA); and mixed model types (TGECKO and TRBA). A detailed description of these model types can be found in the following sections, including implementation details.


GECKO
-----

A GECKO model couples enzyme-catalyzed reactions to protein requirements and places an upper limit on total protein (Sánchez et al., 2017). The f2xba program utilizes a genome-scale metabolic model (FBA model) and several configuration files as input to create an enzyme constraint GECKO model by executing a sequence of steps. Enzyme-catalyzed reactions are split into forward and reverse directions when reversible in the FBA model and per isoenzyme. For instance, the reversible reaction with the identifier ``R_FBA`` in the FBA model, catalyzed by two isoenzymes, is replaced by four reactions with identifiers ``R_FBA_iso1``, ``R_FBA_iso2``, ``R_FBA_iso1_REV`` and ``R_FBA_iso2_REV``. Non-negative protein concentration variables are added, with identifiers ``V_PC_\<uniprotID\>`` and for total modeled protein (``V_PC_total``). These variables are incorporated into the SBML reaction components in units of milligrams per gram of dry weight (mg/gDW). An upper bound on the total modeled protein is configured, and reaction fluxes are coupled to protein requirements according to the following general formulation:

.. math::

   \text{flux} ≤ \text{kcat} \cdot \text{n_AS} \cdot \text{avg_enz_sat} \cdot \frac{[P]}{\text{stoic} \cdot \text{MW}}

The *flux* is expressed in mmol/gDW, *kcat* denotes the turnover number in 1/h, *stoic* signifies the number of protein copies in the catalyzing enzyme, *n_AS* indicates the number of active sites in the enzyme, *\[P]* represents the protein concentration in mg/gDW, *MW* denotes the protein molecular weight in g/mol, and *avg_enz_sat* refers to the average enzyme saturation level. It is noteworthy that the average enzyme saturation can attain non-pysiological values greater than 1.0 to predict measured growth rates. In the context of the optimization problem, the inequality is replaced by an equality. Coupling factors :math:`\text{CC}_{ij}` for protein *i* and reaction *j* are introduced, according to the formula:

.. math::

   \text{CC}_{ij} = \frac{\text{stoic}_{ij} \cdot \text{MW}_i } {(\text{kcat}_{ij} \cdot \text{n_AS}_{ij} \cdot \text{avg_enz_sat})}

Consequently, protein mass balance constraints are incorporated, with identifiers ``C_prot_<uniprotID>`` and total modeled protein (``C_prot_pool``). These constraints are incorporated into the SBML species components. The individual and total protein mass balance constraints are configured using the formula: 

.. math::

   \text{C_prot}_i:   \sum_j{\text{CC}_{ij} \cdot \text{R}_j} = \text{V_PC}_i

.. math::

   \text{C_prot_pool}: \sum_i{\text{V_PC}_i} = \text{V_PC_total}

:math:`R_j` being the flux carried by reaction j. The coupling constraints are added to the SBML reaction components as products and reactants.


ccFBA
-----

The ccFBA model (cost constraint FBA) has been implemented according to the R-package sybilccFBA [1]_.Derivation of a ccFBA model from a GECKO model occurs through the removal of the redundancy provided by isoenzymes. ccFBA retains only these isoenzymes with the lowest protein "cost", respectively having the lowest coupling factor. This results in the creation of significantly smaller models that can exhibit equivalent predictions to those of GECKO models. However, these models demonstrate suboptimal performance when simulating gene deletion studies.


MOMENTmr
--------

The MOMENTmr model (MOMENT considering enzyme promiscuity) has been implemented according to the R-package sybilccFBA [1]_. The derivation of a MOMENTmr model from a GECKO model occurs by setting protein stoichiometry and the number of active sites to 1 in the formula for coupling coefficients. MOMENTmr models can be useful for extended models of organisms when detailed protein stoichiometry in enzyme composition is unknown or not relevant.


MOMENT
------

A MOMENT model (Adadi et al., 2012 <https://doi.org/10.1371/journal.pcbi.1002575>) employs default enzyme composition, akin to MOMENTmr, yet does not consider enzyme promiscuity. In instances where a protein is required to catalyze multiple reactions, this distinction becomes salient. GECKO, ccFBA, and MOMENTmr calculate the total cost of the protein across these reactions, whereas MOMENT only considers the cost of the reaction with the highest cost. Catalysis of the other reactions would not incur a cost to the model. While this model appears more straightforward, it is more intricate in its implementation. Additionally, during optimization, inequality signs on constraints must be configured.


RBA
---

RBA (resource balance constraint) models (`Goelzer et al., 2011 <https://doi.org/https://doi.org/10.1016/j.automatica.2011.02.038>`_) introduce process machines for processing of macromolecules and replacemnt of the fixed biomass reaction by concentration targets. The models are implemented using the formalism of RBApy (`Bulović et al., 2019 <https://doi.org/https://doi.org/10.1016/j.ymben.2019.06.001>`_), with adjustments that result in fully annotated models that are coded and stored in standardized SBML language.

A comparison of the RBA and GECKO models reveals notable distinctions in their implementation. In the RBA model, reactions are divided according to isoenzyme, and no splitting in forward and reverse directions is required. Reactions are coupled at the enzyme level, while in the GECKO model, they are coupled at the protein level. The RBA model contains enzyme concentration variables, while the GECKO model contains protein concentration variables. A one-to-one relationship exists between reactions and modeled enzymes, with enzyme identifiers being derived from reaction identifiers. Enzymes are composed of proteins, which need to be produced by process machines. Enzymes may also require cofactors, which need to be provided by the metabolic network. The expression of protein and RNA masses is carried out in units of average amino acids, while the capacity constraints present in different compartments are expressed in units of mmol amino acids per gram cell dry weight (mmolAA/gDW).

Enzyme efficiency constraints for the forward :math:`\text{C_EF}_i` and reverse :math:`\text{C_ER}_i` directions couple reaction fluxes :math:`\text{R}_i` [mmol/gDWh] with enzyme concentrations :math:`\text{V_EC}_i` [µmol/gDW] using apparent catalytic constants :math:`\text{kcat_app}_i` [:math:`h^{-1}`].

.. math::

  \text{C_EF}_i: \text{R}_i - \text{kcat_app}_i \cdot \text{V_EC}_i \cdot 0.001 \leq 0 
  \text{C_ER}_i: - \text{R}_i - \text{kcat_app}_i \cdot \text{V_EC}_i \cdot 0.001  \leq 0 

The apparent catalytic constants are determined from the turnover numbers :math:`\text{kcat}_i` [:math:`s^{-1}`], the number of active sites :math:`\text{n_AS}_i`, and the average enzyme saturation :math:`\text{avg_enz_sat}`, which is used for all reactions. 

.. math::

  \text{kcat_app}_i = \text{kcat}_i \cdot \text{n_AS}_i  \cdot \text{avg_enz_sat}  \cdot 3600

Macromolecular production reactions are incorporated to synthesize explicitly modeled proteins (e.g., ``R_PROD_b2907``), balance dummy proteins (e.g., ``R_PROD_dummy_protein_im``), and total DNA (``R_PROD_dna``), total mRNA (``R_PROD_mrna``), individual tRNAs (e.g., ``R_PROD_trnaala``) and individual rRNAs (e.g., ``R_PROD_rRNA_16S``), using metabolites produced in the metabolic reaction network. Process machine capacity constraints (e.g., ``C_PMC_pm_translation``) couple these macromolecular production reactions to process machine concentrations (e.g., ``V_PMC_pm_translation``) in the same way as metabolic reactions are coupled to enzyme requirements. Macromolecular degradation reactions are added in a similar way (e.g., ``R_DEGR_mrna``).

It is important to note that macromolecules undergo dilution due to cellular growth. Mass balance constraints for macromolecules (e.g., ``MM_b2907``) ensure that macromolecule production, degradation, and dilution are balanced. The dilution of modeled proteins and ribosomal RNAs is governed by enzyme concentration (e.g., ``V_EC_FBA_iso1``) and process machine concentration (e.g., ``V_PCM_pm_translation``) variables. The dilution of other macromolecules due to growth is managed by specific target concentration variables (e.g., ``V_TMMC_mrna``). The variable ``V_TSMC`` [µmol/gDW] controls the growth dilution of selected metabolites, mainly derived from the biomass pseudo reaction of the FBA model. Compartmental capacity limits are controlled by the variable ``V_TCD`` [mmolAA/gDW]. Density constraints (e.g., ``C_D_im``) serve to regulate the concentrations of enzymes, processes, machines, and macromolecular targets, ensuring that these concentrations do not exceed the capacities of the respective compartments.


TFA
---

The TFA (thermodymanics constraint FBA model) (`Henry et al., 2007 <https://doi.org/10.1529/biophysj.106.093138>`_) has been implemented based on the pyTFA package (`Salvy et al., 2019 <https://doi.org/10.1093/bioinformatics/bty499>`_), with adjustments.

Variables and Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^

The following paragraphs detail the variables and constraints that have been incorporated into the genome-scale metabolic model. To illustrate this, the reaction of fructose-bisphosphate aldolase is examined, which is implemented as a reversible reaction in the iML1515 model with the identifier R_FBA: 'fdp_c -> dhap_c + g3p_c'.  During the extension of the model, the reaction is rendered irreversible, designated as ``R_FBA``: 'fdp_c => dhap_c + g3p_c', and a new reaction catalyzing the reverse direction is incorporated. This reverse reaction is designated as ``R_FBA_REV``: 'dhap_c + g3p_c => fdp_c'. It is noteworthy that the reverse reaction is not included for reactions that have been configured as irreversible in the original model.

Two additional variables are incorporated to represent the transformed Gibbs energy of reaction, designated as ``V_DRG_FBA`` and ``V_DRG0_FBA``, respectively. The former is assigned unlimited bounds, while the latter is constrained to the calculated standard transformed Gibbs energy of reaction plus or minus the estimation error.

The log concentration variables, designated as ``V_LC_fdp_c``, ``V_LC_dhap_c``, and ``V_LC_g3p_c``, are incorporated with default bounds as per compartmental configuration, and the concentrations are expressed in units of mol/L. The calculation formula for the transformed Gibbs energy of reaction is implemented by the constraint ``C_DRG_FBA``: 'V_DRG_FBA = V_DRG0_FBA - 2.48 V_LC_fdp_c + 2.48 V_LC_dhap_c + 2.48 V_LC_g3p_c' (RT = 2.48 kJ/mol).

Two binary variables (values 0 or 1) designated ``V_FU_FBA`` ("forward use") and ``V_RU_FBA`` ("reverse use") are introduced to couple the transformed Gibbs energy of reaction to the flux direction. The implementation of a "simultaneous use" constraint, denoted as ``C_SU_FBA``, ensures that only one of the use variables can take the value "1" This is expressed as "V_FU_FBA + V_RU_FBA ≤ 1." Two Gibbs energy coupling constraints, designated as ``C_GFC_FBA`` and ``C_GRC_FBA``, couple the forward and reverse use variables to the transformed Gibbs energy of reaction with inequalities ‘V_DRG_FBA ≤ 999.99 - 1000 V_FU_FBA’, thereby forcing V_DRG_FBA ≤ 0.01 kJ/mol when V_FU_FBA is active, and 'V_DRG_FBA ≥ 1000 V_RU_FBA - 999.99', forcing V_DRG_FBA ≥ 0.01 kJ/mol when V_RU_FBA is active. In a similar fashion, reactions fluxes in forward and reverse direction are coupled to the forward and reverse use variables via the flux coupling constraints ``C_FFC_FBA``: R_FBA ≤ 1000 V_FU_FBA and ``C_FRC_FBA``: R_FBA_REV ≤ 1000 V_FB_FBA. This configuration can be readily verified by exporting the TFA model to the '.xlsx' format.

Calculation Details 
^^^^^^^^^^^^^^^^^^^

Thermodynamic constraints couple reaction flux directionality with Gibbs energy of reaction. The transformed Gibbs energy values employed in this context are derived through transformations with respect to compartmental pH and ionic strength at the default temperature of 298.15 K (25˚C). Negative values of the transformed Gibbs energy of reaction will drive the reaction in the forward direction, while positive values will drive it in the reverse direction. The incorporation of thermodynamic constraints into a model is a straightforward process that necessitates minimal configuration input. However, it should be noted that the underlying calculations are of a highly complex nature. The f2xba model is aligned with the formulation implemented in the pyTFA package (`Salvy et al., 2019 <https://doi.org/10.1093/bioinformatics/bty499>`_), with the formulas being based on the book by Alberty [2]_.

The natural log of metabolite concentrations is a variable in the optimization problem, with lower and upper bounds defined in the TFA configuration file and potentially further limited prior to optimization. The factor of gas constant times temperature, 'RT', used in subsequent equations, is evaluated to 2.48 kJ/mol at T = 298.15 K.

The transformed Gibbs energy of reaction, denoted by :math:`\Delta_r \text G^{'}`, is calculated from the standard transformed Gibbs energy of reaction, denoted by :math:`\Delta_r \text G^{'0}`, the metabolite concentrations :math:`c_i` [mol/L], and the stoichiometric coefficients :math:`\nu_i` of reactants (negative) and products (positive) using the following equation 4.5-10 [2]_:

.. math::
  
   \Delta_r \text G^{'} = \Delta_r \text G^{'0} + \text{RT} \sum_i {\nu_i \ln c_i}

The group contribution method involves the decomposition of a molecule into cues (groups) for which the Gibbs free energy of formation has been estimated with high confidence. During a reaction, only some of the cues of reactants and products, the net cues, undergo modification, while the majority of the other cues remain unmodified. The TD database contains estimated errors for each of the cues. The estimated errors of the net cues, denoted by :math:`\text{cue_est_error}_j`, are then utilized to ascertain the estimation error for the Gibbs energy of reaction. This estimation error is subsequently employed to establish the bounds of the variable :math:`\Delta_r \text G^{'0}`. The calculation is performed as per pyTFA:

.. math::

   \text{estimation_error} = \sqrt {\sum_j (\nu_j \, \text{cue_est_error}_j)^2}

The standard-transformed Gibbs energy of reaction, denoted by :math:`\Delta_r \text G^{'0}`, is derived from the standard-transformed Gibbs energies of formation, denoted by :math:`\Delta_f \text G_i^{'0}`, of the reactants and products. This derivation is expressed through the following equation 4.4-2 [2]_:

.. math::
 
   \Delta_r \text G^{'0} = \sum_i \nu_i \Delta_f \text G_i^{'0}

In the context of a transport process, it is imperative to incorporate electrical work terms, which are calculated from the membrane potentials, denoted by the symbol :math:`\Delta \varphi_{sd}` (destination minus source potential), and the transported charges, denoted by the symbol :math:`z_{sd}` (source to destination compartment), which can assume positive or negative values. :math:`F` is the Faraday constant (96.485 kJ mol-1 V-1). To derive the equation, we have consulted the work of Jol (`Jol et al., 2010 <https://doi.org/10.1016/j.bpj.2010.09.043>`_).

.. math::

  \Delta_r \text G^{'0} = \sum_i \nu_i \Delta_f \text G_i^{'0} + F \Delta \sum_{sd} \varphi_{ds} z_{sd}

The reaction network implemented by the genome-scale metabolic model consists of biochemical reactions and biochemical reactants (metabolites). At the molecular level, a biochemical reactant, such as ATP, can be regarded as a group (pseudoisomer group) of related chemical species in different protonation states and different complexations with metal ions, such as :math:`\text{ATP}^{4-}`, :math:`\text{HATP}^{3-}`, :math:`\text{H}_2\text{ATP}^{2-}`, :math:`\text{MgATP}^{2-}` or :math:`\text{MgHATP}^-`. The f2xba modeling framework considers different protonation states, but not different complexations with metal ions. 

The Gibbs energy of formation, denoted as :math:`\Delta_f \text G^{'0}(I)`, of a biochemical reactant can be determined from the Gibbs energy of formation of the least protonated chemical species, denoted as :math:`\Delta_f \text {G1}^{'0}(I)`, in the pseudoisomeric group and the contribution of the other chemical species in the pseudoisomeric group, which are considered by the binding polynomial, denoted as :math:`\text P(I)`. It is imperative to note that these values are contingent on the isomeric strength :math:`I` [mol/L] of the compartment, as elucidated in equations 4.5-6 [2]_.
  
.. math::

  \Delta_f \text G^{'0}(I) = \Delta_f \text{G1}^{'0}(I) - \text{RT} \ln(\text P(I))

The determination of the least protonated chemical species is contingent upon the compartmental pH, as defined in the TFA configuration file, in conjunction with the :math:`\text{pKa}_j` values and the electrical charge extracted from the pertinent TD data record. The standard-transformed Gibbs energy of formation for the least protonated species, denoted as :math:`\Delta_f \text{G1}^{'0}`, is derived from the standard Gibbs energy of formation, denoted as :math:`\Delta_f \text G^{0}`, extracted from the corresponding TD database record. The latter is transformed to the compartmental :math:`pH` using equation 4.10-12 [2]_:

.. math::

   \Delta_f \text {G1}^{'0} = \Delta_f \text G^0 + RT \ln{(10)} \sum_j {\text{pKa}_j}

The standard transformed Gibbs energy of formation for this least protonated state, denoted as :math:`\Delta_f \text{G1}^{'0}(I)` at a given ionic strength, is calculated from its value at zero ionic strength, denoted as :math:`\Delta_f \text{G1}^{'0}`, and adjustments with respect to the energy contribution by the number of protons, denoted as :math:`nH`, in the structure and electrical charge, denoted as :math:`z`. Ionic strength :math:`I`, defined as the measure of ion concentration, exerts a significant influence on the activity coefficients employed in equilibrium equations by means of shielding charges. This adjustment is achieved through the utilization of equation 4.4-10 [2]_, which is founded on the extended Debye-Hueckel equation with constants :math:`A = 0.51065 \sqrt\frac{l}{mol}` and :math:`B = 1.6\sqrt\frac{l}{mol}`.

.. math::

   \Delta_f \text{G1}^{'0}(I) = \Delta_f \text{G1}^{'0} + \text{RT} \cdot \ln(10) \cdot \text{nH} \cdot \text{pH}  - \text{RT} \cdot \ln(10) \cdot (z^2 - \text{nH}) \cdot \frac{A \sqrt I}{1 + B \sqrt I}

Ionic strength-adjusted acid dissociation constants, denoted by the symbol :math:`\text{pKa}^{'}(I)`, are essential for determining the binding polynomial. These constants can be calculated from the corresponding :math:`\text{pKa}` values using the following equation 4.10-11 [2]_:

.. math::

   \text{pKa}^{'}(I) = \text{pKa} - \sum_j \nu_j z_j^2 \cdot \frac{A \sqrt I}{1 + B \sqrt I}

The binding polynomial, denoted by :math:`P(I)`, which accounts for the energy contribution of the other chemical species in the pseudoisomer group, is calculated from the equilibrium constants :math:`K_x`, which relate to the ionic strength-adjusted acid dissociation constants (:math:`K_x = 10^{-\text{pKa}_x^{'}}`), with :math:`K_1` having the smallest value (highest pKa) and the compartmental proton concentration :math:`[H^+] = 10^{-pH}`, using equation 4.5-7 [2]_.

.. math::

  P(I) = 1 + \frac {[H^+]} {K_1} + \frac {[H^+]^2} {K_1 \cdot K_2} + \dotsb

The mean number of bound hydrogens, denoted by :math:`\text{avg_h_binding}`, in the pseudoisomeric group can be calculated from the binding polynomial, as outlined in equation 1.3-9 [2]_.

.. math::

  \text{avg_h_binding} = \frac {[H+]} {\text P}  \frac {d \text P} {d[H+]} = \frac{[H+]/K_1 + 2[H+]^2/K_1K_2 + \dotsb } {\text P(I)}


TGECKO
------

The TGECKO model, a thermodynamics constraint GECKO model, is a combination of a TFA and a GECKO model. Similarly, TccFBA, TMOMENTmr, and TMOMENT can be constructed.


TRBA
----

The TRBA model, a thermodynamics constraint RBA model, is a combination of a TFA and a RBA model.


References
----------

.. [1] 
   Desouki, A. A. (2015). sybilccFBA: Cost Constrained FLux Balance Analysis: MetabOlic Modeling with ENzyme kineTics (MOMENT).  
   In CRAN. https://cran.r-project.org/web/packages/sybilccFBA/index.html

.. [2]
   Alberty, R. A. (2003). Thermodynamics of Biochemical Reactions. Massachusetts Institute of Technology Press, Cambridge, MA.
