Supported Models
================

The f2xba modeling framework provides support for a variety of model types, which are grouped into the following categories: enzyme constraint models (GECKO, ccFBA, MOMENT, and MOMENTmr); resource balanced constraint models (RBA); thermodynamics constraint FBA models (TFA); and mixed model types (TGECKO and TRBA). A detailed description of these model types can be found in the following sections, including implementation details.


GECKO
-----

A GECKO model couples enzyme-catalyzed reactions to protein requirements and places an upper limit on total protein (Sánchez et al., 2017). The f2xba program utilizes a genome-scale metabolic model (FBA model) and several configuration files as input to create an enzyme constraint GECKO model by executing a sequence of steps. Enzyme-catalyzed reactions are split into forward and reverse directions when reversible in the FBA model and per isoenzyme. For instance, the reversible reaction with the identifier ``R_FBA`` in the FBA model, catalyzed by two isoenzymes, is replaced by four reactions with identifiers ``R_FBA_iso1``, ``R_FBA_iso2``, ``R_FBA_iso1_REV`` and ``R_FBA_iso2_REV``. Non-negative protein concentration variables are added, with identifiers ``V_PC_<uniprotID>`` and for total modeled protein (``V_PC_total``). These variables are incorporated into the SBML reaction components in units of milligrams per gram of dry weight (mg/gDW). An upper bound on the total modeled protein is configured, and reaction fluxes are coupled to protein requirements according to the following general formulation:

.. math::

   flux \leq kcat \cdot n\_AS \cdot avg\_enz\_sat \cdot \frac{[P]}{stoic \cdot MW}

The *flux* is expressed in mmol/gDW, *kcat* denotes the turnover number in 1/h, *stoic* signifies the number of protein copies in the catalyzing enzyme, *n_AS* indicates the number of active sites in the enzyme, *\[P]* represents the protein concentration in mg/gDW, *MW* denotes the protein molecular weight in g/mol, and *avg_enz_sat* refers to the average enzyme saturation level. A GECKO model not properly parametrized may require an average enzyme saturation with a non-pysiological value greater than 1.0 to predict measured growth rates. In the context of the optimization problem, the inequality is replaced by an equality. Coupling factors :math:`CC_{ij}` for protein *i* and reaction *j* are introduced, according to the formula:

.. math::

   CC_{ij} = \frac{stoic_{ij} \cdot MW_i } {(kcat_{ij} \cdot n\_AS_{ij} \cdot avg\_enz\_sat)}

Consequently, protein mass balance constraints are incorporated, with identifiers ``C_prot_<uniprotID>`` and total modeled protein (``C_prot_pool``). These constraints are incorporated into the SBML species components. The individual and total protein mass balance constraints are configured using the formula: 

.. math::

   C\_prot_i: \sum_j{CC_{ij} \cdot R_j} = V\_PC_i 
   
   C\_prot\_pool: \sum_i{V\_PC_i} = V\_PC\_total


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

RBA (resource balance constraint) models (`Goelzer et al., 2011 <https://doi.org/https://doi.org/10.1016/j.automatica.2011.02.038>`_) introduce process machines for processing of macromolecules and replacement of the fixed biomass reaction by concentration targets. The models are implemented using the formalism of RBApy (`Bulović et al., 2019 <https://doi.org/https://doi.org/10.1016/j.ymben.2019.06.001>`_), with adjustments that result in fully annotated models that are coded and stored in standardized SBML language.

A comparison of the RBA and GECKO models reveals notable distinctions in their implementation. In the RBA model, reactions are divided according to isoenzyme, and no splitting in forward and reverse directions is required. Reactions are coupled at the enzyme level, while in the GECKO model, they are coupled at the protein level. The RBA model contains enzyme concentration variables, while the GECKO model contains protein concentration variables. A one-to-one relationship exists between reactions and modeled enzymes, with enzyme identifiers being derived from reaction identifiers. Enzymes are composed of proteins, which need to be produced by process machines. Enzymes may also require cofactors, which need to be provided by the metabolic network. The expression of protein and RNA masses is carried out in units of average amino acids, while the capacity constraints present in different compartments are expressed in units of mmol amino acids per gram cell dry weight (mmolAA/gDW).

Enzyme efficiency constraints for the forward :math:`C\_EF_i` and reverse :math:`C\_ER_i` directions couple reaction fluxes :math:`R_i` [mmol/gDWh] with enzyme concentrations :math:`V\_EC_i` [µmol/gDW] using apparent catalytic constants :math:`kcat\_app_i` [:math:`h^{-1}`].

.. math::

   C\_EF_i: R_i - kcat\_app_i \cdot V\_EC_i \cdot 0.001 \leq 0
   
   C\_ER_i: - R_i - kcat\_app_i \cdot V\_EC_i \cdot 0.001  \leq 0 


The apparent catalytic constants are determined from the turnover numbers :math:`kcat_i` [:math:`s^{-1}`], the number of active sites :math:`n\_AS_i`, and the average enzyme saturation :math:`avg\_enz\_sat`, which is used for all reactions. 

.. math::

   kcat\_app_i = kcat_i \cdot n\_AS_i  \cdot avg\_enz\_sat  \cdot 3600

Macromolecular production reactions are incorporated to synthesize explicitly modeled proteins (e.g., ``R_PROD_b2907``), balance dummy proteins (e.g., ``R_PROD_dummy_protein_im``), and total DNA (``R_PROD_dna``), total mRNA (``R_PROD_mrna``), individual tRNAs (e.g., ``R_PROD_trnaala``) and individual rRNAs (e.g., ``R_PROD_rRNA_16S``), using metabolites produced in the metabolic reaction network. Process machine capacity constraints (e.g., ``C_PMC_pm_translation``) couple these macromolecular production reactions to process machine concentrations (e.g., ``V_PMC_pm_translation``) in the same way as metabolic reactions are coupled to enzyme requirements. Macromolecular degradation reactions are added in a similar way (e.g., ``R_DEGR_mrna``).

It is important to note that macromolecules undergo dilution due to cellular growth. Mass balance constraints for macromolecules (e.g., ``MM_b2907``) ensure that macromolecule production, degradation, and dilution are balanced. The dilution of modeled proteins and ribosomal RNAs is governed by enzyme concentration (e.g., ``V_EC_FBA_iso1``) and process machine concentration (e.g., ``V_PCM_pm_translation``) variables. The dilution of other macromolecules due to growth is managed by specific target concentration variables (e.g., ``V_TMMC_mrna``). The variable ``V_TSMC`` [µmol/gDW] controls the growth dilution of selected metabolites, mainly derived from the biomass pseudo reaction of the FBA model. Compartmental capacity limits are controlled by the variable ``V_TCD`` [mmolAA/gDW]. Density constraints (e.g., ``C_D_im``) serve to regulate the concentrations of enzymes, processes, machines, and macromolecular targets, ensuring that these concentrations do not exceed the capacities of the respective compartments.


TFA
---

The TFA (thermodymanics constraint FBA model) (`Henry et al., 2007 <https://doi.org/10.1529/biophysj.106.093138>`_) has been implemented based on the pyTFA package (`Salvy et al., 2019 <https://doi.org/10.1093/bioinformatics/bty499>`_), with adjustments.

Variables and Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^

The following paragraphs detail the variables and constraints that have been incorporated into the genome-scale metabolic model. To illustrate this, the reaction of fructose-bisphosphate aldolase is examined, which is implemented as a reversible reaction in the iML1515 model with the identifier R_FBA: 'fdp_c -> dhap_c + g3p_c'.  During the extension of the model, the reaction is rendered irreversible, designated as ``R_FBA``: 'fdp_c => dhap_c + g3p_c', and a new reaction catalyzing the reverse direction is incorporated. This reverse reaction is designated as ``R_FBA_REV``: 'dhap_c + g3p_c => fdp_c'. It is noteworthy that the reverse reaction is not included for reactions that have been configured as irreversible in the original model.

Two additional variables are incorporated to represent the standard transformed and the transformed Gibbs energy of reaction, designated as ``V_DRG0_FBA`` and ``V_DRG_FBA``, respectively. The former is constrained to the calculated standard transformed Gibbs energy of reaction plus and minus the estimation error, while the latter is unconstrained.

The log10 concentration variables, designated as ``V_LC_fdp_c``, ``V_LC_dhap_c``, and ``V_LC_g3p_c``, are incorporated with default bounds as per compartmental configuration, and the concentrations are expressed in units of mol/L. The calculation formula for the transformed Gibbs energy of reaction is implemented by the constraint ``C_DRG_FBA``: 'V_DRG_FBA = V_DRG0_FBA - 5.71 V_LC_fdp_c + 5.71 V_LC_dhap_c + 5.71 V_LC_g3p_c', with RT ln(10) = 5.71 kJ/mol and stoichiometric coefficients of -1 and +1.

Inline with the implementation in pyTFA, two binary variables (values 0 or 1) designated ``V_FU_FBA`` ("forward use") and ``V_RU_FBA`` ("reverse use") are introduced to couple the transformed Gibbs energy of reaction to the flux direction. The implementation of a "simultaneous use" constraint, denoted as ``C_SU_FBA``, ensures that only one of the use variables can take the value "1." This is expressed as 'V_FU_FBA + V_RU_FBA ≤ 1.0'. Two Gibbs energy coupling constraints, designated as ``C_GFC_FBA`` and ``C_GRC_FBA``, couple the forward and reverse use variables to the transformed Gibbs energy of reaction with inequalities 'V_DRG_FBA ≤ 999.99 - 1000 V_FU_FBA', thereby forcing V_DRG_FBA ≤ -0.01 kJ/mol when V_FU_FBA is active, and 'V_DRG_FBA ≥ 1000 V_RU_FBA - 999.99', forcing V_DRG_FBA ≥ 0.01 kJ/mol when V_RU_FBA is active. In a similar fashion, reaction fluxes in forward and reverse direction are coupled to the forward and reverse use variables via the flux coupling constraints ``C_FFC_FBA``: 'R_FBA ≤ 1000 V_FU_FBA' and ``C_FRC_FBA``: 'R_FBA_REV ≤ 1000 V_RU_FBA'. This configuration can be readily verified by exporting the TFA model to the '.xlsx' format.

Calculation Details 
^^^^^^^^^^^^^^^^^^^

Thermodynamic constraints couple reaction flux directionality with the Gibbs energy of reaction. The transformed Gibbs energy values employed in this context are derived through transformations with respect to compartmental pH and ionic strength at the default temperature of 298.15 K (25˚C). Negative values of the transformed Gibbs energy of reaction will drive the reaction in the forward direction, while positive values will drive it in the reverse direction. The incorporation of thermodynamic constraints into a model is a straightforward process that necessitates minimal configuration input. However, it should be noted that the underlying calculations are of a highly complex nature. The f2xba model is aligned with the formulation implemented in the pyTFA package (`Salvy et al., 2019 <https://doi.org/10.1093/bioinformatics/bty499>`_), with adjustments. Standard transformed Gibbs energies of formation are calculated using formulations given in the book by Alberty, 2003 [2]_. Reference has been taken from eQuilibrator-api (`Beber et al., 2022 <https://doi.org/10.1093/nar/gkab1106>`_) and literature (`Haraldsdóttir et al., 2012 <https://doi.org/https://doi.org/10.1016/j.bpj.2012.02.032>`_; `Jol et al., 2010 <https://doi.org/10.1016/j.bpj.2010.09.043>`_) to determine transformed Gibbs reaction energies.

The log10 of metabolite concentrations is a variable in the optimization problem, with lower and upper bounds defined in the TFA configuration file and potentially further constraint prior to optimization. The factor of gas constant times temperature, 'RT', used in subsequent equations, is evaluated to 2.48 kJ/mol at T = 298.15 K.

The transformed Gibbs energy of reaction, denoted by :math:`\Delta_r G^{'}`, is calculated from the standard transformed Gibbs energy of reaction, denoted by :math:`\Delta_r G^{'0}`, the metabolite concentrations :math:`c_i` [mol/L], and the stoichiometric coefficients :math:`\nu_i` of reactants (negative) and products (positive) using the equation (Alberty, 2003 [2]_, Eq. 4.5-10):

.. math::

   \Delta_r G^{'} = \Delta_r G^{'0} + RT \ln(10) \sum_i {\nu_i \log c_i}

The group contribution method (`Mavrovouniotis et al. 1990 <https://doi.org/10.1002/bit.260361013>`_) involves the decomposition of a molecule into functional groups for which the Gibbs free energy of formation has been estimated with high confidence. During a reaction only some of the functional groups of reactants and products, the net groups, undergo modification, while the majority of the other groups remain unmodified. The TD database (`Jankowski et al. 2008 <https://doi.org/10.1529/biophysj.107.124784>`_) contains estimated errors for each of the functional groups. The estimated errors of the net groups, denoted by :math:`cue\_est\_error_j`, are then utilized to ascertain the estimation error for the Gibbs energy of reaction. This estimation error is subsequently employed to establish the bounds of the variable :math:`\Delta_r G^{'0}`. The calculation is performed as per pyTFA:

.. math::

   est\_error = \sqrt {\sum_j (\nu_j \, cue\_est\_error_j)^2}

The standard transformed Gibbs energy of reaction, denoted by :math:`\Delta_r G^{'0}`, is derived from the standard transformed Gibbs energies of formation, denoted by :math:`\Delta_f G_i^{'0}`, of the reactants and products. For non-transport related reactions we used the equation (Alberty, 2003 [2]_, Eq. 4.4-2):

.. math::

   \Delta_r G^{'0} = \sum_i \nu_i \Delta_f G_i^{'0}

Standard transformed Gibbs reaction energies of transport related reactions are calculated by adding two terms. One term accounts for the disparate proton potentials on either side of the membrane. The calculation of this term is based on the net number of hydrogen atoms, including protons, that are transported between the source and destination compartment :math:`nH_{sd}`, and the respective compartmental pH values. The second term encompasses the electrical work involved and is contingent on the membrane potential difference, denoted by the symbol :math:`\Delta \varphi_{sd}` (destination minus source electrical potential), and the net transported charge, denoted by the symbol :math:`z_{sd}` (source to destination compartment). This formulation supports transport reactions, where transported components get modified during transport, e.g. PTS systems, as well as transport reactions involving several membranes, e.g. export of a metabolite across the external membrane driven by proton translocation across the inner membrane. The Faraday constant, denoted by :math:`F` has the value of 96.485 :math:`\frac{kJ}{mol \cdot V}`. The equation was derived by consulting literature (`Haraldsdóttir et al., 2012 <https://doi.org/https://doi.org/10.1016/j.bpj.2012.02.032>`_; `Jol et al., 2010 <https://doi.org/10.1016/j.bpj.2010.09.043>`_) and the implementation in equilibrator-api (`Beber et al., 2022 <https://doi.org/10.1093/nar/gkab1106>`_).

.. math::

   \Delta_r G^{'0} = \sum_i \nu_i \Delta_f G_i^{'0} - RT \ln(10) \sum_{sd} {nH_{sd} (-pH_s + pH_d)} + F \sum_{sd} \Delta \varphi_{sd} z_{sd}

The reaction network implemented by the genome-scale metabolic model consists of biochemical reactions and biochemical reactants (metabolites). At the molecular level a biochemical reactant such as ATP can be regarded as a group (pseudoisomer group) of related chemical species in different protonation states and different complexations with metal ions, such as such as ATP :math:`^{4-}`, HATP :math:`^{3-}`, :math:`H_2` ATP :math:`^{2-}`, MgATP :math:`^{2-}` or MgHATP :math:`^-`. The f2xba modeling framework considers different protonation states, but not different complexations with metal ions. Note: eQuilibrator (`Noor et al., 2013 <https://doi.org/10.1371/journal.pcbi.1003098>`_) is a more elaborate thermodynamics calculation application, considering as well compartmental magnesium concentration levels and employing a differently constructed thermodynamics database. The user is free to determine transformed Gibbs reaction energies and related errors using eQuilibrator or an alternative application and supply these values during TFA model construction, to replace the functionlity implemented in f2xba.

The standard transformed Gibbs energy of formation of a biochemical reactant, denoted as :math:`\Delta_f G^{'0}(I)`, can be determined from the standard transformed Gibbs energy of formation of the least protonated chemical species in the pseudoisomeric group, denoted as :math:`\Delta_f {G1}^{'0}(I)`, and the contribution of the other chemical species in the pseudoisomeric group, which are considered by the binding polynomial, denoted as :math:`P(I)`. It is imperative to note that these values are contingent on the isomeric strength :math:`I` [mol/L] of the compartment, as elucidated in equation (Alberty, 2003 [2]_, Eq. 4.5-6).

.. math::

   \Delta_f G^{'0}(I) = \Delta_f G1^{'0}(I) - RT \ln(P(I))

The standard transformed Gibbs energy of formation of the least protonated state at zero ionic strength, denoted as :math:`\Delta_f G^{'0}`, is calculated from the standard Gibbs energy of formation, denoted as :math:`\Delta_f G^{0}`, extracted from the thermodynamics database and  protonation steps up to pH 9.0. The protonation level of the related species in the TD database is determined from its list of pKa values, its electrical charge in standard conditions and the electrical charge at the most protonated state. Different to pyTFA, we consider positively charged functional groups, e.g. amino groups, when determining the charge for the most protonated state. The transformation from the standard state to the protonation state at pH 9.0 is determind using the equation (Alberty, 2003 [2]_, Eq. 4.10-12):

.. math::

   \Delta_f G1^{'0} = \Delta_f G^0 + RT \ln{(10)} \sum_j {pKa_j}

The standard transformed Gibbs energy of formation for this least protonated state at a given ionic strength, denoted as :math:`\Delta_f G1^{'0}(I)`, is calculated from its value at zero ionic strength and adjusted with respect to the energy contribution by the number of protons in the structure, denoted as :math:`nH`, and electrical charge, denoted as :math:`z`. Ionic strength :math:`I`, defined as the measure of ion concentration, exerts a significant influence on the activity coefficients employed in equilibrium equations by means of shielding charges. This adjustment is achieved through the utilization of below equation (Alberty, 2003 [2]_, Eq. 4.4-10), which is founded on the extended Debye-Hueckel equation with constants :math:`A = 0.51065 \sqrt\frac{l}{mol}` and :math:`B = 1.6\sqrt\frac{l}{mol}`.

.. math::

   \Delta_f G1^{'0}(I) = \Delta_f G1^{'0} + RT \cdot \ln(10) \cdot nH \cdot pH - RT \cdot \ln(10) \cdot (z^2 - nH) \cdot \frac{A \sqrt I}{1 + B \sqrt I}

Ionic strength-adjusted acid dissociation constants, denoted by the symbol :math:`pKa^{'}(I)`, are essential for determining the binding polynomial. These constants can be calculated from the corresponding :math:`pKa` values using the following equation (Alberty, 2003 [2]_, Eq. 4.10-11 ):

.. math::

   pKa^{'}(I) = pKa - \sum_j \nu_j z_j^2 \cdot \frac{A \sqrt I}{1 + B \sqrt I}

The binding polynomial, denoted by :math:`P(I)`, which accounts for the energy contribution of the other chemical species in the pseudoisomer group, is calculated from the equilibrium constants :math:`K_x`, which relate to the ionic strength-adjusted acid dissociation constants (:math:`K_x = 10^{-pKa_x^{'}}`), with :math:`K_1` having the smallest value (highest pKa) and the compartmental proton concentration :math:`[H^+] = 10^{-pH}`, using equation (Alberty, 2003 [2]_, Eq. 4.5-7).

.. math::

   P(I) = 1 + \frac {[H^+]} {K_1} + \frac {[H^+]^2} {K_1 \cdot K_2} + \dotsb

The mean number of bound hydrogens, denoted by :math:`avg\_h\_binding`, in the pseudoisomeric group can be calculated from the binding polynomial, as outlined in equation 1.3-9 [2]_.

.. math::

   avg\_h\_binding = \frac {[H+]} {P}  \frac {dP} {d[H+]} = \frac{\frac{[H+]}{K_1} + \frac{2[H+]^2}{K_1K_2} + \dotsb } {P(I)}


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

