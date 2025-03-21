Configuration Data
==================

This chapter provides information regarding the configuration tables utilized by the f2xba modeling framework. It is recommended that the user proceeds through the tutorials, in which these configuration files are established and employed. Tables are stored in '.xlsx' spreadsheets files with sheet names corresponding to the table names.

Table columns with default or "None" values are not required to be included in the tables. Additional columns are permissible, provided they do not conflict with the specified column headings. A beneficial column could be a "notes" column, which would allow for the annotation of data for future reference.

XBA configuration
-----------------

The XBA configuration file, used in configuration of the XbaModel instance, contains several tables, which are described in the following.


Table ``general``
^^^^^^^^^^^^^^^^^

The table labeled ``general`` contains the primary parameters for extended model configuration. The default turnover numbers can be configured for metabolic and transport reactions. Furthermore, references to online databases, such as the UniProt protein database (``organism_id``), the NCBI genome database (``chromosome2accids``), and the Biocyc organism database (``biocyc_org_prefix``), can be defined as per requirement. Furthermore, the provision of file names to configure reaction specific turnover numbers (``kcats_fname``) and enzyme compositions  (``enzyme_comp_fname``) is permitted, with the objective of replacing default values. Setting the ``cofactor_flag`` to "True" enables the incorporation of cofactors extracted from UniProt entries into enzyme compositions. Finally, a bulk mapping configuration table (``bulk_mappings_fname``) can be used for large-scale model reconfigurations, e.g., for remapping of identifiers used in the model.

.. csv-table:: ``general``: general configuration data
   :header: "Column", "Contents", "Example"

   "parameter", "parameter name", "thermo_data_fname"
   "value", "value", "data/thermo_data.thermodb"

Following parameters can be configured with this table:

.. csv-table:: 
   :header: "Parameter", "Description", "Example"

   "default_metabolic_kcat", "default for metabolic reaction turnover numbers in s-1 (default: 12.5)", 12.5
   "default_transporter_kcat", "default for transporter reaction turnover number in s-1 (default: 50.0)", 100.0
   "organism_dir", "directory with organism specific data files", "data"
   "organism_id", "taxonomic identifier", 83333
   "chromosome2accids", "NCBI genome accession identifiers", "chromosome=U00096.3"
   "biocyc_org_prefix", "BioCyc organism identifier", "ecoli"
   "kcats_fname", "reaction specific turnover numbers", "data/iML1515_predicted_fit_GECKO_kcats.xlsx"
   "enzyme_comp_fname", "enzyme composition", "data/iML1515_enzyme_composition_updated.xlsx"
   "bulk_mappings_fname", " bulk mapping configuration", "data_configs/iJN678_bulk_mappings.xlsx"
   "cofactor_flag", "cofactor use in RBA enzymes (default: False)", True

Notes: RBA models require "chromosome2accids" to reference all chromosomes and plasmids that contain genes used in the model, e.g. "Chr_I=BK006935.2, Chr_II=BK006936.2, Chr_III=BK006937.2, ...".


Table ``modify_attributes``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The table designated ``modify_attributes`` allows the modification of model component attributes, e.g. modifying flux bounds on reaction components.  

.. csv-table:: ``modify_attributes``: modify model component attributes
   :header: "Column", "Contents", "Example"

   "id", "component identifier", "R_ALAt2pp"
   "component", "component type or None", "reaction"
   "attribute", "attribute name", "fbc_lower_bound"
   "value", "new value", "cobra_0_bound"

The "attribute" identifier and the type of "value" depends on the "component" type. Attributes of following components can be modified: "modelAttrs", "compartment", "species", "reaction", "gp", "protein", "enzyme", "ncbi" and "uniprot". Some examples are provided in below table:

.. csv-table:: examples for table entries
   :header: "id", "component", "attribute", "value"

   "R_ALAt2pp", "reaction", fbc_lower_bound", "cobra_0_bound"
   "R_ASPCT", "reaction", "gene_product_assoc", "(G_b4244 and G_b4245)"
   "R_BIOMASS_Ec_iJO1366_core_53p95M", "reaction", "reactant", "M_kdo2lipid4_e=0"
   "R_DNAP", "reaction", "kind", "biomass"
   "enz_b4086_b4087_b4088", "enzyme", "composition", "gene=b4087, stoic=1.0; gene=b4086, stoic=2.0; gene=b4088, stoic=50"
   "P0A6D5", "uniprot", "locus", "b1692"
   None, "modelAttrs", "miriamAnnotation", "bqbiol:hasTaxon, taxonomy/2144189; bqmodel:isDescribedBy, pubmed/30657448"


Table ``remove_gps``
^^^^^^^^^^^^^^^^^^^^

The table designated ``remove_gps`` contains a single column listing the gene products to be removed from the model, e.g. dummy proteins.

.. csv-table:: ``remove_gps``: removal of gene products
   :header: "Column", "Contents", "Example"

   "id", "gene product identifier", "G_s0001"


Table ``add_gps``
^^^^^^^^^^^^^^^^^

The table designated ``add_gps`` facilitates adding gene products to the model, e.g. tRNAs to a RBA model.

.. csv-table:: ``add_gps``: extend the model by adding gene products
   :header: "Column", "Contents", "Example"

   "id", "gene product identifier", "G_b0026"
   "label", "gene locus", "b0026"
   "compartment", "compartment identifier (default: None)", "c"

Note: SBML parameters "miriamAnnotation", "xmlAnnotation", "metaid", "sboterm" and "notes" can be configured as well.


Table ``add_species``
^^^^^^^^^^^^^^^^^^^^^

The table designated ``add_species`` facilitated the addition of species to the model, e.g. cofactors or tRNA metabolites. RBA models distinguish between tRNA metabolites used in reactions, e.g. charging and elongation reactions, and tRNA macromolecules, which need to be synthesized and get diluted by growth and have configured concentration targets.

.. csv-table:: ``add_species``: extend model by adding species
   :header: "Column", "Contents", "Example"

   "id", "species identifier", "M_pqq_p"
   "name", "descriptive name ", "pyrroloquinoline quinone(3−)"
   "compartment", "compartment identifier", "p"
   "miriamAnnotation", "MIRIAM annotation string (default: None)", "bqbiol:is, chebi/CHEBI:58442"
   "fbcCharge", "electrical charge (default: None)", -3
   "fbcChemicalFormula", "chemical formula (default: None)", "C14H3N2O8"

Notes: Additional SBML parameters can be configured with this table, which include "xmlAnnotation", "metaid", "notes", "sboterm", "substanceUnits", and bool values for "constant", "boundaryCondition" and "hasOnlySubstanceUnits".

Table ``add_reactions``
^^^^^^^^^^^^^^^^^^^^^^^

The table designated ``add_reactions`` facilitates the addition of reactions to the model, e.g. tRNA charging reactions used in RBA models. The detailed configuration options can be explored when exporting a FBA or XBA model in tabular format, using the converter `sbmlxdf <https://sbmlxdf.readthedocs.io/en/latest/>`_ or the method XbaModel.export('mymodel.xlsx')

.. csv-table:: ``add_reactions``: extend model by adding reactions
   :header: "Column", "Contents", "Example"

   "id", "reaction identifier", "R_GLYTRS"
   "name", "descriptive name", "Glycyl-tRNA synthetase"
   "fbcGeneProdAssoc", "gene product association (default: None)", "assoc=(G_b3559 and G_b3560)"
   "reactionString", "reaction string (default: None)", "M_gly_c + M_trnagly_c + M_atp_c => M_glytrna_c + M_amp_c + M_ppi_c"
   "fbcLb", " lower flux bound (default: None)", 0
   "fbcUb", "upper flux bound (default: None)", 1000

Note: : Additional SBML parameters can be configured with this table. As an alternative to the provision of a "reactionString", the parameters "reactants", "products" and "reversible" can be configured. In lieu of configuring numerical values for "fbcLb" and "fbcUb", model parameter identifiers can be utilized in "fbcLowerFluxBound" and "fbcUpperFluxBound." Additionally, "miriamAnnotation", "xmlAnnotation", "metaid", "sboterm" and "notes" can be configured.


Table ``add_parameters``
^^^^^^^^^^^^^^^^^^^^^^^^

The table designated ``add_parameters`` facilitates the addition of SBML parameters to the model, e.g. specific variable bounds.

.. csv-table:: ``add_parameters``: add SBML parameters
   :header: "Column", "Contents", "Example"

   "id", "parameter identifier", "R_ATPM_upper_bound"
   "name", "descriptive name", "ATPM upper - bound"
   "value", "numerical value", 13.5
   "units", "reference to a SBML units definition (default: 'dimensionless')", "mmol_per_gDW_per_hr"
   "constant", "SBML constant flag (default: True)", None


Table ``chebi2sid``
^^^^^^^^^^^^^^^^^^^

Enzymes utilized in RBA models may encompass cofactors (refer to the "cofactor_flag" parameter in the "general" table). These cofactors are retrieved from UniProt records as ChEBI identifiers and subsequently mapped to the identifiers of the model species. Cofactors that cannot be mapped will not be included in the model. An automatic mapping is facilitated based on the MIRIAM annotation configured on model species. In instances where this mapping proves unsuccessful, ChEBI identifiers can be mapped using the table ``chebi2sid``."

The table designated ``chebi2sid`` maps cofactor ChEBI identifiers to model species.

.. csv-table:: ``chebi2sid``: cofactor mapping to model species
   :header: "Column", "Contents", "Example"

   "chebi", "number part of ChEBI identifier", 61717
   "sid", "species identifier", "M_pheme_c"


ECM configuration
-----------------

The ECM configuration file, used in configuration of the EcModel instance, contains a single table, designated ``general``.

.. csv-table:: ``general``: general configuration data
   :header: "Column", "Contents", "Example"

   "parameter", "parameter name", "ecm_type"
   "value", "value", "GECKO"

Following parameters can be configured with this table:

.. csv-table:: 
   :header: "Parameter", "Description", "Example"

   "ecm_type", "model type to create: 'GECKO', 'ccFBA', 'MOMENT' or 'MOMENTmr' (default: GECKO)", "GECKO"
   "arm_flag", "flag to create 'arm' reactions (default: False)", None
   "avg_enz_sat", "average enzyme saturation value (default: 0.5)", 0.53
   "p_total", "mass fraction of total protein in cellular dry mass (default: 0.5)", 0.57
   "pm2totpm_val_or_paxdb", "mass fraction of modeled protein to total protein, a numerical value or a PaxDB compliant file", "data/511145-WHOLE_ORGANISM-integrated.txt"

'arm' reactions can be created for reactions catalyzed by multiple enzymes to constrain, through adequate flux bounds, the summary flux post-splitting.


RBA configuration
-----------------

The RBA configuration file, used in configuration of the RbaModel instance, contains several tables, with configuration data based on `RbaPy <https://sysbioinra.github.io/RBApy/installation.html>`_.


Table ``general``
^^^^^^^^^^^^^^^^^

The table ``general`` is mandatory and contains a single parameter.

.. csv-table:: ``general``: general configuration data
   :header: "Column", "Contents", "Example"

   "parameter", "parameter name", "avg_enz_sat"
   "value", "value", 0.41

Following parameters can be configured with this table:

.. csv-table:: 
   :header: "Parameter", "Description", "Example"

   "avg_enz_sat", "average enzyme saturation value", 0.41


Table ``trna2locus``
^^^^^^^^^^^^^^^^^^^^

The table designated ``trna2locus`` is mandatory and contains configuration data for tRNA macromolecules.

.. csv-table:: ``trna2locus``: configuration data for tRNA macromolecules
   :header: "Column", "Contents", "Example"

   "rna_id", "arbitrary tRNA identifier", "trnaala"
   "label", "gene locus", "b0203"
   "compartment", "compartment identifier", "c"
   "biomass_aa", "corresponding amino acid metabolite identifier in biomass reaction (default: None)", "M_ala__L_c"

Note: A single gene locus should be selected for a given tRNA type. The mass distribution of tRNA macromolecules can be determined automatically when the corresponding amino acid metabolite identifier is assigned to “biomass_aa”.

.. _function_definitions:

Table ``functions``
^^^^^^^^^^^^^^^^^^^

The table designated ``functions`` contains constant value and function definitions that are used in the RBA model. Either a ``constant`` value or ``function`` definition must be assigned. Once a function has been defined, its name can be used in aggregates.

.. csv-table:: ``functions``: definitions of constant values and RBA functions
   :header: "Column", "Contents", "Example"

   "function_name", "identifier", "frac_protein_c"
   "constant", "numerical values (default: None)", None
   "function", "RBA function definition (default: None)", "LINEAR_CONSTANT=0.7279, LINEAR_COEF=0.04472, X_MIN=0.2, X_MAX=1.9"

Function definitions are based on `RbaPy <https://sysbioinra.github.io/RBApy/installation.html>`_ function definitions and support following types:

.. csv-table:: ``functions definitions``
   :header: "Type", "Parameters", "Example"

   "constant", "'CONSTANT'", "CONSTANT=4.981"
   "linear", "'LINEAR_CONSTANT', 'LINEAR_COEF', 'X_MIN', 'X_MAX', 'Y_MIN', 'Y_MAX'", "LINEAR_CONSTANT=0.7279, LINEAR_COEF=0.04472, X_MIN=0.2, X_MAX=1.9"
   "michaelisMenten", "'kmax', 'Km', 'Y_MIN'", "kmax=86400.0, Km=0.5, Y_MIN=32400.0"
   "exponential", "'RATE'", "variable=growth_rate, RATE=-0.083333"
   "indicator", "'X_MIN', 'X_MAX'", "X_MIN=0.4, X_MAX=0.6"

Notes: the optional parameter ``variable`` can specify a model variable (default: 'growth_rate'). 
'X_MIN' and 'Y_MIN', if not provided, are set to '-inf',  'X_MAX' and 'Y_MAX' to 'inf'.

Table ``compartments``
^^^^^^^^^^^^^^^^^^^^^^

The table designated ``compartments`` is mandatory and contains RBA compartment specific data, including mapping of reaction cids, compartment roles, synthesis targets for dummy proteins (not explicitly modeled proteins) and compartment density constraints.

.. csv-table:: ``compartments``: RBA compartment specific configuration
   :header: "Column", "Contents", "Example"

   "id", "RBA compartment identifier", "om"
   "name", "descriptive name", "outer_membrane"
   "reaction_cids", "reaction-cids mapped", "e-p"
   "keyword", "assigned role: 'cytoplasm', 'medium', 'uptake' or None", "uptake"
   "translation_target_constant", "dummy protein target, numerical value (default: None)", None
   "translation_target_function", "dummy protein target, RBA function (default: None)", None
   "translation_target_aggregate", "dummy protein target, RBA aggregate (default: None)", "aa_conc, inv_avg_protein_len, frac_protein_om, frac_dummy_protein_om"
   "density_contraint_value_type", "density constraint type: 'value', 'upperBound', 'lowerBound' or None", "upperBound"
   "density_constraint_constant", "density, numerical value (default: None)", None
   "density_constraint_function", "density, RBA function (default: None)", None
   "density_constraint_aggregate", "density, RBA aggregate (default: None)", "aa_conc, frac_protein_om"

Notes: Assign value to either “xxx_constant“, “xxx_functionv or “xxx_aggregate“. “xxx_function“ must comply with the function definitions, see :ref:`function_definitions`. “xxx_aggregate“ must reference function names defined in table 'functions'.  

Table ``targets``
^^^^^^^^^^^^^^^^^

The table ``targets`` is mandatory and contains the target concentrations in mmol/gDW for metabolites and macromolecules to implement dilution by growth. It contains reaction flux targets in mmol/gDWh for metabolic reactions, e.g. non-growth associated maintenance, and production/degradation flux targets from macromolecules. The targets can be grouped by 'target_group'.

.. csv-table:: ``targets``: RBA model target concentrations and fluxes
   :header: "Column", "Contents", "Example"

   "target_group", "arbitrary group name",  "mrna_degradation"
   "target_type", "target type: 'concentrations',  'reactionFluxes', 'productionFluxes' or 'degradationFluxes'", "degradationFluxes"
   "target", "metabolite id, macromolecule id, reaction id or a filter", "mrna"
   "target_value_type", "value type: 'value', 'lowerBound', 'upperBound'", "value"
   "target_constant", "numerical value (default: None)", None
   "target_function", "RBA function (default: None)", "LINEAR_CONSTANT=0.2732, LINEAR_COEF=-0.0376, X_MIN=0.2"
   "target_aggregate", "RBA aggregate (default: None)", None

Notes: Assign value to either “xxx_constant“, “xxx_function“ or “xxx_aggregate“. “xxx_function“ must comply with the function definitions, see :ref:`function_definitions`. “xxx_aggregate“ must reference function names defined in table 'functions'. “target” accepts a filter consisting of a reaction identifier and the keyword 'metabolites', 'amino_acids' or 'dna', to automatically generate individual targets. E.g. 'R_BIOMASS_Ec_iML1515_core_75p37M, metabolites' creates individual metabolite targets with “target_constant” scaled by the stoichiometric coefficients in the biomass reaction. With keyword 'amino_acids' the molar distribution of amino acids is used as scaling factor. 

Table ``processes``
^^^^^^^^^^^^^^^^^^^

The table ``processes`` is mandatory and contains the configuration of the process machines that process macromolecules. The composition of the process machine referenced by ``process`` must be configured in the table 'machineries'. The macromolecules referenced in ``set`` and the processing maps used in ``processing_map`` must be defined in the table 'processing_map'.

.. csv-table:: ``processes``: definition of process machines
   :header: "Column", "Contents", "Example"

   "process", "process name", "pm_translation"
   "name", "descriptive name", "protein synthesis"
   "type", "process type: 'production' or 'degradation'", "production"
   "capacity_constant", "process capacity, numerical value (default: None)", None
   "capacity_function", " RBA function (default: None)", None
   "capacity_aggregate", "RBA aggregate (default: None)", "ribosome_efficiency_MM, fraction_active_ribosomes"
   "processing_map", "processing map", "translation"
   "set", "macromolecule set: 'dna', 'rna', 'protein'", "protein"
   "input_filter", "filter for specific macromolecules (default: None)", None

Notes: The "input_filter" can be utilized to select specific macromolecules as inputs for the process machine. The set 'rna', accepts a comma-separated list of regular expression patterns to select RNA macromolecules by identifier matching, e.g. 'trna, rRNA'. 'proteins' accept the keyword 'signal_peptide' to select proteins with a signal peptide. As an alternative option, a list of RBA compartment identifiers separated by commas can be utilized, or a list of gene loci. Assign value to either “xxx_constant“, “xxx_functionv or “xxx_aggregate“. “xxx_function“ must comply with the function definitions, see :ref:`function_definitions`. “xxx_aggregate“ must reference function names defined in table 'functions'. 

Table ``machineries``
^^^^^^^^^^^^^^^^^^^^^

The mandatory table ``machineries`` contains the composition of the process machines.

.. csv-table:: ``machineries``: process machinery composition
   :header: "Column", "Contents", "Example"

   "process", "process name", "pm_translation"
   "id", "component identifier", "b0023"
   "set", "macromolecule set: 'rna', 'protein' or None for metabolites", "protein"
   "label", "gene locus or species identifier", "b0023"
   "stoic", "stoichiometry, negative value for consumption", -1.0
   "compartment", "RBA compartment", "c"
   "gpid", "gene product identifier or None for metabolites", “G_b0023"


Table ``processing_maps``
^^^^^^^^^^^^^^^^^^^^^^^^^

The table ``processing_maps`` is mandatory and contains two distinct configurations. The first configuration pertains to the composition of macromolecule sets, while the second configuration details the processing of these components. 

.. csv-table:: ``processingMap``: macromolecule components and their processing requirements
   :header: "Column", "Contents", "Example"

   "processingMap", "identifier", "translation"
   "set", "macromolecule set: 'dna', 'rna', 'protein' or None", "protein"
   "component", "one-letter identifier in sequence data or keyword", "A"
   "name", "component name or None", "Alanine"
   "weight", "weight with respect to average amino acid weight or None", 1
   "machinery_cost", "processing cost for this component or None", 1
   "reaction_string", "detailed processing reaction", "M_alatrna_c + 2.0 M_gtp_c + 2.0 M_h2o_c => M_trnaala_c + 2.0 M_gdp_c + 2.0 M_pi_c + 3.0 M_h_c"

Notes: In the context of a macromolecule “set“, a specific "component", denoting the one-letter code utilized in the respective sequence data, may be represented multiple times in the table. Values for the parameters “set“, “name“ and “weight“ have to be assigned only once. The parameter "component" also accepts keywords. The utilization of the keyword ‘constantProcessing’ facilitates the definition of an initial setup reaction for a process machine by assigning a value to the parameter "reaction_string." The keyword ‘cofactor’ enables the allocation of values to "weight," "machinery_cost," and "reaction_string," thereby facilitating the specification of cofactor treatment. Finally, the keyword ‘amino_acids’ facilitates the allocation of values to the "machinery_cost" and "reaction_string" categories, a process that is applicable to all amino acid components present in proteins.


TFA configuration
-----------------

The TFA configuration file, used in configuration of the TfaModel instance, contains several tables, which are described in the following.


Table ``general``
^^^^^^^^^^^^^^^^^

The table ``general`` is mandatory as the file name of the TD database must be specified in the ``thermo_data_fname`` parameter. The corresponding file must have the same structure as the file ‘thermo_data.thermodb’ used in the pyTFA package.

.. csv-table:: ``general``: general configuration data
   :header: "Column", "Contents", "Example"

   "parameter", "parameter name", "thermo_data_fname"
   "value", "value", "data/thermo_data.thermodb"

Following parameters can be configured with this table:

.. csv-table:: 
   :header: "Parameter", "Description", "Example"

   "thermo_data_fname", "parameter name", "data/thermo_data.thermodb"
   "mid_regex_pattern", "regular expression pattern to extract metabolite identifier from species identifier", "^M_(\w+)_\w+$"

Table ``td_compartments``
^^^^^^^^^^^^^^^^^^^^^^^^^

The table designated ``td_compartments`` is mandatory and contains the thermodynamics data related to the compartments defined in the model. The number of columns is determined by the number of compartments, with membrane potential columns ``<cid>_mV`` added for each compartment.

.. csv-table:: ``td_compartments``: compartment specific thermodynamics configuration
   :header: "Column", "Contents", "Example"

   "cid", "compartment identifier", "c"
   "ph", "compartmental pH", 7.5
   "ionic_strength_M", "ionic strength in mol/L", 0.25
   "c_min_M", "default minimal metabolite concentration in mol/L", 1.0e-8
   "c_max_M", "default maximal metabolite concentration in mol/L", 0.05
   "<cid>_mV", "membrane potentials in mV (<cid> el. potential - own el. potential)", 150.0


Table ``modify_td_sids``
^^^^^^^^^^^^^^^^^^^^^^^^

The optional table ``modify_td_sids`` facilitates the hard linking of selected metabolites to specific TD data records, thereby overruling the automated matching procedure. Use of this table requires the parameter “mid_regex_pattern” in the table “general” to be configured.

.. csv-table:: ``modify_td_sids``: link metabolite identifiers to TD database records
   :header: "Column", "Contents", "Example"

   "mid", "metabolite identifier, without compartment prefix", "h2o"
   "td_sid", " TD database record identifier", "cpd00001"


Table ``modify_thermo_data``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The table entitled ``modify_thermo_data`` is optional and can be utilized to modify data in the TD database, which has been found to be inconsistent. Attributes employed by f2xba include ‘charge_std’, ‘formula’, ‘deltaGf_std’ and ‘name’.

.. csv-table:: ``modify_thermo_data``: modify metabolite records in TD database
   :header: "Column", "Contents", "Example"

   "td_sid", "TD database record identifier", "cpd00637"
   "attribute", "record attribute", "charge_std"
   "value", "value", 1


Table ``modify_drg0_bounds``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The table entitled ``modify_drg0_bounds`` is optional. Typically, this table is generated automatically as a byproduct of parameter relaxation, and it encompasses adjustments to the lower and upper bounds of the standard transformed Gibbs energy of reaction variables.

.. csv-table:: ``modify_drg0_bounds``: modify lower and upper bounds of variables
   :header: "Column", "Contents", "Example"

   "id", "variable identifier", "V_DRG0_MECDPS"
   "component", "set to 'reaction'", "reaction"
   "attribute", "either 'fbc_lower_bound' or 'fbc_upper_bound'", "fbc_lower_bound"
   "value", "numerical value", 83.3918


Bulk configuration
------------------

This bulk mappings configuration file is used for large scale XbaModel modifications, which includes adding or updating MIRIAM annotations to species, reaction and gene product components. The bulk mappings configuration file can be referenced in the XBA configuration file (table “general”, parameter “bulk_mappings_fname”). The configuration data could be generated based on supplementary information of poorly annotated metabolic models or based on database research.


Table ``species``
^^^^^^^^^^^^^^^^^

The ``species`` table is utilized for the purpose of bulk updating species components within the model. For instance, this facilitates model annotation.

.. csv-table:: ``species``: bulk update of species 
   :header: "Column", "Contents", "Example"

   "id", "species identifier", "M_f6p_c"
   "MA:bigg.metabolite", "add BiGG identifier (default:None)", "f6p"
   "MA:kegg.compound", "add KEGG compound annotation, used for turnover number prediction (default: None)", "C00085"
   "MA:seed.compound", "add SEED compound annotation, used in TD modelling (default: None)", "cpd00072"
   "MA:chebi", "add ChEBI annotation, used for RBA enzyme cofactors (default: None)", "CHEBI:10375"

Note: SBML parameters "compartment", "fbcCharge", "fbcChemicalFormula", "substanceUnits", "constant", "boundaryCondition", "hasOnlySubstanceUnits", "miriamAnnotation", "xmlAnnotation", "metaid", "sboterm" and "notes" can be configured as well.


Table ``reactions``
^^^^^^^^^^^^^^^^^^^

The table entitled ``reactions`` .

The ``reactions`` table is utilized for the purpose of bulk updating reaction components within the model. For instance, this facilitates model annotation.

.. csv-table:: ``reactions``: 
   :header: "Column", "Contents", "Example"

   "id", "variable identifier", "R_PGI"
   "name", "set reaction name (default: None)", "Glucose-6-phosphate isomerase"
   "MA:ec-code", "add EC code annotation(default: None)", "5.3.1.9"

Notes: SBML parameters "reactants", "products", "reversible", "fbcGeneProdAssoc", "fbcLowerFluxBound", "fbcUpperFluxBound", "miriamAnnotation", "xmlAnnotation", "metaid", "sboterm" and "notes" can be configured as well.


Table ``fbcGeneProducts``
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``fbcGeneProducts`` table is employed for the purpose of bulk updating gene product components within the model. Gene locus identifiers provided in the ``label`` column are used to generate new or substitute existing gene products.

.. csv-table:: ``fbcGeneProducts``: bulk modification of gene products
   :header: "Column", "Contents", "Example"

   "id", "gene product identifier", "G_MMSYN1_0445"
   "label", "gene locus identifier (default: None)", "JCVISYN3A_0445"
   "MA:uniprot", "update/create uniprot identifier in MIRIAM annotation (default: None)", "AVX54806.1"

Note: SBML parameters "miriamAnnotation", "xmlAnnotation", "metaid", "sboterm" and "notes" can be configured as well.


Table ``groups``
^^^^^^^^^^^^^^^^

The ``groups`` table contains configuration data that is used to replace the SBML groups component (SBML Level 3 groups package) in the model. In the event that the SBML groups component is not present, it is added to the model. 

.. csv-table:: ``groups``: 
   :header: "Column", "Contents", "Example"

   "id", "group identifier", "g1"
   "name", "descriptive name", "Central metabolism"
   "kind", "nature of the group :'partonomy', 'classification' or 'collection' (default: 'partonomy')", None
   "members", "list of reaction identifiers", "idRef=R_PGI; idRef=R_PFK; idRef=R_FBA; ..."


Turnover numbers
----------------

A template file containing default turnover numbers for each enzyme-catalyzed reaction can be extracted from the XbaModel instance by calling the function “xba_model.export_kcats()”. The turnover number configuration file, which can be referenced in the XBA configuration file (table “general”, parameter “kcats_fname”), contains a single table designated ``kcats`` with the following data:

.. csv-table:: ``kcats``: configuration of turnover numbers
   :header: "Column", "Contents", "Example"

   "key", "reaction identifier", "R_FLVR_iso1"
   "rid", "net reaction identifier (= FBA reaction id)", "R_FLVR"
   "dirxn", "reaction direction (1: forward, -1: reverse)", 1
   "enzyme", "enzyme identifier", "enz_b2763_b2764"
   "kcat_per_s", "turnover number in s-1 (per active site)", 10.27


Enzyme composition
------------------

A template file containing default enzyme compositions can be extracted by calling the function “xba_model.export_enz_composition()”. The enzyme composition file, which can be referenced in the XBA configuration file (table “general”, parameter “enzyme_comp_fname”), contains a single table designated ``enzymes`` with the following data:
 
.. csv-table:: ``enzymes``: configuration of enzyme composition
   :header: "Column", "Contents", "Example"

   "eid", "enzyme identifier", "enz_b2763_b2764"
   "name", "enzyme name", "assimilatory sulfite reductase (NADPH)"
   "composition", "proteins (designated by gene label) and stoichiometry", "gene=b2763, stoic=4.0; gene=b2764, stoic=8.0"
   "active_sites", "number of active sites", 4


Proteomics
----------

Following the optimization of the model across a singular or multiple conditions, the data in the ``proteomics`` table can be utilized to create correlation reports between predicted and measured protein levels, for protein correlation plots, and to add information to results related to proteins. Proteins are designated by their gene ``locus`` in the first column of the table. Columns ``<condition>``, with condition names used during optimization, contain the values in units of mg measured protein to g total protein (mpmf; milli protein mass fraction). 

The ``locus`` field 

.. csv-table:: ``proteomics``: measured protein concentrations in mg/gP
   :header: "Column", "Contents", "Example"

   "locus", "gene locus identifier", "b2763"
   "uniprot", "UniProt identifier", "P17846"
   "description", "protein name", "Sulfite reductase [NADPH] hemoprotein beta-component"
   "gene_name", "gene name", "cysI"
   "mw_Da", "protein molecular weight in g/mol", 63939.93
   "avg_mpmf", "mean value of mpmf across conditions", 0.9764
   "rank", "rank after sorting avg_mpmf", 204   
   "<condition>", "mg measured protein per gram total protein (mpmf)", 1.3147


Thermodynamics Database
-----------------------

The thermodynamics database is required to generate models with thermodynamics constraints. The format of this database corresponds to the TD database used by `pyTFA <https://pytfa.readthedocs.io/en/latest/thermoDB.html>`_ package. The file contents is a pickled Python dictionary compressed by zlib. This Python dictionary contains 4 entries:

.. csv-table:: thermodynamics database
   :header: "Key", "Contents", "Example"

   "name", "database name", "DB_AlbertyUpdate"
   "units", "units ('kcal/mol' or 'kJ/mol')", "kcal/mol"
   "metabolites", "thermodynamics data related to metabolites", "Python dictionary, see below"
   "cues", "thermodynamics data related to cues", "Python dictionary, see below"


Metabolites data
^^^^^^^^^^^^^^^^

Metabolites data is stored under the key ``metabolites``. It is a dictionary with keys set to the metabolite identifier (SEED compound id) and values stored in dictionaries as given below.

.. csv-table:: ``metabolites``: thermodynamics data related to metabolites
   :header: "Key", "Contents", "Example"

   "id", "TD data record identifier", "cpd00002"
   "name", "compound name", "ATP"
   "other_names", "list of alternative names", "['ATP', 'Adenosine 5-triphosphate', 'atp']"
   "nH_std", "number of hydrogen atoms in standard condition", 13
   "mass_std", "molecular mass (g/mol) of compound in standard condition", 504.0
   "formula", "chemical formula in standard condition", "C10H13N5O13P3"
   "deltaGf_std", "Gibbs energy of formation in standard condition", -673.85
   "error", "error flag", "Nil"
   "charge_std", "electrical charge in standard condition", -3
   "deltaGf_err", "estimated error", 3.0431
   "struct_cues", "groups", "{'RWWNW': 1, 'RWCdblWW': 1, 'WNH2': 1, 'HeteroAromatic': 2, ...}"
   "pKa", "list of pKa values", "[0.84, 1.83, 4.68, 7.6, 13.03]"

Cues data
^^^^^^^^^

Cues data is stored under the key ``cues``. It is a dictionary with keys set to the cue identifier containing the data in a dictionary as given below

.. csv-table:: ``cues``: thermodynamics data related to cues
   :header: "Key", "Contents", "Example"

   "id", "cues record identifier", "RWWNW"
   "names", "list of names", ['RWWNW']
   "formula", "electrical formula, if applicable", "N1"
   "charge", "electrical charge", 0
   "energy", "energy contribution", 22.1
   "error", "estimated error", 0.617
   "datfile", "thermodynamics data file, if applicable", "RWWNW.gds"
   "small", "boolean flag", False
