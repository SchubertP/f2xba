
Support functions
=================

.. currentmodule:: f2xba

Managing configuration files
----------------------------

Support functions to load and write configuration files used in extended model generation.

.. autofunction:: f2xba.utils.mapping_utils.load_parameter_file

.. autofunction:: f2xba.utils.mapping_utils.write_parameter_file


Calculate molecular weights
---------------------------

Calculate molecular weights for DNA, RNA, proteins and metabolites.

.. autofunction:: f2xba.utils.calc_mw.calc_mw_from_formula

.. autofunction:: f2xba.utils.calc_mw.protein_mw_from_aa_comp

.. autofunction:: f2xba.utils.calc_mw.rna_mw_from_nt_comp

.. autofunction:: f2xba.utils.calc_mw.dsdna_mw_from_dnt_comp


SGKO
----

Support functions for single gene knockout analysis.

.. autofunction:: f2xba.utils.sgko_utils.confusion_matrix

.. autofunction:: f2xba.utils.sgko_utils.export_gene_predictions


