Model optimization
------------------

.. currentmodule:: f2xba

.. autoclass:: f2xba.optimization.optimize.Solution

.. autoclass:: f2xba.optimization.optimize.Optimize
   :members: __init__, set_medium, optimize, get_variable_bounds, set_variable_bounds, set_objective, set_tfa_metab_concentrations, get_tx_metab_genes, pfba, fva, gene_deletions

.. autoclass:: FbaOptimization
   :members: __init__, update_relaxation

.. autoclass:: EcmOptimization
   :members: __init__, configure_moment_model_constraints, scale_kcats, unscale_kcats

.. autoclass:: GeckoFitKcats
   :members: __init__, process_data, update_kcats

.. autoclass:: RbaOptimization
   :members: __init__, set_medium_conc, solve








   