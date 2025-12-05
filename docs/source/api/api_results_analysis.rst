Results analysis
================

FBA model optimization determines values for only one type of optimization
variables, namely reaction fluxes. Optimization of extended models determines
values for different types of optimization variables. Variables include reaction
fluxes, which may be split per catalyzing isoenzyme and direction, but also include
variables for protein, enzyme and metabolite concentration. Meaningful presentation of
optimization results become more challenging. Further, predicted values may be
compared to measured values, results may be enriched with information,
and data may be converted for downstream processing.

.. toctree::

    api_fba_results
    api_ecm_results
    api_rba_results


Results
-------

Base class for optimization results analysis.
Analysis is independent of the interface used to retrieve the optimization solutions.

.. currentmodule:: f2xba

.. autoclass:: f2xba.optimization.results.Results
   :members:
   :member-order: bysource

