
.. image:: ./images/ccb_logo-2.1.1_p_2100x2970.png
   :width: 105 px
   :height: 148 px
   :align: right
   :target: https://www.cs.hhu.de/lehrstuehle-und-arbeitsgruppen/computational-cell-biology

Welcome to f2xba's documentation
================================

.. image:: ./images/f2xba.png
   :align: center


f2xba modelling framework: from FBA to extended genome-scale modelling
----------------------------------------------------------------------

In the domain of systems biology, the **f2xba** modeling framework has been developed for the purpose of generating a variety of extended genome-scale metabolic model types using simple and consistent workflows. This modeling framework was developed at the research group for `Computational Cell Biology (CCB) <https://www.cs.hhu.de/en/research-groups/computational-cell-biology>`_  at Heinrich-Heine-University Düsseldorf, Germany.

The CCB research group has developed a suite of `software tools <https://www.cs.hhu.de/en/research-groups/computational-cell-biology/software-contributions>`_ to facilitate genome-scale metabolic modeling. Sybil is an R package that utilizes genome-scale metabolic network optimization through the use of flux balance analysis (FBA)-based methods. SybilccFBA is an extension designed to enhance the optimization of enzyme constraint models. `TurNuP <https://turnup.cs.hhu.de/Kcat>`_ is a machine learning model that predicts turnover numbers, which are required to parametrize extended genome-scale models. `smblxdf <https://sbmlxdf.readthedocs.io/en/latest/>`_  is a Python package that converts between SBML coded genome-scale metabolic models and tabular formats. It is used to create and modify SBML coded models, as well as to access model information.

Extended model types 
--------------------

f2xba support generation of enzyme constraint models, such as GECKO (`Sánchez et al., 2017 <https://doi.org/https://doi.org/10.15252/msb.20167411>`_), ccFBA [1]_, MOMENT (`Adadi et al., 2012 <https://doi.org/10.1371/journal.pcbi.1002575>`_) and MOMENTmr [1]_, resource balance constraint RBA models (`Bulović et al., 2019 <https://doi.org/https://doi.org/10.1016/j.ymben.2019.06.001>`_; `Goelzer et al., 2011 <https://doi.org/https://doi.org/10.1016/j.automatica.2011.02.038>`_), and thermodynamics constraint models, such as TFA (`Henry et al., 2007 <https://doi.org/10.1529/biophysj.106.093138>`_; `Salvy et al., 2019 <https://doi.org/10.1093/bioinformatics/bty499>`_) and TGECKO (thermodynamic GECKO) and TRBA (thermodynamic RBA). These advanced model types, which have been developed in recent years, are based on existing genome-scale metabolic models used for FBA (flux balance analysis), a methodology that has been utilized for decades (`Watson, 1986 <https://doi.org/10.1093/bioinformatics/2.1.23>`_). Genome-scale metabolic models can be obtained from databases such as the BiGG models database (`King, Lu, et al., 2015 <https://doi.org/10.1093/nar/gkv1049>`_), or retrieved from publications.

Relevance of extended modelling
-------------------------------
The advent of high-throughput data has led to a growing importance of these extended models. Fundamentally, FBA can be regarded as a predictor of the macroscopic behavior of metabolic networks, while extended models offer insights into the intricate functioning of these networks. Extended models contain considerably more parameters. While some of these additional parameters require definition, the majority are automatically retrieved from online databases and tools, including NCBI, UniProt, BioCyc, and TurNuP (`Kroll et al., 2023 <https://doi.org/10.1038/s41467-023-39840-4>`_). The development of these extended models and the enhancement of their parameters can be facilitated through simple and consistent workflows. Furthermore, the sharing of configuration data among different model types is encouraged. All extended models are exported in stand-alone SBML (Systems Biology Markup Language) coded files (`Hucka et al., 2019 <https://doi.org/10.1515/jib-2019-0021>`_) to facilitate model sharing and processing by downstream tools, such as cobrapy (`Ebrahim et al., 2013 <https://doi.org/10.1186/1752-0509-7-74>`_). Additionally, the f2xba modeling framework provides optimization support via cobrapy or gurobipy interfaces. Optimization results are structured and enriched with additional data. This includes tables for each variable type, correlation plots, and exports to `Esher <https://escher.github.io>`_ (`King, Dräger, et al., 2015 <https://doi.org/10.1371/journal.pcbi.1004321>`_). This facilitates interpretation of model predictions and supports workflows for model parameter adjustments.

Integrated solution
-------------------

Research groups have already developed tools to support extended genome-scale modeling. These tools have been implemented in various programming environments, each exhibiting a distinct approach to model parametrization, generation, and optimization. However, none of these tools generate stand-alone models coded in SBML. ccFBA and MOMENT modeling is supported by the R package `sybilccFBA <https://cran.r-project.org/src/contrib/Archive/sybilccFBA/>`_, GECKO modeling by the MATLAB package `geckomat <https://github.com/SysBioChalmers/GECKO/tree/main/src>`_, RBA modeling by the Python package `RBApy <https://sysbioinra.github.io/RBApy/installation.html>`_, and thermodynamics modeling by the Python package `pyTFA <https://pytfa.readthedocs.io/en/latest/index.html>`_. f2xba is the first integrated tool to support model generation of various extended model types within a single programming environment, compatible model parametrizations, shareable configuration files, and consistent workflows for both model generation and optimization. The resulting models are exported to files and are fully compliant with the SBML standard. Furthermore, all annotation data from the original genome-scale (FBA) model is carried over. Depending on the availability of organism-specific data and actual requirements, different extended model types and differently parametrized versions of a target organism can be generated with relative ease. It is our hope that the f2xba modeling framework will support the community in actively using these extended model types, which have been published in the previous few years.

Tutorials
---------

The documentation includes a set of tutorials with detailed descriptions, where different types of extended models are created based on the most recent genome-scale metabolic network reconstruction of *Escherichia coli*, iML1515 (`Monk et al., 2017 <https://doi.org/10.1038/nbt.3956>`_). Similar jupyter notebooks are available upon request for the generation of extended models based on yeast9 (`Zhang et al., 2024 <https://doi.org/10.1038/s44320-024-00060-7>`_), *Saccharomyces cerevisiae*, iJN678 (`Nogales et al., 2012 <https://doi.org/10.1073/pnas.1117907109>`_), *Synechocystis* sp. PCC 6803, and MMSYN (`Breuer et al., 2019 <https://doi.org/10.7554/eLife.36842>`_), the synthetic cell JCVI-Syn3A based on *Mycoplasma mycoides capri*.


Outlook
----------
Growth balance analysis (GBA) (`Dourado & Lercher, 2020 <https://doi.org/10.1038/s41467-020-14751-w>`_) modeling is an active research project in CCB. In GBA models, reaction fluxes are coupled with protein requirements using non-linear kinetic functions, where enzyme saturation depends on variable metabolite concentrations. We have previously demonstrated the generation of small, schematic GBA models in SBML, the loading of these models from SBML, and the optimization of them using non-linear solvers. However, the optimization of genome-scale GBA models remains challenging. Once this optimization problem is resolved, f2xba could be extended to support GBA model generation, e.g., by extending GECKO or RBA configuration data, and GBA model optimization, either using nonlinear optimization features available in gurobi 12 or using a dedicated nonlinear solver like IPOP.


References:

.. [1] 
   Desouki, A. A. (2015). sybilccFBA: Cost Constrained FLux Balance Analysis: MetabOlic Modeling with ENzyme kineTics (MOMENT).  
   In CRAN. https://cran.r-project.org/web/packages/sybilccFBA/index.html



