Online Resources
================

These classes provide access to online resources. Usually, data will be accessed during the
creation of extended models. Organism specific parameters need to be defined in the
XBA configuration file, sheet 'general'. Data will be accessed from the specified organism directory.
Data will be downloaded from the online resource, only if data cannot be found locally. Data access to
BioCyc online resources may require a subscription.

Data resources can be accessed as well from outside a model generation workflow.

UniProt
-------

Access to protein related information available on UniProt knowledge base.

.. currentmodule:: f2xba

.. autoclass:: f2xba.uniprot.uniprot_data.UniprotData
   :members:
   :member-order: bysource


NCBI
----

Access to genome information available on NCBI data base.

.. currentmodule:: f2xba

.. autoclass:: f2xba.ncbi.ncbi_data.NcbiData
   :members:
   :member-order: bysource


BioCyc
------

Access to enzyme information available on BioCyc data base. A subscription may be required to access
information from certain organism. Access to BioCyc organism specific information should only be configured
in the XBA configuration file, if the data is of high quality for the organism in question.

.. currentmodule:: f2xba

.. autoclass:: f2xba.biocyc.biocyc_data.BiocycData
   :members:
   :member-order: bysource



