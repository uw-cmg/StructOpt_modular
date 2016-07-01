Core Concepts
#############

Cost Function
=============

Individual
==========

Structure Types
---------------

Crystal
"""""""
Implemented by default.

Has periodic boundary conditions along all dimensions. The entire model is relaxed.

Cluster
"""""""
Not implemented.

Does not have periodic boundary conditions. The entire model is relaxed.

Defect
""""""
Not implemented.

There are three layers to this structure type. The outer layer contains the fixed atoms that are never moved but are used as constrains in the relaxations and fitnesses. The middle layer contains atoms that are part of the crystal structure but that will be mutated, crossed-over, and relaxed. The inner layer contains the defect and it will also be mutated cross-over, and relaxed.

.. figure:: static/images/defect_regions.png
   :align: center

Surface
"""""""
Not implemented.


Population
==========

Crossovers
==========

Selection Schemes
=================

Mutations
=========

Predators
=========

Fingerprinters
==============

Relaxations
===========

Fitnesses
=========

