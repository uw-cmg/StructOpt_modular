StructOpt User Documentation
############################

StructOpt is a structure optimization framework that incorportates multiple forward simulation techniques into its minimization scheme. It is designed with modularity in mind, and forward simulation techniques that accept an atomic model as input and output a fitness value can be integrated into this framework. [We need much more introduction about reverse structure optimization in general here.]

This documentation serves as both a user and developer guide for StructOpt. StructOpt is a general optimizer for Python targeted at identifying stable atomic structures. StructOpt was designed to provide both flexibility and simplicity.

While there are many optional parameters to allow the user great customizability, there are also many default settings that allow for faster more simple runs.

The `Examples <>`_ section provides multiple examples of basic StructOpt configurations. Details on input and output structures can be found in the `Input and Output Syntax <>`_ section. Details on the many options currently available for StructOpt are provided in the `Functionalities <>`_ section of this document.  An explanation of commonly generated errors and troubleshooting advice is provided in the section entitled `Troubleshooting <>`_.


Structure Types
===============

Crystal
-------
Not implemented.

Has periodic boundary conditions along all dimensions. The entire model is relaxed.

Cluster
-------
Not implemented.

Does not have periodic boundary conditions. The entire model is relaxed.

Defect
------
Not implemented.

There are three layers to this structure type. The outer layer contains the fixed atoms that are never moved but are used as constrains in the relaxations and fitnesses. The middle layer contains atoms that are part of the crystal structure but that will be mutated, crossed-over, and relaxed. The inner layer contains the defect and it will also be mutated cross-over, and relaxed.

.. figure:: static/images/defect_regions.png
   :align: center

Surface
-------
Not implemented.

