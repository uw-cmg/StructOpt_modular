Algorithm Workflow
##################

General Overview
================

StructOpt uses a Genetic Algorithm to optimize a set of atomic structures according to a customizable objective function (aka cost function).


Genetic Algorithm
-----------------

A genetic algorithm utilizes a population of structures rather than a single individual. A genetic algorithm, or evolutionary algorithm, is conceptually similar to genetic Darwinism where animals are replaced by "individuals" (in the case of StructOpt a "individual" is an atomic model). A population of atomic models is first generated. Given a population, pairs of individuals are mated (aka crossed over) by selecting different aspects of each model and pasting them into each other. Crossovers always produce two children, one for each section of the models combined together. The offspring are added to the population. After the mating scheme has finished, single individuals can "mutate" (i.e. moving atoms in a unique way) to add new genes to the population's gene pool. After the atoms have been moved via crossovers and mutations, the structures are relaxed. Finally, each structure is run though a series of "fitness" evaluations to determine how "fit to survive" it is, and the population is then reduced to its original size based on a number of optional selection criteria. This process is repeated many times.

In summary:

0. Generate initial structures
1. Locally relax structures
2. Calculate fitness values (e.g. energies) of each structure
3. Remove some individuals from the population based on their fitness value
4. Perform crossovers and selected individuals to generate offspring for the next generation
5. Perform mutations on the selected individuals in the current population and offspring for the next generation
6. Repeat steps 1-5 until the convergence criteria are met

Relevant References
+++++++++++++++++++

* Crystals
  * Artem Oganov
  * Alex Zunger
  * Scott Woodley
  * Richard Catlow
* Clusters
  * Roy L. Johnston
  * Bernd Hartke
  * David Deaven
* Surfaces
  * Cristian V. Ciobanu
  * Kai-Ming Ho



Interaction Between Crossovers/Mutations and Relaxation/Fitness Evalutions
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Crossovers and mutations modify individuals. During relaxation and fitness evaluations, only individuals that have been modified by mutations and crossovers are computed. This avoids recomputing relaxations and fitnesses for individuals who were unchanged during the generation's crossover/mutation/selection scheme.

During crossovers, the offspring are collected into a list. After all crossovers have been completed, these offspring are added to the entire population. Each individual in the entire population then has a chance to be mutated. There will theremore be a variable number of modified individuals that will need to be relaxed and fitnessed at each generation. The number of modified individuals can only be predicted by using the probability of mutation and crossover; however, the number of modified individuals will likely never be more than twice the size of the original population at the start of the generation (although this depends on how the crossover scheme selects individuals to be crossed).

Open questions:
+++++++++++++++

Should every individual be able to be crossed with every other individual?

Or should crossovers only be done pair-wise as a function of their fitness value?

I think we should have different crossover schemes to select which individuals are crossed.

We also need to keep in mind genetic diversity. If duplicate individuals are allowed in the population, this limits genetic diversity because there are fewer spots available for unique individuals.

What does it mean to have the chance of crossover be 100%? Is that all individuals mate with all other individuals, or at most len(population) children can be generated?

Should we put a cap on the number of children that can be generated?


Crossovers and Selection Schemes
++++++++++++++++++++++++++++++++

Predators and Selection Schemes
+++++++++++++++++++++++++++++++

