Algorithm Workflow
##################

Interaction Between Crossovers/Mutations and Relaxation/Fitness Evalutions
==========================================================================

Crossovers and mutations modify individuals. During relaxation and fitness evaluations, only individuals that have been modified by mutations and crossovers are computed. This avoids recomputing relaxations and fitnesses for individuals who were unchanged during the generation's crossover/mutation/selection scheme.

During crossovers, the offspring are collected into a list. After all crossovers have been completed, these offspring are added to the entire population. Each individual in the entire population then has a chance to be mutated. There will theremore be a variable number of modified individuals that will need to be relaxed and fitnessed at each generation. The number of modified individuals can only be predicted by using the probability of mutation and crossover; however, the number of modified individuals will likely never be more than twice the size of the original population at the start of the generation (although this depends on how the crossover scheme selects individuals to be crossed).

Open questions:
---------------

Should every individual be able to be crossed with every other individual?

Or should crossovers only be done pair-wise as a function of their fitness value?

I think we should have different crossover schemes to select which individuals are crossed.

We also need to keep in mind genetic diversity. If duplicate individuals are allowed in the population, this limits genetic diversity because there are fewer spots available for unique individuals.

What does it mean to have the chance of crossover be 100%? Is that all individuals mate with all other individuals, or at most len(population) children can be generated?

Should we put a cap on the number of children that can be generated?


Crossovers and Selection Schemes
================================

Predators and Selection Schemes
===============================

