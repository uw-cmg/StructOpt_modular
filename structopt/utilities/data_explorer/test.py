import sys

from core import DataExplorer


dir = sys.argv[1]

de = DataExplorer(dir)
print(de.parameters)

id_ = 149
indiv = de.generations[-1][id_]
print("History:")
indiv.print_history(de.generations)

crossovers, mutations = de.generations[-1][id_].count_modifiers(de.generations)
print("Crossovers:", crossovers)
print("Mutations:", mutations)

indiv.load_structure()
print(indiv.get_positions()[:5])
print(indiv.relaxations)
