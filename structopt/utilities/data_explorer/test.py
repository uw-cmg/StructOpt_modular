
from core import Generations

generations = Generations('genealogy.log')

id10 = generations[-1][10]
print("History:")
id10.print_history(generations)

crossovers, mutations = generations[-1][10].count_modifiers(generations)
print("Crossovers:", crossovers)
print("Mutations:", mutations)

