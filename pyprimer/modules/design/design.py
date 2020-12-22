from deap import creator, base, tools, algorithms
from pyprimer.modules import PPC
import numpy as np
from pyprimer.utils.sequence import PCRPrimer, Sequence, READ_MODES
from tqdm import trange

class PrimerDesigner:
    def __init__(self, ppc: PPC):
        self.ppc = PPC
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)

        self.toolbox = base.Toolbox()

        self.toolbox.register("nuclei", np.random.choice, list("ATGC"))
        self.toolbox.register("individual", tools.initRepeat, creator.Individual, self.toolbox.nuclei, n=40)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

        def eval_population(population):
            pop = np.asarray(population)
            dps = test_pcr.get_primer_metrics(pop[:, :20], pop[:, 20:], ["", ""], nCores=8)
            return dps,

        self.toolbox.register("evaluate", eval_population)
        self.toolbox.register("mate", tools.cxTwoPoint)
        self.toolbox.register("mutate", tools.mutShuffleIndexes, indpb=0.05)
        self.toolbox.register("select", tools.selTournament, tournsize=3)

        self.population = self.toolbox.population(n=200)
        self.NGEN=100

    def design(self):
        with trange(self.NGEN) as t:
            for gen in t:
                offspring = algorithms.varAnd(self.population, self.toolbox, cxpb=0.2, mutpb=0.2)
                fits = self.toolbox.evaluate(offspring)
                for fit, ind in zip(fits, offspring):
                    ind.fitness.values = fit

                self.population = self.toolbox.select(offspring, k=len(self.population))
                t.set_postfix(max_fitness=np.max(fits), min_fitness=np.min(fits))

        top = tools.selBest(self.population, k=3)
        print(["".join(x) for x in top])

if __name__=="__main__":
    data_dir = "/bi/aim/scratch/afrolova/COVID19/github/pyprimer/data"

    test_sequence = Sequence(READ_MODES.FASTA)
    sequences_df = test_sequence.describe_sequences(f"{data_dir}/merged.fasta")

    test_pcr = PPC(None, sequences_df.sample(100))
    designer = PrimerDesigner(test_pcr)

    try:
        designer.design()
    except KeyboardInterrupt:
        print(["".join(x) for x in designer.population])