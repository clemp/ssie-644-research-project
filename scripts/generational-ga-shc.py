import itertools
import fnmatch
import numpy as np
import random
import re
import csv
from math import log2



def generate_alL_schemas(nbits: int) -> list:
    """Generates a list of all wildcard schemas of length `nbits`.

    Args:
        nbits (int): number of bits in the binary string

    Returns:
        list: all possible combinations of 0,1,* schemas of length `nbits`
    """
    alphabets = ["0", "1", "*"]
    # generate list of all schemas
    all_schemas = [''.join(i) for i in itertools.product(alphabets, repeat = nbits)]
    # subset of schemas containing wildcards (avoid duplication)
    wildcard_schemas = [s for s in all_schemas if '*' in s]
    # remove the schema containing only wildcards '***..'
    wildcard_schemas = wildcard_schemas[:-1]
    return wildcard_schemas

def generate_schemas_from_solution(bitstring, original=None, index=0) -> list:
    """Given a bitstring solution, this recursive function returns all wildcard allele schemas present in that solution.

    Args:
        solution (str): a bitstring solution, eg "010011011".
    
    Returns:
        list: a list of all wildcard allele schemas present in the solution.
    """
    if original is None:
        original = bitstring
    if index == len(bitstring):
        if bitstring == original:
            return []
        else:
            return [bitstring]
    else:
        combinations = []
        combinations += generate_schemas_from_solution(bitstring, original, index + 1)
        bitstring = bitstring[:index] + '*' + bitstring[index+1:]
        combinations += generate_schemas_from_solution(bitstring, original, index + 1)
        return combinations
    # schemata = []
    # for i, bit in enumerate(solution):
    #     # if this is the first bit, initialize the schemata list
    #     if i == 0:
    #         schemata.append(bit)
    #         schemata.append("*")
    #     else:

def count_schemata(schemata: list, counts_dict: dict) -> dict:
    for element in schemata:
        if element in counts_dict.keys():
            counts_dict[element] += 1
        else:
            counts_dict[element] = 1
    return counts_dict    

def generate_schemas_from_population(population:list, schemata_set:list) -> list:
    """Given a population of `nbit` length bitstring solutions, this returns all wildcard allele schemas present in that population.

    Args:
        population (list): a list of `nbit` length bitstring solutions.

    Returns:
        list: a list of all wildcard allele schemas present in the population.
    """
    if len(set([len(s) for s in population])) != 1:
        raise Exception("Population contains solutions of different lengths. All solutions must be the same length bitstring.")
    
    for solution in population:
        schemata_matches = [s for s in schemata_set if fnmatch.fnmatch(solution, s)]
        print("solution:", solution)
        print("wildcard matches:", schemata_matches)
        # wildcard_schemata = 0
    pass
def order(allele:str) -> int:
    # number of non wildcard bits
    return sum(1 for bit in allele if bit != "*")

def length(schema:str) -> int:
    # find first non-wildcard
    for i, bit in enumerate(schema):
        if bit != "*":
            start = i
            break
    # find last non-wildcard
    for i, bit in enumerate(reversed(schema)):
        if bit != "*":
            end = (len(schema)-1)-i
            break
        
    # start = next(bit for bit in allele if bit != "*")
    # [*,0,0,*,1,*]
    # [^*]
    # end = next(bit for bit in reversed(allele) if bit != "*")
    l = end - start
    return l

def precision(bounds:dict, solution_length:int) -> float:
    """
    {
        "x": [-3,3],
        "y": [-2,2]
    }    
    """
    precisions = {}
    for var in bounds.keys():
        lower = bounds[var][0]
        upper = bounds[var][1]
        
        precision = (upper - lower) / (pow(2, solution_length) - 1)

        precisions[var] = precisions

    return precisions

def encode(solution:tuple, nbits:int, cutpoint:int, bounds:nbits) -> tuple:
    """Convert (x, y) solution into a bitstring.

    Args:
        solution (tuple): _description_
        bounds (list): dict with a key value for each variable and the bounds (domain).

    Returns:
        tuple: (x, y) solution decoded from bitstring.
    """
    bv = binary_value(solution)
    a,b = bounds
    precision = (b - a) / (pow(2, len(solution)) - 1)
    return bv * precision + a

def decode(bitstring:str, cutpoint:int, bounds:list) -> float:
    """
    Encodes a bitstring into an (x, y) solution.
    x and y each have domain boundaries
    
    Args:
        solution (str): a bitstring solution, eg "010011011".
        cutpoint (int): the point in the solution which to split into x, y

    Returns:
        tuple: (x, y)
    """
    # TODO: all kinds of error handling
    xbit = bitstring[:cutpoint]
    ybit = bitstring[cutpoint:]

    xbounds = bounds["x"]
    ybounds = bounds["y"]

    # Int value of binary string.
    # x
    xbv = binary_value(bitstring=xbit)
    xlength = len(xbit)
    xlower = bounds["x"][0] # lower bound
    xupper = bounds["x"][1] # upper bound
    xprecision = (xupper - xlower) / (pow(2, xlength)-1)
    
    # y
    ybv = binary_value(bitstring=ybit)
    ylength = len(ybit)
    ylower = bounds["y"][0] # lower bound
    yupper = bounds["y"][1] # upper bound
    yprecision = (yupper - ylower) / (pow(2, ylength)-1)
    
    xdecoded = xbv * xprecision + xlower
    ydecoded = ybv * yprecision + ylower
    
    return tuple((xdecoded, ydecoded))
    

def generate_solution(length:int) -> str:
    """
    
    """
    # wildcards = re.finditer("[*]", allele)
    # solution = allele
    
    # initialize solution
    solution = ""

    # generate solution from allele, randomly assigning bits for wildcards
    for bit in range(length):
        if random.random() >= 0.5:
            solution += "1"   
        else:
            solution += "0"        
    return solution

def generate_solution_from_schema(schema:str) -> str:
    """
    
    """
    # wildcards = re.finditer("[*]", allele)
    # solution = allele
    
    # initialize solution
    solution = ""

    # generate solution from allele, randomly assigning bits for wildcards
    for bit in schema:
        if bit == "*":
            if random.random() >= 0.5:
                solution += "1"   
            else:
                solution += "0"
        else:
            solution += bit
        
    return solution

def binary_value(bitstring:str) -> int:
    """
    
    """
    value = 0
    for i, bit in enumerate(bitstring):
        if int(bit) == 1:
            value += pow(2, len(bitstring)-1-i)
        else:
            pass
    return value

def crossover(parents:list) -> str:
    """
    Randomly select a cutpoint and swap tails 
    between two parents.
    """
    # make sure there's only two parents
    if len(parents) > 2:
        raise Exception

    p1 = parents[0]
    p2 = parents[1]

    # make sure the parent chromosomes are the same length
    if len(p1) != len(p2):
        raise Exception

    # index to cut the parent chromosomes at
    cutpoint = random.randint(1, len(parents[0])-1)
    if cutpoint == 0:
        pass

    # make the crossover at the cutpoint for both children.
    child1 = parents[0][:cutpoint]
    child1 += parents[1][cutpoint:]

    child2 = parents[1][:cutpoint]
    child2 += parents[0][cutpoint:]

    return [child1, child2]    

def mutate(chromosome:str, pm:float) -> str:
    """
    Bit-wise mutation of prob `pm` on a bitstring.
    """
    solution = ""
    for bit in chromosome:
        if random.random() <= pm:
            if bit == "0":
                solution += "1"
            elif bit == "1":
                solution += "0"
        else:
            solution += bit
    return solution

def objective(solution):
    """Six-hump camelback objective function.

    Args:
        solution (_type_): _description_

    Returns:
        _type_: _description_
    """
    # need to validate solution format here. 
    x,y = solution

    z = (4 - 2.1 * pow(x, 2) + (pow(x, 4) / 3)) * pow(x, 2) + x * y + (-1*4 + 4 * pow(y, 2)) * pow(y, 2)
        
    return z

def obj_booth(solution: tuple) -> float:
    """Booth function.

    Note bounds are [-10,10].

    optimal z = 0.

    Args:
        solution (tuple): (x, y) solution to the function

    Returns:
        float: z value of solution
    """
    x,y = solution
    
    z = pow((x + 2*y - 7), 2) + pow((2*x + y - 5), 2)
    return z
def GA(
        bounds, # dictionary of domain bounds for each variable x, y
        Pc = 1.00, # crossover probability
        Pm = 0.01, # mutation probability
        n = 250, # population size. the number of solutions per generation.
        nbits = 16, # number of bits in a chromosome / solution
                     # x: first 8 bits, y: last 8 bits
        xycutpoint = 8, # first 4 bits are x, the rest are y
        gmax = 100 # maximum number of generations. 
    ):
    xbounds = bounds["x"]
    ybounds = bounds["y"]

    # BEGIN GENETIC ALGORITHM
    # Initialize a generation 0 population.
    population = []

    # fill it with solutions
    for j in range(n): # population size
        solution = generate_solution(nbits)
        population.append(solution)

    # initialize final_solution and final_fitness variables
    final_solution = False
    final_fitness = False
    
    # create data structure to hold records from each generation
    # these will be used to construct the final dataset that's returned
    # by the genetic algorithm.
    records = []
    for g in range(gmax):
        print("Generation: \t", g)
        # Create structure to hold state of this generation
        state = {}
        # Initialize list of fitnesses for this population
        generation_fitnesses = []

        for solution in population:
            # Decode into x, y points for fitness testing
            bitx = solution[:xycutpoint]
            bity = solution[xycutpoint:]
            x, y = decode(bitstring=solution, cutpoint=xycutpoint, bounds=bounds)
            # Calculate fitness
            # fitness = objective((x,y))
            fitness = obj_booth((x,y))

            # Store results
            generation_fitnesses.append((solution, fitness))

        # Rank the initial population (gen 0) to identify the starting optimal solution
        generation_fitnesses = []
        for solution in population:
            # Decode into x,y points for fitness testing
            bitx = solution[:xycutpoint]
            bity = solution[xycutpoint:]
            
            x, y = decode(bitstring=solution, cutpoint=xycutpoint, bounds=bounds)

            # encoded = encode(solution=solution, cutpoint=4, bounds=bounds)
            # fitness = objective((x,y))
            fitness = obj_booth((x,y))
            
            # Add to list that will be sorted later.
            generation_fitnesses.append((solution, (x,y), fitness))

        # TODO: implement tournament selection method.
        # Create the next generation
        generation_fitnesses = sorted(generation_fitnesses, key=lambda x: x[1])
        population = [x[0] for x in generation_fitnesses]
        num_solutions = len(generation_fitnesses)
        
        # rank and create selection probabilities.
        sum_ranks = (num_solutions * (num_solutions + 1) / 2)
        ranks = [round((num_solutions-i)/sum_ranks, 8) for i in range(num_solutions)]

        # normalize the ranks (numpy floating decimal issue)
        ranks = [p/sum(ranks) for p in ranks]
        
        # generation best solution
        generation_solution = generation_fitnesses[0]

        # generation best fitness
        generation_fitness = generation_solution[1]
        
        # if this is the first run, intialize the final solution.
        if final_fitness is False or final_solution is False:
            final_fitness = generation_fitness
            final_solution = generation_solution
        # calculate final fitness
        elif generation_fitness < final_fitness:
            final_fitness = generation_fitness
            final_solution = generation_solution                

        # Create a new population for the next generation.
        new_population = []
        # TODO: create a network model to represent the parent-child relationships
        # generation_fitnesses = []
        for i in range(int(n/2)):
            parents = np.random.choice(population, 2, p=ranks)

            # reproduction steps
            # crossover the parents to create children of the parents with prob Pc
            if random.random() <= Pc:
                children = crossover(parents)
            else:
                children = parents
            # mutate the children with prob Pm (per bit)
            for i, child in enumerate(children):
                mutation = mutate(child, Pc)
                children[i] = mutation

                # add the child the the new population this generation
                new_population.append(children[i])
        # Save data from this run to a data structure
        # TODO: save data to a data structure
        # Replace the previous population with the new population
        # for the next generation.
        print("random seed:", seed)
        print("population id:", idx)
        print("Pc", Pc)
        print("Pm:", Pm)
        print("n:", n)
        print("max G", gmax)
        # print("final solution bitstring", vars["final_solution_bit"])
        print("final solution (x,y)", final_solution)
        # print("final_fitness:", vars["final_fitness"])
        # print("worst fitness", vars["worst_fitness"])
        print("----")

        population = new_population.copy()
        # End this generation run. 
        # Store results, enter next generations.
        final_fitnesses.append(final_fitness)
        fitnesses.append(generation_fitness)
        generations.append(generation_fitnesses)
        simulation_dict = {
            "final_fitness": final_fitness,
            "final_solution_bit": final_solution[0],
            "final_solution_xy": final_solution[1],
            "final_fitnesses": final_fitnesses,
            "fitnesses": fitnesses,
            "generations": generations,
            "worst_fitness": max(fitnesses)
        }
    return simulation_dict
    
if __name__ == "__main__":
    # wildcard_schemas = generate_schemas(nbits=nbits)
    # Inialize a random seed then change the random seed 10 times.
    # random_seeds = [random.randint(0,1000) for i in range(10)]
    seed = 2923408
    random.seed(seed)

    # Six-hump camelback bounds
    # xbounds = [-3,3]
    # ybounds = [-2,2]

    # Booth's equation bounds
    xbounds = [-10,10]
    ybounds = [-10,10]

    bounds = {
            "x": xbounds,
            "y": ybounds
    }

    n = 100   # population size
    Pc = 1.00 # crossover probability
    Pm = 0.01  # mutation probability
    gmax = 100 # number of generations
    nbits = 8 # number of bits in a chromosome / solution
    xycutpoint = 4 # point in the solution to split into x, y

    test_pop = ["0000000000000000", "0101010101010101", "1100110011001100"]
    # schemata = generate_alL_schemas(nbits=nbits)
    # generate_schemas_from_population(population=test_pop, schemata_set=schemata)
    counts = {}
    for s in test_pop:
        test_schemata = generate_schemas_from_solution(bitstring=s)
        counts = count_schemata(schemata=test_schemata, counts_dict=counts)
    pass
        
    # BEGIN GENETIC ALGORITHM
    # Data structure to store different configurations of initial populations
    # initial_populations = []
    # for i in range(6):
    #     seed = random.randint(0,10000)
    #     random.seed(seed)
    population = [] # initialize a generation 0 population.
    # fill it with solutions
    for j in range(n): # population size
        solution = generate_solution(nbits)
        population.append(solution)
    # initial_populations.append(population)
    # Total number of distinct solutions created
    # len(set([y for x in initial_populations for y in x]))

    # Get ready to write some data.
    # with open("SHC-GA-output.csv", 'a') as f:
    #     writer = csv.writer(f, dialect='excel')
    #     headers = ["RANDOM_SEED", "INITIAL_POPULATION_ID", "Pc", "Pm", "POPULATION_SIZE", "MAX_GENERATIONS", "FINAL_SOLUTION_BITS", "FINAL_SOLUTION_XY", "WORST_FITNESS", "FINAL_FITNESS"]
    #     writer.writerow(headers)

        # SET PARAMS
        # n = sim_configs["0"]["n"]
        # Pc = sim_configs["0"]["Pc"]
        # Pm = sim_configs["0"]["Pm"]
        
        # gmax = 100
        # nbits = 8
        # xycutpoint = 4

    # 1. ITERATE THROUGH INITIAL POPULATION CONFIGURATIONS
    
        # Run the GA algo for 5 separate initial populations.
        # for idx in range(len(initial_populations)):
        #     # random seed
        #     population = initial_populations[idx]
            # vars = GA(
            #     init_population=population,
            #     bounds=bounds,
            #     Pc=Pc,
            #     Pm=Pm,
            #     xycutpoint=xycutpoint,
            #     gmax=gmax
            # )
            
            # GENETIC ALGORITHM #
            # xbounds = bounds["x"]
            # ybounds = bounds["y"]
            # xycutpoint = 4
            # gmax = 100

        # initialize final_solution and final_fitness variables
        final_solution = False
        final_fitness = False

        # Iterate through `gmax` number of generations
        for g in range(gmax):
            print("Generation: \t", g)
            # Calculate entropy
            # Calculate probability of each occurence
            # probabilities = dict((x, population.count(x)/len(population)) for x in set(population))
            # prob_space = [v for v in probabilities.values()]
            # self_information = -1*sum([p*log2(p) for p in prob_space])
            # print("population entropy: \t", self_information)
            # Initialize list of fitnesses for this population
            generation_fitnesses = []
            # Initialize schema counts as 0 for all schemas for this generation.
            # schema_counts = { k: 0 for k in wildcard_schemas }
            # schema_fitnesses = { k: [] for k in wildcard_schemas }

            # # Quantify num. schema H at step `g`
            for solution in population:
            #     count = [s for s in wildcard_schemas if fnmatch.fnmatch(solution, s)]
            #     for c in count:
            #         schema_counts[c] += 1
                # Decode into x, y points for fitness testing
                bitx = solution[:xycutpoint]
                bity = solution[xycutpoint:]
                x, y = decode(bitstring=solution, cutpoint=xycutpoint, bounds=bounds)
                # Calculate fitness
                # fitness = objective((x,y))
                fitness = obj_booth((x,y))
        
                # Store results
                generation_fitnesses.append((solution, fitness))
            # # Average fitness of this generation
            # generation_avg_fitness = np.average([s[1] for s in generation_fitnesses])
            # # Average fitness of solutions containing each schema
            # for schema in schema_counts.keys():
            #     schema_population = [s for s in generation_fitnesses if fnmatch.fnmatch(s[0], schema)]
            #     # Average fitness
            #     schema_fitnesses[schema] = np.average([s[1] for s in schema_population])
            #     # for debugging
            #     if schema_counts[schema] > 15:
            #         pass
            #     # Schema theorem formula
            #     # schema_counts[schema] * (1 - Pc * (length(schema)/nbits-1) - Pm * order(schema))
            #     Mht =  schema_counts[schema] * (1 - Pc * (length(schema)/nbits-1) - Pm * order(schema))
            #     # if schema == '00*00***':
            #     #     print(schema, ": \t", Mht)
            #     pass
            # TODO: implement tournament selection method.
            # Create the next generation
            generation_fitnesses = sorted(generation_fitnesses, key=lambda x: x[1])
            population = [x[0] for x in generation_fitnesses]
            num_solutions = len(generation_fitnesses)
            
            # rank and create selection probabilities.
            sum_ranks = (num_solutions * (num_solutions + 1) / 2)
            ranks = [round((num_solutions-i)/sum_ranks, 8) for i in range(num_solutions)]

            # normalize the ranks (numpy floating decimal issue)
            ranks = [p/sum(ranks) for p in ranks]
            
            # generation best solution
            generation_solution = generation_fitnesses[0]

            # generation best fitness
            generation_fitness = generation_solution[1]
            
            # if this is the first run, intialize the final solution.
            if final_fitness is False or final_solution is False:
                final_fitness = generation_fitness
                final_solution = generation_solution
            # calculate final fitness
            elif generation_fitness < final_fitness:
                final_fitness = generation_fitness
                final_solution = generation_solution                

            # Create a new population for the next generation.
            new_population = []
            # TODO: create a network model to represent the parent-child relationships
            # generation_fitnesses = []
            for i in range(int(n/2)):
                parents = np.random.choice(population, 2, p=ranks)

                # reproduction steps
                # crossover the parents to create children of the parents with prob Pc
                if random.random() <= Pc:
                    children = crossover(parents)
                else:
                    children = parents
                # mutate the children with prob Pm (per bit)
                for i, child in enumerate(children):
                    mutation = mutate(child, Pc)
                    children[i] = mutation

                    # add the child the the new population this generation
                    new_population.append(children[i])
            # Save data from this run to a data structure
            # TODO: save data to a data structure
            # Replace the previous population with the new population
            # for the next generation.
            print("random seed:", seed)
            print("population id:", idx)
            print("Pc", Pc)
            print("Pm:", Pm)
            print("n:", n)
            print("max G", gmax)
            # print("final solution bitstring", vars["final_solution_bit"])
            print("final solution (x,y)", final_solution)
            # print("final_fitness:", vars["final_fitness"])
            # print("worst fitness", vars["worst_fitness"])
            print("----")

            population = new_population.copy()
        
            

                # row = [
                #     seed,
                #     idx,
                #     Pc,
                #     Pm,
                #     n,
                #     gmax,
                #     vars["final_solution_bit"],
                #     vars["final_solution_xy"],
                #     vars["worst_fitness"],
                #     vars["final_fitness"]
                # ]
                # writer.writerow(row)

            # 2. EXERIMENT WITH CROSSOVER, MUTATION PROBS AND POPULATION SIZE
            # use the last `initial solution` from the list above.
            # if idx == len(initial_populations)-1:
            #     for config_id in range(1, len(sim_configs)):
            #         # Create the population for this configuration
            #         population = [] # initialize a population.
            #         # fill it with solutions
            #         n = sim_configs[str(config_id)]["n"]
            #         for j in range(n): # population size
            #             solution = generate_solution(nbits)
            #             population.append(solution)
            #         Pc = sim_configs[str(config_id)]["Pc"]
            #         Pm = sim_configs[str(config_id)]["Pm"]

            #         vars = GA(
            #             init_population=population,
            #             bounds=bounds,
            #             Pc=Pc,
            #             Pm=Pm,
            #             xycutpoint=xycutpoint,
            #             gmax=gmax
            #         )

            #         print("random seed:", seed)
            #         print("population id:", idx)
            #         print("Pc", Pc)
            #         print("Pm:", Pm)
            #         print("n:", n)
            #         print("max G", gmax)
            #         print("final solution bitstring", vars["final_solution_bit"])
            #         print("final solution (x,y)", vars["final_solution_xy"])
            #         print("final_fitness:", vars["final_fitness"])
            #         print("worst fitness", vars["worst_fitness"])
            #         print("----")

            #         row = [
            #             seed,
            #             idx,
            #             Pc,
            #             Pm,
            #             n,
            #             gmax,
            #             vars["final_solution_bit"],
            #             vars["final_solution_xy"],
            #             vars["worst_fitness"],
            #             vars["final_fitness"]
            #         ]
                    # writer.writerow(row)