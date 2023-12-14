import itertools
import fnmatch
import numpy as np
import random
import re
import csv
from math import log2
import hashlib
import json
import uuid
import pickle
import os
from utils import *

def generate_id():
    return str(uuid.uuid4())

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

# def encode(solution:tuple, nbits:int, cutpoint:int, bounds:nbits) -> tuple:
#     """Convert (x, y) solution into a bitstring.

#     Args:
#         solution (tuple): _description_
#         bounds (list): dict with a key value for each variable and the bounds (domain).

#     Returns:
#         tuple: (x, y) solution decoded from bitstring.
#     """
#     bv = binary_value(solution)
#     a,b = bounds
#     precision = (b - a) / (pow(2, len(solution)) - 1)
#     return bv * precision + a

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

def crossover(parents:list, strategy="uniform") -> str:
    """
    Randomly select a cutpoint and swap tails 
    between two parents.
    """
    # make sure the strategy is one of the available strategies
    if strategy not in ["single point", "two point", "uniform"]:
        raise Exception
    
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

def obj_shc(solution):
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

def tournament_selection(population:list, k:int) -> str:
    """Select a solution from the population using tournament selection.

    Args:
        population (list): a list of solutions, eg ["010011011",..].
        k (int): the number of solutions to select from the population to compete in the tournament.

    Returns:
        str: the selected solution
    """
    # randomly select k solutions from the population
    # sort the solutions by fitness
    # return the best solution
    competitors = random.sample(population, k)
    winner = min(competitors, key=lambda x: x[2])
    return winner

def GA(
        bounds, # dictionary of domain bounds for each variable x, y
        Pc = 1.00, # crossover probability
        Pm = 0.01, # mutation probability
        n = 24, # population size. the number of solutions per generation.
        nbits = 16, # number of bits in a chromosome / solution
                    # x: first 8 bits, y: last 8 bits
        xycutpoint = 8, # first 4 bits are x, the rest are y
        gmax = 100 # maximum number of generations. 
    ):
    xbounds = bounds["x"]
    ybounds = bounds["y"]

    # BEGIN GENETIC ALGORITHM
    # Initialize a generation 0 population.
    population = [generate_solution(nbits) for _ in range(n)]

    # initialize final_solution and final_fitness variables
    # final_solution = False
    # final_fitness = False
    
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

        # Decode into x, y points for fitness testing
        for solution in population:
            bitx = solution[:xycutpoint]
            bity = solution[xycutpoint:]
            x, y = decode(bitstring=solution, cutpoint=xycutpoint, bounds=bounds)
            # Calculate fitness
            # fitness = objective((x,y))
            fitness = obj_shc((x,y))

            # Store results
            generation_fitnesses.append((solution, (x,y), fitness))

        # Select parents for crossover using tournament selection
        parents = [tournament_selection(generation_fitnesses, 4) for _ in range(n)]
        
        # Create the next generation using uniform crossover and mutation
        new_population = []
        # for i in range(0, len(parents)-1):
        for i in range(n):
            # Take only the chromosome from the parent data structure
            parent1 = random.sample(parents, 1)[0][0]
            parent2 = random.sample(parents, 1)[0][0]
            # parent1 = parents[i][0]
            # parent2 = parents[i + 1][0]

            # Perform uniform crossover
            crossover_point = random.randint(1, 11)
            child = parent1[:crossover_point] + parent2[crossover_point:]

            # Mutate the child
            child = mutate(child, Pm)

            new_population.append(child)
        # Record the best fitness of the generation
        generation_best_solution = min(generation_fitnesses, key=lambda x: x[2])
        generation_best_fitness = generation_best_solution[2]
        
        # Record best fitness for all generations so far
        if g == 0:
            final_solution = generation_best_solution
            final_fitness = generation_best_fitness
        else:
            if generation_best_fitness < final_fitness:
                final_solution = generation_best_solution
                final_fitness = generation_best_fitness
            else:
                pass
        # Save data from this run to a data structure
        state["generation"] = g
        state["population"] = population
        state["generation_fitnesses"] = generation_fitnesses
        state["best_solution"] = generation_best_solution
        state["best_fitness"] = generation_best_fitness
        state["final_solution"] = final_solution
        state["final_fitness"] = final_fitness

        # Append to the records list
        records.append(state)
        
        # Replace the population for the next generation
        population = new_population.copy()
        # End this generation run. 
        
    return records
    
if __name__ == "__main__":
    # wildcard_schemas = generate_schemas(nbits=nbits)
    # Inialize a random seed then change the random seed 10 times.
    # random_seeds = [random.randint(0,1000) for i in range(10)]
    seed = 289908
    random.seed(seed)

    # Six-hump camelback bounds
    xbounds = [-3,3]
    ybounds = [-2,2]

    # Booth's equation bounds
    # xbounds = [-10,10]
    # ybounds = [-10,10]

    bounds = {
            "x": xbounds,
            "y": ybounds
    }

    n = 50   # population size
    Pc = 1.00 # crossover probability
    Pm = 0.01  # mutation probability
    gmax = 50 # number of generations
    nbits = 16 # number of bits in a chromosome / solution
    xycutpoint = 8 # point in the solution to split into x, y

    # test_pop = ["0000000000000000", "0101010101010101", "1100110011001100"]
    # # schemata = generate_alL_schemas(nbits=nbits)
    # # generate_schemas_from_population(population=test_pop, schemata_set=schemata)
    # counts = {}
    # for s in test_pop:
    #     test_schemata = generate_schemas_from_solution(bitstring=s)
    #     counts = count_schemata(schemata=test_schemata, counts_dict=counts)
    # pass
        
    # BEGIN GENETIC ALGORITHM
    simdata = GA(bounds, Pc, Pm, n, nbits, xycutpoint, gmax)
    params = {
        "seed": seed,
        "n": n,
        "Pc": Pc,
        "Pm": Pm,
        "gmax": gmax,
        "nbits": nbits,
        "xycutpoint": xycutpoint,
        "objective": "six-hump camelback"
    }
    results = {
        "params": params,
        "data": simdata
    }
    idx = generate_id()
    output_fname = "./data/objective/six-hump-camelback/schemata/" + idx + "/01-raw.pickle"
    
    print("Writing data for simulation: ", idx)
    print("Filename: ", output_fname)
    
    dir_name = os.path.dirname(output_fname)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    with open(output_fname, 'wb') as f:
        pickle.dump(results, f)
    # # Create pickle file
    # with open(f"./data/{idx}.pickle", "wb") as f:
    #     pickle.dump(results, f)
    
  # Write params file
    headers = results['params'].keys()

    with open(f'./data/configs/{idx}.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerow(results['params'])
    
    
            # Calculate entropy
            # Calculate probability of each occurence
            # probabilities = dict((x, population.count(x)/len(population)) for x in set(population))
            # prob_space = [v for v in probabilities.values()]
            # self_information = -1*sum([p*log2(p) for p in prob_space])
            # print("population entropy: \t", self_information)
           
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