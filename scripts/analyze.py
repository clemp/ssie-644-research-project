import pandas as pd
import numpy as np
import pickle
import csv
from utils import *

def read_pickle_file(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data

def read_csv_file(file_path):
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        data = list(reader)
    return data

def generate_schemas(generation_fitnesses):
    schemas = []
    for generation in generation_fitnesses:
        for solution in generation:
            # Generate schema from solution
            schema = generate_schema(solution)
            schemas.append(schema)
    return schemas

def generate_schema(solution):
    # Generate schema from solution
    pass

def stats_generation_schema_fitnesses(generation_schemata_fitnesses):
    """Return stats on the schema fitnesses for each generation.

    Min fitness, max fitness, avg fitness, median fitness, std dev fitness
    by schema group (e.g. 000*, 0*0*, 0**0, etc.)

    Args:
        generation_schemata_fitnesses (_type_): _description_
    """
    if not generation_schemata_fitnesses:
        return {}

    simulation_stats = {}
# TODO: implement this for each generation to create the list of average fitnesses for each schema group
    for generation, state in enumerate(generation_schemata_fitnesses.values()):
        generation_min_fitnesses = []
        generation_max_fitnesses = []
        generation_median_fitnesses = []
        generation_stats = {}
        generation_schemata = generation_schemata_fitnesses[generation].keys()
    # list of average fitnesses across all schema groups in this generation
        fitnesses = [state[schema]["avg_fitness"] for schema in generation_schemata]
        # TODO: implement this for each generation to create the list of average fitnesses for each schema group
        # This is input to exploratory stats like min, max, median, std dev
        # [generation_schemata_fitnesses[99][schema]["avg_fitness"] for schema in [generation_schemata_fitnesses[g] for g in generation_schemata_fitnesses][99]]
        generation_min_fitnesses.append(np.min(fitnesses))
        generation_max_fitnesses.append(np.max(fitnesses))
        generation_median_fitnesses.append(np.median(fitnesses))

        generation_stats =  {
            'min_fitness': generation_min_fitnesses,
            'max_fitness': generation_max_fitnesses,
            # 'avg_fitness': avg_fitnesses,
            'median_fitness': generation_median_fitnesses,
            # 'std_dev_fitness': std_dev_fitnesses
        }
        simulation_stats[generation] = generation_stats
    return simulation_stats

if __name__ == '__main__':
    # Read in data
    fname = "82375c64-86b4-43e7-a5d9-b98be57e8d41"
    generation_fitnesses = read_pickle_file("./data/" + fname + ".pickle")
    # schemas = generate_schemas(generation_fitnesses)
    schemata_fitnesses = {}
    generation_schemata_fitnesses = {}
    # each schema is a key, with a second dictionary with keys: "avg_fitness", "count"
    # { "000*": {"avg_fitness": 0.2, "count": 5}, "0*0*": {"avg_fitness": 0.4, "count": 2}, "0**0": {"avg_fitness": 0.7, "count": 5} }
    # generate_population_from_schemata("000*0*10000*0*1")

    for generation in generation_fitnesses["data"]:
        schemata = generate_schemas_from_population(generation["population"])
        counts = count_schemata(schemata, {})
        # print(counts)
        # Initialize a dictionary to store the average fitness of each schema this generation
        schemata_fitnesses = {}
        # Calculate average schema fitness by solution
        ls = [length(s) for s in schemata]
        for solution in generation["generation_fitnesses"]:
            solution_schemata = generate_schemas_from_solution(solution[0]) 
            for schema in solution_schemata:
                if schema in schemata_fitnesses.keys():
                    schemata_fitnesses[schema]["sum_fitness"] += solution[2]
                    schemata_fitnesses[schema]["count"] += 1
                    schemata_fitnesses[schema]["avg_fitness"] = schemata_fitnesses[schema]["sum_fitness"] / schemata_fitnesses[schema]["count"]
                else: # initialize this new schema
                    schemata_fitnesses[schema] = {"sum_fitness": solution[2], "avg_fitness": solution[2], "count": 1}
        generation_schemata_fitnesses[generation["generation"]] = schemata_fitnesses
        # print(generation_schemata_fitnesses)
    stats = stats_generation_schema_fitnesses(generation_schemata_fitnesses)
    print(stats)
    # Write data to file
    # with open('./data/schemas.pkl', 'wb') as file:
    #     pickle.dump(schemas, file)