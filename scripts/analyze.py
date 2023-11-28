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

if __name__ == '__main__':
    # Read in data
    fname = "7d040fa6-8f8e-4693-9326-998727b29bc7"
    generation_fitnesses = read_pickle_file("./data/" + fname + ".pickle")
    # schemas = generate_schemas(generation_fitnesses)
    schemata = generate_schemas_from_population(generation_fitnesses["data"][0]["population"])
    counts = count_schemata(schemata, {})
    pass    
    # Write data to file
    with open('./data/schemas.pkl', 'wb') as file:
        pickle.dump(schemas, file)