import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import pickle
# sys.path.append("..")
import csv
from utils import *
import math

def calculate_sample_size(confidence_level, margin_of_error, population_size):
        z = 1.96  # Z-score for 95% confidence level
        p = 0.5  # Assumed proportion (0.5 for maximum sample size)
        q = 1 - p  # Complement of the assumed proportion
        e = margin_of_error

        sample_size = math.ceil((z**2 * p * q) / (e**2))
        adjusted_sample_size = math.ceil(sample_size / (1 + ((sample_size - 1) / population_size)))
        return adjusted_sample_size

def read_pickle_file(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data

def read_csv_file(file_path):
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        data = list(reader)
    return data

if __name__ == "__main__":
    # Read in data
    fname = "01d842b2-f369-4681-b1d3-a29f6f9095a5"
    print("Reading in data from: ", fname)
    df = read_pickle_file("./data/objective/six-hump-camelback/schemata/" + fname + "/03-full-population.pickle")
    params = pd.read_csv("./data/configs/" + fname + ".csv")

    print("# rows: \t", df.shape[0])
    print("# columns: \t", df.shape[1])

    
    print("Calculating average fitness of population each generation...")
    average_fitness_pop = df.groupby('generation')['fitness'].mean()
    average_fitness_pop_df = pd.DataFrame({'generation': average_fitness_pop.index, 'average_fitness': average_fitness_pop.values})

    print("Calculating average fitness of wildcards each generation...")
    average_fitness_wildcard = df.groupby(['generation', 'wildcard', 'defining_length', 'order'])['fitness'].mean().reset_index()
    # print("Unique schemata: \t", len(df["wildcard"].unique()))
    print("Counting wildcard counts per generation...")
    average_fitness_wildcard['wildcard_count'] = df.groupby(['generation', 'wildcard'])['wildcard'].transform('count')

    print("Merging average fitness of wildcards with average fitness of population...")
    average_fitness_wildcard = average_fitness_wildcard.merge(average_fitness_pop_df[['generation', 'average_fitness']], on=['generation'], suffixes=('_wildcard', '_pop'))
    
    print("Calculating Holland's Schema Theorem values...")
    Pc = params["Pc"]
    Pm = params["Pm"]
    l = params["nbits"]

    holland_eq_df = average_fitness_wildcard.copy()
    holland_eq_df["rhs_selection"] = holland_eq_df["wildcard_count"] * abs(holland_eq_df["average_fitness"] - -1.0316) / abs(holland_eq_df["fitness"] - -1.0316 + 0.0001) # small term to avoid division by zero
    holland_eq_df["rhs_crossover"] = holland_eq_df.apply(lambda x: (Pc * x["defining_length"] / (int(l.iloc[0]) - 1)) , axis=1)
    holland_eq_df["rhs_mutation"] = holland_eq_df.apply(lambda x: Pm * x["order"], axis=1)
    holland_eq_df["rhs_value"] = holland_eq_df.apply(lambda x: x["rhs_selection"] * (1 - x["rhs_crossover"] - x["rhs_mutation"]), axis=1)

    # Write the sampled dataframe
    print("Writing data to file...")
    
    output_fname = "./data/objective/six-hump-camelback/schemata/" + fname + "/04-hollands-schema-formulas.pickle"
    print("Filename: ", output_fname)
    dir_name = os.path.dirname(output_fname)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    with open(output_fname, 'wb') as f:
        pickle.dump(holland_eq_df, f)

