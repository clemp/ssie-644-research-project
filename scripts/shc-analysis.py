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
    df = read_pickle_file("./data/objective/six-hump-camelback/schemata/" + fname + "/02-transformed.pickle")
    

    print("# rows: \t", df.shape[0])
    print("# columns: \t", df.shape[1])

    df.head()

    # print("Unique schemata: \t", len(df["wildcard"].unique()))
    
    # Create a dataframe of the unique schemata
    # This is to more efficiently calculate defining length and order of each 
    print("Calculating defining length and order of each schema...")
    unique_wildcards = df["wildcard"].unique()
    df_unique_wildcards = pd.DataFrame(unique_wildcards, columns=["wildcard"])

    # Calculate defining length and order of each schema
    df_unique_wildcards["defining_length"] = df_unique_wildcards["wildcard"].apply(lambda x: defining_length(x))
    df_unique_wildcards["order"] = df_unique_wildcards["wildcard"].apply(lambda x: order(x))

    print("Joining back to original dataframe...")
    # Now join back to the original dataframe
    df = df.merge(df_unique_wildcards, on="wildcard")

    print("Calculating record counts by generation...")
    # Calculate the count of records per generation
    df_count_per_generation = df.groupby("generation").size().reset_index(name="count")
    
    print("Calculating sample size per generation...")
    # Add a column to df_count_per_generation dataframe that is the output of "calculate_sample_size"
    df_count_per_generation["sample_size"] = df_count_per_generation["count"].apply(lambda x: calculate_sample_size(confidence_level=0.95, margin_of_error=0.05, population_size=x))
    
    # Create a new dataframe to store the sampled data
    df_sampled = pd.DataFrame()

    # Iterate over each generation
    for generation in df["generation"].unique():
        # Sample the specified number of rows for the current generation
        df_generation_sampled = df[df["generation"] == generation].sample(n=df_count_per_generation.loc[df_count_per_generation["generation"] == generation, "sample_size"].values[0], replace=False)
        
        # Concatenate the sampled data to the new dataframe
        df_sampled = pd.concat([df_sampled, df_generation_sampled])

    # Reset the index of the new dataframe
    df_sampled.reset_index(drop=True, inplace=True)

    # Write the sampled dataframe
    print("Writing data to file...")
    
    output_fname = "./data/objective/six-hump-camelback/schemata/" + fname + "/03-sampled.pickle"
    print("Filename: ", output_fname)
    dir_name = os.path.dirname(output_fname)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    with open(output_fname, 'wb') as f:
        pickle.dump(df_sampled, f)
    print(df_sampled)

