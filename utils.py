"""
Utility functions that main.py uses. Here we define the individual processes for sequence reading, SBT insertion, SBT
querying, and SBT saving. In addition, the processes also report metrics such as time elapsed, false positive rate, and
memory used. Different hash and similarity functions are also included here for use.
"""
import math
import time
import os
import pandas as pd
from SBT.SBT import SBT
import random
from collections import defaultdict


# Print parameters being used
def print_params(p):
    for key, value in p.items():
        print(key.ljust(20), value)
    print()


# Read sequences from data file
# @profile
def read_sequences(file_names, sequence_len, dictionary):
    sequences = {}
    start = time.time()
    for file_name in file_names:
        # Read test sequences and insert
        f = open(file_name, 'r')
        sequences[file_name] = f.readline()[:sequence_len]
        f.close()
    end = time.time()
    dictionary["read_time"] = end - start
    print("Read Time           ", dictionary["read_time"])
    return sequences


# Insert sequences into SBT using potentially different clustering methods
# @profile
def insert_sequences(sbt, sequences, bits_to_check, dictionary, method="Greedy"):
    start = time.time()
    if method == "Cluster1":
        sbt.insert_cluster_sequences1(sequences=sequences.values(), experiment_names=sequences.keys(),
                                      bits_to_check=bits_to_check)
    elif method == "Cluster2":
        sbt.insert_cluster_sequences2(sequences=sequences.values(), experiment_names=sequences.keys(),
                                      bits_to_check=bits_to_check)
    else:
        for name, sequence in sequences.items():
            sbt.insert_sequence(sequence=sequence, experiment_name=name)
    end = time.time()
    dictionary["insert_time"] = end - start
    print("Insert Time         ", dictionary["insert_time"])


# Query from SBT and report results
# mode in ("Normal", "Fast", "Faster")
# repeat: number of times to run queries
# @profile
def query_sequences(sbt, all_sequences, dictionary, num_queries, query_size, method="Normal", boyer_moore="False"):
    queries = []
    hits = defaultdict(list)
    for name, sequence in all_sequences.items():
        idx = random.randint(0, len(sequence) - query_size)
        query = sequence[idx:idx+query_size]
        queries += [query]
        if not boyer_moore:
            hits[query] += [name]
    # Report Boyer-Moore time to get an idea of how fast SBT runs and to verify hits
    if boyer_moore:
        start = time.time()
        for query in queries:
            for name, sequence in all_sequences.items():
                if query in sequence:
                    hits[query] += [name]
        end = time.time()
        dictionary["boyer_moore_time"] = end - start
        print("Boyer-Moore Time    ", dictionary["boyer_moore_time"])

    # Begin querying sequences
    start = time.time()
    total_positives = 0
    true_positives = 0
    false_negatives = 0
    total_negatives = 0
    queries_done = 0
    while queries_done < num_queries:
        for query in queries:
            queries_done += 1
            if method == "Fast":
                results = sbt.fast_query_sequence(sequence=query)
            else:
                results = sbt.query_sequence(sequence=query)
            total_positives += len(results)
            total_negatives += len(all_sequences) - len(results)
            for name in hits[query]:
                true_positives += name in results
                false_negatives += name not in results
    end = time.time()
    dictionary["query_time"] = end - start
    dictionary["true_positives"] = true_positives
    dictionary["true_negatives"] = total_negatives - false_negatives
    dictionary["false_positives"] = total_positives - true_positives
    dictionary["false_negatives"] = false_negatives
    dictionary["false_positive_rate"] = 0
    if dictionary["false_positives"] + dictionary["true_negatives"] > 0:
        dictionary["false_positive_rate"] = dictionary["false_positives"] / (dictionary["false_positives"] +
                                                                             dictionary["true_negatives"])
    print("Query Time          ", dictionary["query_time"])
    print("True Positives      ", dictionary["true_positives"])
    print("True Negatives      ", dictionary["true_negatives"])
    print("False Positives     ", dictionary["false_positives"])
    print("False Negatives     ", dictionary["false_negatives"])
    print("False Positive Rate ", dictionary["false_positive_rate"])


# Load SBT and report size
# @profile
def load_sbt(file_name):
    start = time.time()
    sbt = SBT.load(file_name)
    end = time.time()
    print("Load Time           ", end - start)
    return sbt


# Save SBT and report size
# @profile
def save_sbt(sbt, file_name, dictionary):
    start = time.time()
    sbt.save(file_name)
    end = time.time()
    dictionary["save_time"] = end - start
    dictionary["sbt_size"] = os.stat(file_name)[6]
    print("Save Time           ", dictionary["save_time"])
    print("SBT Size (Bytes)    ", dictionary["sbt_size"])


# Print Graph Itself
def print_graph(sbt, print_sbt, print_type):
    if print_sbt:
        if print_type == "Bits":
            print(sbt.graphviz_bits())
        if print_type == "Names":
            print(sbt.graphviz_names())


# Save experiment results into csv
def save_experiment_results(dictionary, benchmark_name, pandas_location):
    pd.DataFrame(dictionary).to_csv(pandas_location + str(benchmark_name).replace('<', '').replace('>', '') + '.csv')
    print()
    print()


# Main.py as a function so that we can pipeline our experiments
def main(p):
    # Report Parameters
    print_params(p)

    # Create SBT
    sbt = SBT(k=p["k"], bloom_filter_length=p["bloom_filter_length"], hash_functions=p["hash_functions"],
              threshold=p["threshold"], similarity_function=p["similarity_function"], sbt_type=p["sbt_type"],
              hash_fraction=p["hash_fraction"])

    # Read Sequences
    sequences = read_sequences(file_names=[p['sequence_prefix'] + str(i) for i in range(p["num_sequences"])],
                               sequence_len=p["sequence_len"], dictionary=p)

    # Insert sequences into SBT
    insert_sequences(sbt=sbt, sequences=sequences, bits_to_check=p["bits_to_check"], method=p["insert_method"],
                     dictionary=p)

    # Query from SBT and report results
    query_sequences(sbt=sbt, all_sequences=sequences, method=p["query_method"], num_queries=p["num_queries"],
                    dictionary=p, query_size=p["query_size"], boyer_moore=p["boyer_moore"])

    # Save SBT
    save_sbt(sbt=sbt, file_name=p["sbt_location"] + "sbt_" + p["sbt_type"],
             dictionary=p)

    # Print Graph Itself
    print_graph(sbt, p["print_sbt"], p["print_type"])

    # Save experiment results
    save_experiment_results(p, p["benchmark_name"], p["pandas_location"])


hash_dict = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}


# Hash function (generates between [0, 2^32])
def hash_lcg(s: str):
    x = 0
    for c in s[:16]:
        x << 2
        x += hash_dict[c]
    return x


# Bit Similarity Function (Small pertubations to break ties)
def hamming(a, b):
    return -sum(a ^ b) + random.random() * 1e-9


# Bit Similarity Function (hamming to break ties)
def and_hamming(a, b):
    return sum(a & b) - sum(a ^ b) * 1e-9


# Bit Similarity Function
def cosine(a, b):
    d = math.sqrt(sum(a) * sum(b))
    return (sum(a & b) / d if d else 0) + random.random() * 1e-9


# Bit Similarity Function
def jaccard(a, b):
    return 1 - sum(a & b) / sum(a | b) + random.random() * 1e-9


# Bit Similarity Function
def manhattan(a, b):
    return sum(a) + sum(b) - 2 * sum(a & b) + random.random() * 1e-9


# Bit Similarity Function
def euclidian(a, b):
    return math.sqrt(manhattan(a, b)) + random.random() * 1e-9


# Bit Similarity Function
def dice(a, b):
    return 2 * sum(a & b) / (sum(a) + sum(b)) + random.random() * 1e-9


# Bit Similarity Function
def tanimoto(a, b):
    c = sum(a & b)
    return c / (sum(a) + sum(b) + c) + random.random() * 1e-9
