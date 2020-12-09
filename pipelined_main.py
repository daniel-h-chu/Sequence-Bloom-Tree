"""
Pipelined main that runs main multiple times but with different parameters
"""
from utils import *
import copy


experiments = [
    # {"key": "bloom_filter_length", "values": [500000, 750000, 1000000, 1500000, 2000000]},
    # {"key": "num_sequences", "values":  [100, 250, 500]},
    # {"key": "query_size", "values":  [250, 500, 1000]},
    # {"key": "hash_fraction", "values":  [0.5, 0.75, 1]},
    # {"key": "threshold", "values":  [0.5, 0.75, 1]},
]

# double_experiments = [
#     {"keys": ("node_class", "insert_method"),
#      "values0": ["Base", "SSBT", "HowDet"],
#      "values1": ["Greedy", "Cluster1", "Cluster2"]},
#     {"keys": ("node_class", "query_method"),
#      "values0": ["Base", "SSBT", "HowDet"],
#      "values1": ["Normal", "Fast"]},
#     {"keys": ("node_class", "similarity_function"),
#      "values0": ["Base", "SSBT", "HowDet"],
#      "values1": [hamming, cosine, jaccard]},
# ]

double_experiments = [
    # {"keys": ("node_class", "query_method"),
    #  "values0": ["HowDet"],
    #  "values1": ["Normal", "Fast"]},
    {"keys": ("node_class", "similarity_function"),
     "values0": ["Base", "SSBT", "HowDet"],
     "values1": [hamming, cosine, jaccard]},
]

# Default Parameters
default_parameters = {
    "bloom_filter_length": 1000000,             # m - Size of bloom filters
    "k": 25,                                   # k - Size of kmer
    "bits_to_check": 1000,                     # b' - Number of bits to check when clustering for insertions
    "num_sequences": 250,                      # n - Number of sequences to insert
    "threshold": 0.9,                          # theta - Proportion of kmers that must hit in order to return a node

    "sequence_len": 1000000,                  # Size of each sequence inserted
    "query_size": 500,                         # Size of query sequence
    "num_queries": 500,                        # Number of queries to perform

    "node_class": "Base",                      # SBT Type ("Base", "SSBT", "HowDet")
    "insert_method": "Cluster2",                 # SBT Insertion Method - ("Greedy", "Cluster1", "Cluster2")
    "query_method": "Fast",                  # SBT Query Method - ("Normal", "Fast")

    "similarity_function": hamming,            # Similarity metric to compare filters - (hamming, cosine, jaccard, etc)
    "hash_functions": [hash],                  # h - Function to hash kmers
    "hash_fraction": 1,                        # Simulate partial hash function

    "print_sbt": False,                        # Print SBT graph
    "print_type": "Bits",                      # What to print in SBT nodes - ("Bits", "Names")

    "sequence_prefix": "fasta/sim",       # Prefix of genome files (e.g. test_data/sequence1, test_data/sequence2)
    "pandas_location": "experiment_results/",  # Whether to save to pandas file
    "sbt_location": "sbt_data/",               # Local to store sbt
    "benchmark_name": "test_benchmark"         # Name of benchmark
}

for experiment in experiments:
    for value in experiment["values"]:
        p = copy.deepcopy(default_parameters)
        p["benchmark_name"] = experiment["key"] + str(value)
        p[experiment["key"]] = value
        main(p)

for experiment in double_experiments:
    for value0 in experiment["values0"]:
        for value1 in experiment["values1"]:
            p = copy.deepcopy(default_parameters)
            p["benchmark_name"] = experiment["keys"][0] + str(value0) + experiment["keys"][1] + str(value1)
            p[experiment["keys"][0]] = value0
            p[experiment["keys"][1]] = value1
            main(p)
