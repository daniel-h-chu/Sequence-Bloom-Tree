"""
Main file from which benchmarking, SBT insertions, SBT queries, and parameterization can be done. The file contains only
utility functions that cleanly read test sequences, insert them into the SBT, query, and report benchmarking statistics.
Nothing here should be changed other than the parameters listed in dictionary p. Feel free to comment in/out one of the
last two lines, which print graphviz versions of the SBT.
"""
from utils import *

# Parameters
p = {
    "bloom_filter_length": 1000000,         # m - Size of bloom filters
    "k": 25,                                # k - Size of kmer
    "bits_to_check": 1000,                  # b' - Number of bits to check when clustering for insertions
    "num_sequences": 250,                   # n - Number of sequences to insert
    "threshold": 0.9,                       # theta - Proportion of kmers that must hit in order to return a node

    "sequence_len": 1000000,                # Size of each sequence inserted
    "query_size": 500,                      # Size of query sequence
    "num_queries": 500,                     # Number of queries to perform

    "sbt_type": "Base",                     # SBT Type ("Base", "SSBT", "HowDe")
    "insert_method": "Cluster2",            # SBT Insertion Method - ("Greedy", "Cluster1", "Cluster2")
    "query_method": "Fast",                 # SBT Query Method - ("Normal", "Fast")

    "similarity_function": hamming,         # Similarity metric to compare filters - (hamming, cosine, jaccard, etc)
    "hash_functions": [hash],               # h - Function to hash kmers
    "hash_fraction": 1,                     # Simulate partial hash function

    "print_sbt": False,                     # Print SBT graph
    "print_type": "Bits",                   # What to print in SBT nodes - ("Bits", "Names")

    "sequence_prefix": "fasta/sim",         # Prefix of genome files (e.g. test_data/sequence1, test_data/sequence2)
    "pandas_location": "experiment_results/",  # Whether to save to pandas file
    "sbt_location": "sbt_data/",            # Local to store sbt
    "benchmark_name": "test_benchmark",     # Name of benchmark

    "boyer_moore": False,                   # Use Boyer-Moore to benchmark against SBT and to verify hits
}

main(p)
