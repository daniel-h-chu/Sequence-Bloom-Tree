"""
Main file from which benchmarking, SBT insertions, SBT queries, and parameterization can be done. The file contains only
utility functions that cleanly read test sequences, insert them into the SBT, query, and report benchmarking statistics.
Nothing here should be changed other than the parameters listed in dictionary p. Feel free to comment in/out one of the
last two lines, which print graphviz versions of the SBT.
"""
from utils import *

# Parameters
p = {
    "bloom_filter_length": 200000,      # m - Size of bloom filters
    "k": 30,                            # k - Size of kmer
    "bits_to_check": 1000,              # b' - Number of bits to check when clustering for insertions
    "num_sequences": 100,               # n - Number of sequences to insert
    "threshold": 0.7,                   # theta - Proportion of kmers that must hit in order to return an experiment

    "sequence_len": 100000,             # Size of each sequence inserted
    "query_size": 1000,                 # Size of query sequence
    "num_queries": 1000,                # Number of queries to perform

    "node_class": "SSBT",               # SBT Type ("Base", "SSBT")
    "insert_method": "Cluster2",        # SBT Insertion Method - ("Greedy", "Cluster1", "Cluster2")
    "query_method": "Fast",             # SBT Query Method - ("Normal", "Fast")

    "similarity_function": cosine,     # Similarity metric to compare filters - (hamming, cosine, jaccard, etc)
    "hash_functions": [hash],           # h - Function to hash kmers (python's hash is the only good hash function)
    "hash_fraction": 1                  # Simulate partial hash function (Note SSBT works well if hash_fraction = 1)
}

# Report Parameters
print_params(p)

# Create SBT
sbt = SBT(k=p["k"], bloom_filter_length=p["bloom_filter_length"], hash_functions=p["hash_functions"],
          threshold=p["threshold"], similarity_function=p["similarity_function"], node_class=p["node_class"],
          hash_fraction=p["hash_fraction"])

# Read Sequences
sequences = read_sequences(file_names=['testData/sequence' + str(i) for i in range(p["num_sequences"])],
                           sequence_len=p["sequence_len"])

# Insert sequences into SBT
insert_sequences(sbt=sbt, sequences=sequences, bits_to_check=p["bits_to_check"], method=p["insert_method"])

# Query from SBT and report results
query_sequences(sbt=sbt, q_sequences={name: sequence[-p["query_size"]:] for name, sequence in sequences.items()},
                all_sequences=sequences, method=p["query_method"], repeat=int(p["num_queries"] / p["num_sequences"]))

# Save SBT
save_sbt(sbt=sbt, file_name="testData/sbt_" + p["node_class"])

# Print Expected Bytes
print("Reads Size (Bytes)  ", int(p["sequence_len"] * p["num_sequences"] / 4))

# Print Graph Itself
# print(sbt.graphviz_names())
# print(sbt.graphviz_bits())
