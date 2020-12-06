from utils import *
# https://github.com/pythonprofilers/memory_profiler

# Parameters
p = {
    "bloom_filter_length": 300000,            # m
    "sequence_len": 100000,
    "k": 30,                                  # k
    "query_size": 1000,                       # size of query
    "bits_to_check": 200,                     # b'
    "num_sequences": 100,                     # n
    "num_queries": 1000,
    "threshold": 0.9,                         # theta

    "node_class": "Base",                     # SBT Type ("Base", "SSBT")
    "insert_method": "Cluster2",              # ("Greedy", "Cluster1", "Cluster2")
    "query_method": "Normal",                 # ("Normal", "Fast", "Faster")

    "similarity_function": cosine,            # similarity()
    "hash_functions": [hash],                 # h
    "hash_fraction": 1                         # partial hash function simulation
}

# Report Parameters
print_params(p)
# Create SBT
sbt = SBT(p["k"], p["bloom_filter_length"], p["hash_functions"], p["threshold"], p["similarity_function"],
          p["node_class"])
# Read Sequences
sequences = read_sequences(file_names=['testData/sequence' + str(i) for i in range(p["num_sequences"])],
                           sequence_len=p["sequence_len"])
# Insert sequences into SBT
insert_sequences(sbt, sequences, p["bits_to_check"], p["insert_method"])
# Query from SBT and report results
query_sequences(sbt, {name: sequence[-p["query_size"]:] for name, sequence in sequences.items()}, sequences,
                p["query_method"], int(p["num_queries"] / p["num_sequences"]))
# Save SBT
save_sbt(sbt, "testData/sbt_" + p["node_class"])
# Print Expected Bytes
print("Reads Size (Bytes)  ", int(p["sequence_len"] * p["num_sequences"] / 4))
# Print Graph Itself
# print(sbt.graphviz_names())
