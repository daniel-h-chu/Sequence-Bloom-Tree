from SBT import SBT
from utils import *
# https://github.com/pythonprofilers/memory_profiler

# Parameters
sequence_len = 100000
num_sequences = 100                    # n
k = 30                                 # k
bloom_filter_length = 300000           # m
hash_functions = [hash]                # h
threshold = 0.4                        # theta
similarity_function = hamming          # similarity()
query_size = 1000                      # size of query
bits_to_check = 200                    # b'
# node_class = "BaseNode"                 # Base SBT
node_class = "SSBTNode"                 # Split SBT

# Create SBT
sbt = SBT(k, bloom_filter_length, hash_functions, threshold, similarity_function, node_class)
# Read Sequences
sequences = read_sequences(file_names=['testData/sequence' + str(i) for i in range(num_sequences)], sequence_len=sequence_len)
# Insert sequences into SBT
# insert_sequences(sbt, sequences)                         # Greedy Method
# insert_cluster_sequences1(sbt, sequences, bits_to_check) # Clustering Method 1
insert_cluster_sequences2(sbt, sequences, bits_to_check)   # Clustering Method 2
# Query from SBT and report results
query_sequences(sbt, {name: sequence[:query_size] for name, sequence in sequences.items()}, sequences)
# Save SBT
save_sbt(sbt, "savedSBT/sbt_" + node_class)
# Print Expected Bytes
print("Reads Size (Bytes)  ", int(sequence_len * num_sequences / 4))

# print(sbt.graphviz_names())