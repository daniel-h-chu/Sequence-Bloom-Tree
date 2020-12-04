from SBT import SBT
from utils import *
# https://github.com/pythonprofilers/memory_profiler


# Parameters
sequence_len = 100000
num_sequences = 100                     # n
k = 300                                 # k
bloom_filter_length = 40000             # m
hash_functions = [hash]                 # h
threshold = 0.99                        # theta
similarity_function = hamming           # similarity()
query_size = 400                        # size of query

# Create SBT
sbt = SBT(k, bloom_filter_length, hash_functions, threshold, similarity_function)
# Read Sequences
sequences = read_sequences(file_names=['testData/sequence' + str(i) for i in range(num_sequences)])
# Insert sequences into SBT
insert_sequences(sbt, sequences)
# Query from SBT and report results
query_sequences(sbt, {name: sequence[:query_size] for name, sequence in sequences.items()}, sequences)
# Save SBT
save_sbt(sbt, "savedSBT/sbt")
# Print Expected Bytes
print("Reads Size (Bytes)  ", int(sequence_len * num_sequences / 4))
