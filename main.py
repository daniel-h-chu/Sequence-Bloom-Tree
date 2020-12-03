from SBT import SBT
import math
import time


# Similarity Functions
def hamming(a, b):
    return sum([x ^ y for x, y in zip(a, b)])


# Similarity Function
def cosine(a, b):
    return sum([x * y for x, y in zip(a, b)]) / \
           math.sqrt(sum([x * x for x in a]) * sum([y * y for y in b]))


# Parameters
sequence_len = 20000
num_sequences = 100                     # n
k = 30                                  # k
bloom_filter_length = 8000              # m
hash_functions = [hash]                 # h
threshold = 0.99                        # theta
similarity_function = hamming           # similarity()
query_size = 200                        # size of query

# Create SBT
sbt = SBT(k, bloom_filter_length, hash_functions, threshold, similarity_function)

# Read sequences
sequences = {}
start_read = time.time()
for i in range(num_sequences):
    # Read test sequences and insert
    f = open('testData/sequence' + str(i), 'r')
    sequences[str(i)] = f.readline()[:sequence_len]
    f.close()

# Insert sequences into SBT
start_insert = time.time()
for name, sequence in sequences.items():
    sbt.insert_sequence(sequence=sequence, experiment_name=name)

# Do stuff: Start here ------------------------------------------------------------

# Print out bit arrays in SBT
# sbt.print()

# Print out Graphviz form of SBT with bits
# print(sbt.graphviz_bits().source)

# Print out Graphviz form of SBT with experiment names
# print(sbt.graphviz_names().source)

# Query from SBT and report results
start_query = time.time()
total_positives = 0
true_positives = 0
for name, sequence in sequences.items():
    results = sbt.query_sequence(sequence=sequence[:query_size])
    total_positives += len(results)
    true_positives += name in results
    # print(name, results)

# Time measurements
end_query = time.time()
for sequence in sequences.values():
    for sequence2 in sequences.values():
        if sequence[:query_size] in sequence2:
            pass
end_inefficient_query = time.time()

print("Read Time", start_insert - start_read)
print("Insert Time", start_query - start_insert)
print("Query Time", end_query - start_query)
print("Inefficient Query", end_inefficient_query - end_query)
print("Precision", true_positives / total_positives)
