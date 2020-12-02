from SBT import SBT
import math


# Similarity Functions
def hamming(a, b):
    return sum([x ^ y for x, y in zip(a, b)])


# Similarity Function
def cosine(a, b):
    return sum([x * y for x, y in zip(a, b)]) / \
           math.sqrt(sum([x * x for x in a]) * sum([y * y for y in b]))


# Parameters
sequence_len = 10000
num_sequences = 100                     # n
k = 30                                  # k
bloom_filter_length = 3000               # m
hash_functions = [hash]                 # h
threshold = 0.99                         # theta
similarity_function = hamming           # similarity()

# Create SBT
sbt = SBT(k, bloom_filter_length, hash_functions, threshold, similarity_function)

# Insert sequences into SBT
for i in range(num_sequences):
    # Read test sequences and insert
    f = open('testData/sequence' + str(i), 'r')
    sbt.insert_sequence(sequence=f.readline()[:sequence_len], experiment_name=str(i))

# Do stuff: Start here ------------------------------------------------------------

# Print out bit arrays in SBT
# sbt.print()

# Print out Graphviz form of SBT with bits
# print(sbt.graphviz_bits().source)

# Print out Graphviz form of SBT with experiment names
# print(sbt.graphviz_names().source)

# Query from SBT and report results
query_size = 200
for i in range(num_sequences):
    # Read test text
    f = open('testData/sequence' + str(i), 'r')
    print(i, sbt.query_sequence(sequence=f.readline()[:query_size]))
