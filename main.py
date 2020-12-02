from SBT import SBT
import math


def hamming(a, b):
    return sum([x ^ y for x, y in zip(a, b)])


def cosine(a, b):
    return sum([x * y for x, y in zip(a, b)]) / \
           math.sqrt(sum([x * x for x in a]) * sum([y * y for y in b]))


# Parameters
num_files = 8                          # n
str_length = 30

k = 25                                  # k
bloom_filter_length = 10                # m
hash_functions = [hash]                 # h
threshold = 0.5                         # theta
similarity_function = hamming           # similarity()

sbt = SBT(k, bloom_filter_length, hash_functions, threshold, similarity_function)

for file_num in range(num_files):
    # Read test text and insert
    f = open('testData/sequence' + str(file_num), 'r')
    sbt.insert_sequence(f.readline()[:str_length], str(file_num))

# sbt.print()
print(sbt.graphviz_bits().source)
# for file_num in range(num_files):
#     # Read test text and insert
#     f = open('testData/sequence' + str(file_num), 'r')
#     print(file_num, sbt.query_sequence(f.readline()[:35]))

