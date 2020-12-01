from SBT import SBT
from BloomFilter import BloomFilter
import math


def hamming(a, b):
    return sum([x ^ y for x, y in zip(a, b)])


def cosine(a, b):
    return sum([x * y for x, y in zip(a, b)]) / \
           math.sqrt(sum([x * x for x in a]) * sum([y * y for y in b]))


# Parameters
num_files = 100
str_length = 1000

k = 25

bloom_filter_length = 10
hash_functions = [hash]
threshold = 0.5
similarity_function = hamming

sbt = SBT(bloom_filter_length, hash_functions, threshold, similarity_function)

for file_num in range(num_files):
    # Read test text and create bloom filter
    f = open('testData/string' + str(file_num), 'r')
    string = f.readline()
    bloom_filter = BloomFilter(bloom_filter_length, hash_functions)

    for kmer_index in range(0, len(string) - k + 1):
        kmer = string[kmer_index:kmer_index + k]
        bloom_filter.insert(kmer)

    sbt.insert(bloom_filter)
