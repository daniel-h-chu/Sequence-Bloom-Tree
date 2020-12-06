from SBT.SBT import SBT
from SBT.BaseNode import BaseNode
from utils import *
from bitarray import bitarray

sequence_len = 100000
num_sequences = 100                     # n
k = 300                                 # k
bloom_filter_length = 10                # m
hash_functions = [hash]                 # h
threshold = 0.99                        # theta
similarity_function = hamming           # similarity()
query_size = 400                        # size of query

# Create SBT
sbt = SBT(k, bloom_filter_length, hash_functions, threshold, similarity_function, "BaseNode")
a1 = bitarray([0, 0, 0, 0, 0, 0, 0, 1, 1, 1])
a2 = bitarray([1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
a3 = bitarray([0, 0, 0, 0, 0, 0, 1, 1, 1, 1])
a4 = bitarray([1, 1, 1, 1, 0, 0, 0, 0, 0, 0])
a5 = bitarray([1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
sbt.insert_node(BaseNode(bloom_filter_length, hash_functions, similarity_function, "node", a1))
sbt.insert_node(BaseNode(bloom_filter_length, hash_functions, similarity_function, "node", a2))
sbt.insert_node(BaseNode(bloom_filter_length, hash_functions, similarity_function, "node", a3))
sbt.insert_node(BaseNode(bloom_filter_length, hash_functions, similarity_function, "node", a4))
sbt.insert_node(BaseNode(bloom_filter_length, hash_functions, similarity_function, "node", a5))
print(sbt.graphviz_bits())
