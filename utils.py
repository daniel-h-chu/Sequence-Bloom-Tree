import math
import time
import os
from SBT.SBT import SBT
from memory_profiler import profile


def print_params(p):
    for key, value in p.items():
        print(key.ljust(20), value)
    print()


# Read sequences
# @profile
def read_sequences(file_names, sequence_len):
    sequences = {}
    start = time.time()
    for file_name in file_names:
        # Read test sequences and insert
        f = open(file_name, 'r')
        sequences[file_name] = f.readline()[:sequence_len]
        f.close()
    end = time.time()
    print("Read Time           ", end - start)
    return sequences


# Insert sequences into SBT
# @profile
def insert_sequences(sbt, sequences, bits_to_check, method="Greedy"):
    start = time.time()
    if method == "Cluster1":
        sbt.insert_cluster_sequences1(sequences=sequences.values(), experiment_names=sequences.keys(),
                                      bits_to_check=bits_to_check)
    elif method == "Cluster2":
        sbt.insert_cluster_sequences2(sequences=sequences.values(), experiment_names=sequences.keys(),
                                      bits_to_check=bits_to_check)
    else:
        for name, sequence in sequences.items():
            sbt.insert_sequence(sequence=sequence, experiment_name=name)
    end = time.time()
    print("Insert Time         ", end - start)


# Query from SBT and report results
# mode in ("Normal", "Fast", "Faster")
# repeat: number of times to run queries
# @profile
def query_sequences(sbt, q_sequences, all_sequences, method="Normal", repeat=1):
    start = time.time()
    for q_sequence in q_sequences.values():
        for sequence in all_sequences.values():
            if q_sequence in sequence:
                pass
    end = time.time()
    print("Boyer-Moore Time    ", (end - start) * repeat)

    start = time.time()
    total_positives = 0
    true_positives = 0
    false_negatives = 0
    for name, q_sequence in list(q_sequences.items()) * repeat:
        if method == "Faster":
            results = sbt.faster_query_sequence(sequence=q_sequence)
        elif method == "Fast":
            results = sbt.fast_query_sequence(sequence=q_sequence)
        else:
            results = sbt.query_sequence(sequence=q_sequence)
        total_positives += len(results)
        true_positives += name in results
        false_negatives += name not in results
    end = time.time()
    print("Query Time          ", end - start)
    print("Precision           ", true_positives / total_positives)
    print("False Negatives Oops", false_negatives)


# Load SBT and report size
# @profile
def load_sbt(file_name):
    start = time.time()
    sbt = SBT.load(file_name)
    end = time.time()
    print("Load Time           ", end - start)
    return sbt


# Save SBT and report size
# @profile
def save_sbt(sbt, file_name):
    start = time.time()
    sbt.save(file_name)
    end = time.time()
    print("Save Time           ", end - start)
    print("SBT Size (Bytes)    ", os.stat(file_name)[6])


hash_dict = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}


# Hash function (generates between [0, 2^32])
def hash_lcg(s: str):
    x = 0
    for c in s[:16]:
        x << 2
        x += hash_dict[c]
    return x



# Bit Similarity Function
def hamming(a, b):
    return -sum(a ^ b)


# Bit Similarity Function
def cosine(a, b):
    return sum(a & b) / math.sqrt(sum(a) * sum(b))


# Bit Similarity Function
def jaccard(a, b):
    return 1 - sum(a & b) / sum(a | b)


# Bit Similarity Function
def manhattan(a, b):
    return sum(a) + sum(b) - 2 * sum(a & b)


# Bit Similarity Function
def euclidian(a, b):
    return math.sqrt(manhattan(a, b))


# Bit Similarity Function
def dice(a, b):
    return 2 * sum(a & b) / (sum(a) + sum(b))


# Bit Similarity Function
def tanimoto(a, b):
    c = sum(a & b)
    return c / (sum(a) + sum(b) + c)

