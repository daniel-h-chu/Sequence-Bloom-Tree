from Node import Node
from graphviz import Digraph


class SBT(object):
    def __init__(self, k, bloom_filter_length, hash_functions, threshold, similarity_function):
        self.k = k
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.threshold = threshold
        self.similarity_function = similarity_function
        self.root = None

    def insert_experiment(self, sequences: list, experiment_name=None):
        node = Node(self.bloom_filter_length, self.hash_functions, self.similarity_function, experiment_name)

        for sequence in sequences:
            for kmer_index in range(0, len(sequence) - self.k + 1):
                kmer = sequence[kmer_index:kmer_index + self.k]
                node.insert_kmer(kmer)

        self.insert_node(node)

    def insert_sequence(self, sequence: str, experiment_name=None):
        node = Node(self.bloom_filter_length, self.hash_functions, self.similarity_function, experiment_name)

        for kmer_index in range(0, len(sequence) - self.k + 1):
            kmer = sequence[kmer_index:kmer_index + self.k]
            node.insert_kmer(kmer)

        self.insert_node(node)

    def insert_node(self, node):
        if self.root is None:
            self.root = node
        else:
            self.root.insert_experiment(node)

    def query_sequence(self, sequence: str):
        kmers = [sequence[kmer_index:kmer_index + self.k] for kmer_index in range(0, len(sequence) - self.k + 1)]
        return self.root.query_experiment(kmers, self.threshold)

    def print(self):
        self.root.print()

    def graphviz_names(self):
        graph = Digraph()
        return self.root.graphviz(graph, False)

    def graphviz_bits(self):
        graph = Digraph()
        return self.root.graphviz(graph, True)
