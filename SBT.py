from Node import Node
from graphviz import Digraph
import pickle


# Sequence Bloom Tree
# Holds parameters for Bloom Filters and querying, also holds the root Node/Bloom Filter
class SBT(object):
    def __init__(self, k, bloom_filter_length, hash_functions, threshold, similarity_function):
        self.k = k
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.threshold = threshold
        self.similarity_function = similarity_function
        self.root = None

    # Creates a node for an experiment (list of sequences) and inserts it into the SBT
    def insert_experiment(self, sequences: list, experiment_name=None):
        node = Node(self.bloom_filter_length, self.hash_functions, self.similarity_function, experiment_name)

        for sequence in sequences:  # Iterate through sequences
            for kmer_index in range(0, len(sequence) - self.k + 1):  # Iterate through k-mers
                kmer = sequence[kmer_index:kmer_index + self.k]
                node.insert_kmer(kmer)

        self.insert_node(node)

    # Creates a node for a single sequence and inserts it into the SBT
    def insert_sequence(self, sequence: str, experiment_name=None):
        node = Node(self.bloom_filter_length, self.hash_functions, self.similarity_function, experiment_name)

        for kmer_index in range(0, len(sequence) - self.k + 1):  # Iterate through k-mers
            kmer = sequence[kmer_index:kmer_index + self.k]
            node.insert_kmer(kmer)

        self.insert_node(node)

    # Insert a pre-generated node into the SBT
    def insert_node(self, node):
        if self.root is None:
            self.root = node
        else:
            self.root.insert_experiment(node)

    # Search the SBT for a single sequence/read
    def query_sequence(self, sequence: str):
        # Break sequence into individual kmers
        kmers = [sequence[kmer_index:kmer_index + self.k] for kmer_index in range(0, len(sequence) - self.k + 1)]
        # Determine absolute threshold (theta * # kmers) and begin query
        return self.root.query_experiment(kmers=kmers, absolute_threshold=self.threshold * len(kmers))

    # Print the experiment names and bits of every node in the SBT
    def print(self):
        self.root.print()

    # Obtain the experiment names of every node in the SBT in graphviz format
    def graphviz_names(self):
        graph = Digraph()
        self.root.graphviz(graph=graph, bits=False)
        return graph

    # Obtain the experiment names of every node in the SBT in graphviz format
    def graphviz_bits(self):
        graph = Digraph()
        self.root.graphviz(graph=graph, bits=True)
        return graph

    def save(self, file_name):
        pickle.dump(self, open(file_name, "wb"))

    @staticmethod
    def load(file_name):
        return pickle.load(open(file_name, "rb"))