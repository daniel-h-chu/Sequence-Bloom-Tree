from Node import Node


class SBT(object):
    def __init__(self, k, bloom_filter_length, hash_functions, threshold, similarity_function):
        self.k = k
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.threshold = threshold
        self.similarity_function = similarity_function
        self.root = Node(self.bloom_filter_length, self.hash_functions, similarity_function)

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
        self.root.insert_child_node(node)

    def query_sequence(self, node):

        pass


