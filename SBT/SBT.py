from graphviz import Digraph
from SBT.SSBTNode import SSBTNode
from SBT.BaseNode import BaseNode
import pickle
import numpy as np
from collections import defaultdict


# Sequence Bloom Tree
# Holds parameters for Bloom Filters and querying, also holds the root Node/Bloom Filter
class SBT(object):
    def __init__(self, k, bloom_filter_length, hash_functions, threshold, similarity_function, node_class="Base",
                 hash_fraction=1):
        self.k = k
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.threshold = threshold
        self.similarity_function = similarity_function
        if node_class is not "Base" and node_class is not "SSBT":
            raise ValueError("Node class should be Base or SSBT")
        self.NodeClass = SSBTNode if node_class is "SSBT" else BaseNode
        self.hash_fraction = hash_fraction
        self.root = None

    def node_from_sequence(self, sequence: str, experiment_name):
        node = self.NodeClass(self.bloom_filter_length, self.hash_functions, self.similarity_function, experiment_name)

        if self.hash_fraction == 1:
            for kmer_index in range(0, len(sequence) - self.k + 1):  # Iterate through k-mers
                kmer = sequence[kmer_index:kmer_index + self.k]
                node.insert_kmer(kmer)
            return node
        to_hash = np.random.random(len(sequence) - self.k + 1) < self.hash_fraction
        for kmer_index in range(0, len(sequence) - self.k + 1):  # Iterate through k-mers
            if to_hash[kmer_index]:
                kmer = sequence[kmer_index:kmer_index + self.k]
                node.insert_kmer(kmer)
        return node

    # Creates a node for a single sequence and inserts it into the SBT
    def insert_sequence(self, sequence: str, experiment_name: str):
        self.insert_node(self.node_from_sequence(sequence, experiment_name))

    # Insert a pre-generated node into the SBT
    def insert_node(self, node):
        if self.root is None:
            self.root = node
        else:
            self.root.insert_experiment(node)

    # Insert multiple sequences using clustering (AllSome Paper)
    def insert_cluster_sequences1(self, sequences: list, experiment_names: list, bits_to_check):
        nodes = []
        if self.root is not None:
            nodes.append(self.root)
        for sequence, name in zip(sequences, experiment_names):
            nodes.append(self.node_from_sequence(sequence, name))
        # Iterate through all nodes, select the two that are the most similar and then create a parent node from them
        while len(nodes) > 1:
            max_similarity = -np.inf
            max_pair = ()
            for idx1 in range(len(nodes)):
                for idx2 in range(idx1 + 1, len(nodes)):
                    similarity = nodes[idx1].similarity(nodes[idx2], bits_to_check)
                    if similarity > max_similarity:
                        max_similarity = similarity
                        max_pair = (idx1, idx2)
            node = self.NodeClass.from_children(nodes[max_pair[0]], nodes[max_pair[1]])
            nodes.remove(node.left_child)
            nodes.remove(node.right_child)
            nodes.append(node)
        self.root = nodes[0]

    # Insert multiple sequences using clustering (AllSome Paper)
    def insert_cluster_sequences2(self, sequences: list, experiment_names: list, bits_to_check):
        nodes = []
        if self.root is not None:
            nodes.append(self.root)
        for sequence, name in zip(sequences, experiment_names):
            nodes.append(self.node_from_sequence(sequence, name))
        # Iterate through all nodes, select the two that are the most similar and then create a parent node from them
        while len(nodes) > 1:
            similarities = [[0] * len(nodes) for _ in range(len(nodes))]
            for idx1 in range(len(nodes)):
                for idx2 in range(idx1 + 1, len(nodes)):
                    similarities[idx1][idx2] = nodes[idx1].similarity(nodes[idx2], bits_to_check)
                    similarities[idx2][idx1] = similarities[idx1][idx2]
            parent_nodes = []
            unmatched = {i for i in range(len(nodes))}
            while len(unmatched) > 1:
                max_similarity = -np.inf
                max_pair = ()
                for idx1 in unmatched:
                    for idx2 in unmatched:
                        if idx1 == idx2:
                            continue
                        if similarities[idx1][idx2] > max_similarity:
                            max_similarity = similarities[idx1][idx2]
                            max_pair = (idx1, idx2)
                node = self.NodeClass.from_children(nodes[max_pair[0]], nodes[max_pair[1]])
                unmatched.remove(max_pair[0])
                unmatched.remove(max_pair[1])
                parent_nodes.append(node)
            nodes = [nodes[idx] for idx in unmatched]
            nodes.extend(parent_nodes)
        self.root = nodes[0]

    # Search the SBT for a single sequence/read
    def query_sequence(self, sequence: str):
        # Break sequence into individual kmers
        kmers = [sequence[kmer_index:kmer_index + self.k] for kmer_index in range(0, len(sequence) - self.k + 1)]
        # Determine absolute threshold (theta * # kmers) and begin query
        return self.root.query_experiment(kmers=kmers, absolute_threshold=self.threshold * len(kmers))

    # Faster lookup if 1 or fewer hash functions are used
    def fast_query_sequence(self, sequence: str):
        if len(self.hash_functions) > 1:
            raise ValueError("Cannot use query method if more than 1 hash function is employed")
        filter_index_dict = defaultdict(int)
        for kmer_index in range(0, len(sequence) - self.k + 1):
            kmer = sequence[kmer_index:kmer_index + self.k]
            filter_index_dict[self.hash_functions[0](kmer) % self.bloom_filter_length] += 1
        return self.root.fast_query_experiment(filter_index_dict=filter_index_dict,
                                               filter_indices=filter_index_dict.keys(),
                                               absolute_threshold=self.threshold * (len(sequence) - self.k + 1),
                                               total_kmers=len(sequence) - self.k + 1)

    # Even faster lookup if 1 or fewer hash functions are used (at the cost of lower precision)
    def faster_query_sequence(self, sequence: str):
        if len(self.hash_functions) > 1:
            raise ValueError("Cannot use query method if more than 1 hash function is employed")
        filter_indices = []
        for kmer_index in range(0, len(sequence) - self.k + 1):
            kmer = sequence[kmer_index:kmer_index + self.k]
            filter_indices.append(self.hash_functions[0](kmer) % self.bloom_filter_length)
        return self.root.faster_query_experiment(filter_indices=filter_indices,
                                                 absolute_threshold=self.threshold * (len(sequence) - self.k + 1))

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
