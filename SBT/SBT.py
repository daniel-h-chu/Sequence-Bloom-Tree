"""
Sequence Bloom Tree implementation. The SBT relies on one of two Node implementation (i.e. using BaseNode will result
in a basic SBT awhile using SSBTNode will result in a Split-SBT). The SBT contains functions for insertion (and
different insertion algorithms), querying (and different querying algorithms), and I/O (printing, saving, loading)
"""
from graphviz import Digraph
from SBT.SSBTNode import SSBTNode
from SBT.BaseNode import BaseNode
import pickle
import numpy as np


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

    """ Creates a SBT Node from a sequence by breaking down the sequence into kmers and then inserting the kmers using
     the node's implemented insert_kmer() method. If hash_fraction < 1, then some kmers are randomly chosen to not be
     inserted. The node also is labeled with the experiment_name """
    def node_from_sequence(self, sequence: str, experiment_name):
        node = self.NodeClass(self.bloom_filter_length, self.hash_functions, self.similarity_function, experiment_name)

        if self.hash_fraction == 1:  # Just insert all kmers
            for kmer_index in range(0, len(sequence) - self.k + 1):  # Iterate through k-mers
                kmer = sequence[kmer_index:kmer_index + self.k]
                node.insert_kmer(kmer)
            return node
        # Hash only some of the kmers
        kmers_to_insert = np.random.random(len(sequence) - self.k + 1) < self.hash_fraction
        for kmer_index in range(0, len(sequence) - self.k + 1):  # Iterate through k-mers
            if kmers_to_insert[kmer_index]:
                kmer = sequence[kmer_index:kmer_index + self.k]
                node.insert_kmer(kmer)
        return node

    """ Creates a node for a single sequence and inserts it into the SBT using the given experiment_name """
    def insert_sequence(self, sequence: str, experiment_name: str):
        self.insert_node(self.node_from_sequence(sequence, experiment_name))

    """ Insert a pre-generated node into the SBT """
    def insert_node(self, node):
        if self.root is None:
            self.root = node
        else:
            self.root.insert_experiment(node)

    """ Clustering Method 1"""
    """ Inserts a list of sequences using clustering heuristics described in the AllSome Paper. Essentially, we first
    create an SBT of one node for each sequence. Then we check the first (bits_to_check) bits of each sequence and 
    compute pairwise similarities between all sequences. The SBT of sequences with the highest similarities are then 
    joined together as children of a new parent node. The parent node goes back into the group of SBTs we compute the
    pairwise similarities on and can themselves be joined and parented. We repeat until there is only one SBT remaining.
     At that point, the last SBT remaining becomes the root node. """
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
            # Select most similar pair
            for idx1 in range(len(nodes)):
                for idx2 in range(idx1 + 1, len(nodes)):
                    similarity = nodes[idx1].similarity(nodes[idx2], bits_to_check)
                    if similarity > max_similarity:
                        max_similarity = similarity
                        max_pair = (idx1, idx2)
            # Parent the most similar pair
            node = self.NodeClass.from_children(nodes[max_pair[0]], nodes[max_pair[1]])
            nodes.remove(node.left_child)
            nodes.remove(node.right_child)
            nodes.append(node)
        self.root = nodes[0]

    """ Clustering Method 2"""
    """ Inserts a list of sequences using clustering heuristics that were created by us. This heuristic is similar to
    the AllSome heuristic, but the main difference is that we don't allow newly combined SBTs to be considered in
    similarity calculations until all the nodes have been paired once. Then we pair together the parents until all the
    parents have been paired once. We continue until we remain with one node. This ensures the height of the SBT is
    reasonable and also runs faster than the first method """
    def insert_cluster_sequences2(self, sequences: list, experiment_names: list, bits_to_check):
        nodes = []
        if self.root is not None:
            nodes.append(self.root)
        for sequence, name in zip(sequences, experiment_names):
            nodes.append(self.node_from_sequence(sequence, name))
        # Iterate through all nodes, select the two that are the most similar and then create a parent node from them
        while len(nodes) > 1:
            similarities = [[0] * len(nodes) for _ in range(len(nodes))]
            # Compute similairities
            for idx1 in range(len(nodes)):
                for idx2 in range(idx1 + 1, len(nodes)):
                    similarities[idx1][idx2] = nodes[idx1].similarity(nodes[idx2], bits_to_check)
                    similarities[idx2][idx1] = similarities[idx1][idx2]
            parent_nodes = []
            unmatched = {i for i in range(len(nodes))}
            # Matches up pairs until there is 0 or 1 remaining
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
            # Assign parents to now be matched
            nodes = [nodes[idx] for idx in unmatched]
            nodes.extend(parent_nodes)
        self.root = nodes[0]

    """ Generic SBT querying algorithm. This involves checking each kmer as we walk down the tree. """
    def query_sequence(self, sequence: str):
        # Break sequence into individual kmers
        kmers = [sequence[kmer_index:kmer_index + self.k] for kmer_index in range(0, len(sequence) - self.k + 1)]
        # Determine absolute threshold (theta * # kmers) and begin query
        return self.root.query_experiment(kmers=kmers, absolute_threshold=self.threshold * len(kmers))

    """ Fast querying algorithm. We keep track of the indices that the kmers hash so that we don't
     have to hash our kmers every time we search a node. This only works when we have 1 or fewer hash functions. """
    def fast_query_sequence(self, sequence: str):
        if len(self.hash_functions) > 1:
            raise ValueError("Cannot use query method if more than 1 hash function is employed")
        # Determine what indices kmers get mapped to
        filter_indices = []
        for kmer_index in range(0, len(sequence) - self.k + 1):
            kmer = sequence[kmer_index:kmer_index + self.k]
            filter_indices.append(self.hash_functions[0](kmer) % self.bloom_filter_length)
        return self.root.fast_query_experiment(filter_indices=filter_indices,
                                               absolute_threshold=self.threshold * (len(sequence) - self.k + 1))

    """ Print the experiment names and bits of every node in the SBT """
    def print(self):
        self.root.print()

    """ Obtain the experiment names of every node in the SBT in graphviz format """
    def graphviz_names(self):
        graph = Digraph()
        self.root.graphviz(graph=graph, bits=False)
        return graph

    """ Obtain the experiment names of every node in the SBT in graphviz format """
    def graphviz_bits(self):
        graph = Digraph()
        self.root.graphviz(graph=graph, bits=True)
        return graph

    """ Save SBT to a pickle file """
    def save(self, file_name):
        pickle.dump(self, open(file_name, "wb"))

    """ Load SBT from a pickle file """
    @staticmethod
    def load(file_name):
        return pickle.load(open(file_name, "rb"))
