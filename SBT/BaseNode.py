import random
from bitarray import bitarray


# Node
# Holds Bloom Filter parameters
# Each Node holds a single Bloom Filter (boolean array)
class BaseNode(object):
    count = 0  # How many Nodes have been created

    def __init__(self, bloom_filter_length, hash_functions, similarity_function, experiment_name, bloom_filter=None,
                 hash_fraction=1):
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.similarity_function = similarity_function
        self.experiment_name = experiment_name
        self.hash_fraction = hash_fraction
        # Initialize Bloom Filters to all 0s
        self.bloom_filter = bitarray('0') * bloom_filter_length if bloom_filter is None else bloom_filter
        self.left_child = None
        self.right_child = None
        # Give node an id
        self.id = str(BaseNode.count)
        BaseNode.count += 1

    @staticmethod
    def from_children(left_child, right_child):
        node = left_child.copy()
        node.bloom_filter |= right_child.bloom_filter
        node.left_child = left_child
        node.right_child = right_child
        return node

    # Bloom filter function: insert kmer
    def insert_kmer(self, kmer):
        for hash_function in self.hash_functions:
            self.bloom_filter[hash_function(kmer) % self.bloom_filter_length] = True

    # Bloom filter function: query kmer
    def query_kmer(self, kmer):
        for hash_function in self.hash_functions:  # Check if any bits are 0, if so return false
            if not self.bloom_filter[hash_function(kmer) % self.bloom_filter_length]:
                return False
        return True

    # Bloom filter function: return similarity between this node's bloom filter and another node's
    def similarity(self, node, bits_to_check=None):
        if bits_to_check is None:
            return self.similarity_function(self.bloom_filter, node.bloom_filter) + random.random() * 1e-9
        return self.similarity_function(self.bloom_filter[:bits_to_check], node.bloom_filter[:bits_to_check]) + \
            random.random() * 1e-9  # Small pertubation to break ties

    # Deep copy node (except for left and right children)
    def copy(self):
        return BaseNode(self.bloom_filter_length, self.hash_functions, self.similarity_function, self.experiment_name,
                        self.bloom_filter.copy())

    # Insert node based on bloom filter similarity and child presence
    def insert_experiment(self, node):
        # 0 children - copy current node into left child and insert into right child
        if self.left_child is None:
            self.left_child = self.copy()
            self.experiment_name = "I" + str(self.id)  # Label inner nodes
            self.right_child = node
        # 2 children - iterate into the more similar child
        else:
            left_similarity = self.left_child.similarity(node)
            right_similarity = self.right_child.similarity(node)
            if left_similarity > right_similarity:
                self.left_child.insert_experiment(node)
            else:
                self.right_child.insert_experiment(node)
        # Union bloom filter
        self.bloom_filter |= node.bloom_filter

    # Insert node based on presence of kmers and an absolute threshold
    def query_experiment(self, kmers: list, absolute_threshold):
        hits = []
        num_misses = 0
        for kmer in kmers:  # Check if kmer is present
            if self.query_kmer(kmer):
                hits += [kmer]
                if len(hits) >= absolute_threshold:  # Search children since enough hits
                    if self.left_child is None:  # Return since this node is a leaf/exp
                        return [self.experiment_name]
                    return self.left_child.query_experiment(hits, absolute_threshold) + \
                        self.right_child.query_experiment(hits, absolute_threshold)
            else:
                num_misses += 1
                if num_misses > len(kmers) - absolute_threshold:  # Stop since too many misses
                    return []
        print('oops')
        return []  # Should not reach here

    # Faster query (consider a dict of indices and counts to check in each filter)
    def fast_query_experiment(self, filter_index_dict, filter_indices, absolute_threshold, total_kmers):
        filter_indices_hits = []
        num_hits = 0
        num_misses = 0
        for index in filter_indices:  # Check if kmer is present
            if self.bloom_filter[index]:
                filter_indices_hits += [index]
                num_hits += filter_index_dict[index]
                if num_hits >= absolute_threshold:
                    if self.left_child is None:  # Return since this node is a leaf/exp
                        return [self.experiment_name]
                    return self.left_child.fast_query_experiment(filter_index_dict, filter_indices_hits,
                                                                 absolute_threshold, num_hits) + \
                        self.right_child.fast_query_experiment(filter_index_dict, filter_indices_hits,
                                                               absolute_threshold, num_hits)
            else:
                num_misses += filter_index_dict[index]
                if num_misses > total_kmers - absolute_threshold:  # Stop since too many misses
                    return []
        print('oops')
        return []  # Should not reach here

    # Faster query (only consider a set of indices to check in each filter)
    def faster_query_experiment(self, filter_indices, absolute_threshold):
        hits = []
        num_misses = 0
        for index in filter_indices:  # Check if kmer is present
            if self.bloom_filter[index]:
                hits += [index]
                if len(hits) >= absolute_threshold:
                    if self.left_child is None:  # Return since this node is a leaf/exp
                        return [self.experiment_name]
                    return self.left_child.faster_query_experiment(hits, absolute_threshold) + \
                        self.right_child.faster_query_experiment(hits, absolute_threshold)
            else:
                num_misses += 1
                if num_misses > len(filter_indices) - absolute_threshold:  # Stop since too many misses
                    return []
        print('oops')
        return []  # Should not reach here

    # Print experiment name and the bits of the bloom filter, then iterate on children
    def print(self):
        print(self.experiment_name, '\t', ''.join(map(str, map(int, self.bloom_filter))))
        if self.left_child is not None:
            self.left_child.print()
            self.right_child.print()

    # Obtain experiment name or bits of the bloom filter, then add to a graphviz Graph, then iterate on children
    def graphviz(self, graph, bits):
        graph.node(self.id, ''.join(map(str, map(int, self.bloom_filter))) if bits else self.experiment_name)
        if self.left_child is not None:
            graph.edge(tail_name=self.id, head_name=self.left_child.id)
            graph.edge(tail_name=self.id, head_name=self.right_child.id)
            self.left_child.graphviz(graph, bits)
            self.right_child.graphviz(graph, bits)
