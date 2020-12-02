import random


# Node
# Holds Bloom Filter parameters
# Each Node holds a single Bloom Filter (boolean array)
class Node(object):
    count = 0  # How many Nodes have been created

    def __init__(self, bloom_filter_length, hash_functions, similarity_function, experiment_name,
                 bloom_filter=None):
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.similarity_function = similarity_function
        self.experiment_name = experiment_name
        # Initialize Bloom Filters to all 0s
        self.bloom_filter = [False] * bloom_filter_length if bloom_filter is None else bloom_filter
        self.left_child = None
        self.right_child = None
        # Give node an id
        self.id = str(Node.count)
        Node.count += 1

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

    # Bloom filter function: union this node's bloom filter with other node bloom filter and store in this bloom filter
    def union_bloom_filter(self, node):
        for i in range(self.bloom_filter_length):
            self.bloom_filter[i] = self.bloom_filter[i] | node.bloom_filter[i]

    # Bloom filter function: return similarity between this node's bloom filter and another node's
    def bloom_filter_similarity(self, node):
        return self.similarity_function(self.bloom_filter, node.bloom_filter)

    # Deep copy node (except for left and right children)
    def copy(self):
        return Node(self.bloom_filter_length, self.hash_functions, self.similarity_function, self.experiment_name,
                    [bit for bit in self.bloom_filter])

    # Insert node based on bloom filter similarity and child presence
    def insert_experiment(self, node):
        # 0 children - copy current node into left child and insert into right child
        if self.left_child is None and self.right_child is None:
            self.left_child = self.copy()
            self.experiment_name = "I" + str(self.id)  # Label inner nodes
            self.right_child = node
        # 1 child - insert into other child
        elif self.left_child is None:
            self.left_child = node
        elif self.right_child is None:
            self.right_child = node
        # 2 children - iterate into the more similar child
        else:
            left_similarity = self.left_child.bloom_filter_similarity(node)
            right_similarity = self.right_child.bloom_filter_similarity(node)
            if left_similarity > right_similarity:
                self.left_child.insert_experiment(node)
            elif right_similarity > left_similarity:
                self.right_child.insert_experiment(node)
            elif random.random() > 0.5:  # Randomly choose direction if equal so that tree is balanced
                self.left_child.insert_experiment(node)
            else:
                self.right_child.insert_experiment(node)
        # Union bloom filter
        self.union_bloom_filter(node)

    # Insert node based on presence of kmers and an absolute threshold
    def query_experiment(self, kmers: list, absolute_threshold):
        hits = []
        num_misses = 0
        for kmer in kmers:  # Check if kmer is present
            if self.query_kmer(kmer):
                hits += [kmer]
            else:
                num_misses += 1
            if len(hits) >= absolute_threshold:  # Search children since enough hits
                if self.left_child is None and self.right_child is None:  # Return since this node is a leaf/exp
                    return [self.experiment_name]
                experiment_hits = []
                if self.left_child is not None:  # Iterate search in children only on kmer hits
                    experiment_hits += self.left_child.query_experiment(hits, absolute_threshold)
                if self.right_child is not None:  # Iterate search in children only on kmer hits
                    experiment_hits += self.right_child.query_experiment(hits, absolute_threshold)
                return experiment_hits
            if num_misses > len(kmers) - absolute_threshold:  # Stop since too many misses
                return []
        print('oops')
        return []  # Should not reach here

    # Print experiment name and the bits of the bloom filter, then iterate on children
    def print(self):
        print(self.experiment_name, '\t', ''.join(map(str,map(int, self.bloom_filter))))
        if self.left_child is not None:
            self.left_child.print()
        if self.right_child is not None:
            self.right_child.print()

    # Obtain experiment name or bits of the bloom filter, then add to a graphviz Graph, then iterate on children
    def graphviz(self, graph, bits):
        graph.node(self.id, ''.join(map(str, map(int, self.bloom_filter))) if bits else self.experiment_name)
        if self.left_child is not None:
            graph.edge(tail_name=self.id, head_name=self.left_child.id)
            self.left_child.graphviz(graph, bits)
        if self.right_child is not None:
            graph.edge(tail_name=self.id, head_name=self.right_child.id)
            self.right_child.graphviz(graph, bits)
