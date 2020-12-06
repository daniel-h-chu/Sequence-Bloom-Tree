from bitarray import bitarray


def hamming(a, b):
    return -sum(a ^ b)


def ssbt_similarity_function(a, b):
    return sum(a & b) + hamming(a, b) * 1e-9


# Node
# Holds Bloom Filter parameters
# Each Node holds a single Bloom Filter (boolean array)
class SSBTNode(object):
    count = 0  # How many Nodes have been created

    def __init__(self, bloom_filter_length, hash_functions, similarity_function, experiment_name, sim_filter=None,
                 hash_fraction=1):
        self.bloom_filter_length = bloom_filter_length
        self.hash_function = hash_functions[0]  # Only need 1 hash function
        self.similarity_function = ssbt_similarity_function  # Small pertubations to break ties
        self.experiment_name = experiment_name
        self.hash_fraction = hash_fraction
        # Initialize Bloom Filters to all 0s
        self.sim_filter = bitarray('0') * bloom_filter_length if sim_filter is None else sim_filter
        self.rem_filter = None
        self.left_child = None
        self.right_child = None
        # Give node an id
        self.id = str(SSBTNode.count)
        SSBTNode.count += 1

    @staticmethod
    def from_children(left_child, right_child):
        node = left_child.copy()
        node.experiment_name = "I" + str(node.id)  # Label inner nodes
        # Set new node filters
        node.sim_filter &= right_child.sim_filter
        node.rem_filter = left_child.sim_filter | right_child.sim_filter
        if left_child.rem_filter is not None:
            node.rem_filter = node.rem_filter | left_child.rem_filter
        if left_child.rem_filter is not None:
            node.rem_filter = node.rem_filter | right_child.rem_filter
        node.rem_filter &= ~node.sim_filter
        # Remove redundant bits
        left_child.sim_filter &= ~node.sim_filter
        right_child.sim_filter &= ~node.sim_filter
        # Assign children
        node.left_child = left_child
        node.right_child = right_child
        return node

    # Bloom filter function: insert kmer
    def insert_kmer(self, kmer):
        self.sim_filter[self.hash_function(kmer) % self.bloom_filter_length] = True

    # Bloom filter function: insert kmer
    def query_kmer_sim(self, kmer):
        return self.sim_filter[self.hash_function(kmer) % self.bloom_filter_length]

    # Bloom filter function: insert kmer
    def query_kmer_rem(self, kmer):
        return self.rem_filter[self.hash_function(kmer) % self.bloom_filter_length]

    # Bloom filter function: return similarity between this node's bloom filter and another node's
    def similarity(self, node, bits_to_check=None):
        if bits_to_check is None:
            return self.similarity_function(self.sim_filter, node.sim_filter)
        return self.similarity_function(self.sim_filter[:bits_to_check], node.sim_filter[:bits_to_check])

    # Deep copy node (except for left and right children)
    def copy(self):
        return SSBTNode(self.bloom_filter_length, [self.hash_function], self.similarity_function, self.experiment_name,
                        self.sim_filter.copy())

    # Insert node based on bloom filter similarity and child presence
    def insert_experiment(self, node):
        # 0 children - copy current node into left child and insert into right child
        if self.left_child is None:
            self.left_child = self.copy()
            self.experiment_name = "I" + str(self.id)  # Label inner nodes
            self.right_child = node
            new_sim_filter = self.sim_filter & node.sim_filter
            self.rem_filter = (self.sim_filter & ~new_sim_filter) | (node.sim_filter & ~new_sim_filter)
            self.sim_filter = new_sim_filter
            # These might not be necessary
            self.left_child.sim_filter &= ~self.sim_filter
            self.right_child.sim_filter &= ~self.sim_filter
        # 2 children - iterate into the more similar child
        else:
            new_sim_filter = self.sim_filter & node.sim_filter  # If node is 1 then sim filter remains 1
            new_rem_filter = self.rem_filter | (self.sim_filter ^ node.sim_filter)  # Not all or nothing filter
            new_node_filter = ~self.sim_filter & node.sim_filter  # If sim filter is 1, then don't have to carry sim
            self.left_child.sim_filter |= (self.sim_filter & ~node.sim_filter)
            self.right_child.sim_filter |= (self.sim_filter & ~node.sim_filter)
            self.sim_filter = new_sim_filter
            self.rem_filter = new_rem_filter
            node.sim_filter = new_node_filter

            left_similarity = self.left_child.similarity(node)
            right_similarity = self.right_child.similarity(node)
            if left_similarity > right_similarity:
                self.left_child.insert_experiment(node)
            else:
                self.right_child.insert_experiment(node)

    # Insert node based on presence of kmers and an absolute threshold
    def query_experiment(self, kmers: list, absolute_threshold):
        partial_hits = []
        complete_hits = 0
        complete_misses = 0
        for kmer in kmers:  # Check if kmer is present
            if self.query_kmer_sim(kmer):  # Complete hit - all children have
                complete_hits += 1
                if complete_hits >= absolute_threshold:  # Enough hits to return all children
                    return self.iter_children()
            elif self.rem_filter is not None and self.query_kmer_rem(kmer):  # Partial hit:some children have,some don't
                partial_hits += [kmer]
            else:  # Complete miss - no children have
                complete_misses += 1
                if complete_misses > len(kmers) - absolute_threshold:  # Stop since too many misses
                    return []
        # Search children since not enough hits but not enough misses only on kmer partial hits
        return self.left_child.query_experiment(partial_hits, absolute_threshold - complete_hits) + \
            self.right_child.query_experiment(partial_hits, absolute_threshold - complete_hits)
    
    # Faster query (consider a dict of indices and counts to check in each filter)
    def fast_query_experiment(self, filter_index_dict, filter_indices, absolute_threshold, total_kmers):
        partial_hits = []
        complete_hits = 0
        complete_misses = 0
        for index in filter_indices:  # Check if kmer is present
            if self.sim_filter[index]:  # Complete hit - all children have
                complete_hits += filter_index_dict[index]
                if complete_hits >= absolute_threshold:  # Enough hits to return all children
                    return self.iter_children()
            elif self.rem_filter is not None and self.rem_filter[index]:  # Partial hit - some children have, some don't
                partial_hits += [index]
            else:  # Complete miss - no children have
                complete_misses += filter_index_dict[index]
                if complete_misses > total_kmers - absolute_threshold:  # Stop since too many misses
                    return []
        # Search children since not enough hits but not enough misses only on kmer partial hits
        return self.left_child.fast_query_experiment(filter_index_dict, partial_hits, absolute_threshold -
                                                     complete_hits, total_kmers - complete_hits - complete_misses) + \
            self.right_child.fast_query_experiment(filter_index_dict, partial_hits, absolute_threshold -
                                                   complete_hits, total_kmers - complete_hits - complete_misses)
        # Faster query (consider a dict of indices and counts to check in each filter)

    # Faster query (only consider a set of indices to check in each filter)
    def faster_query_experiment(self, filter_indices, absolute_threshold):
        partial_hits = []
        complete_hits = 0
        complete_misses = 0
        for index in filter_indices:  # Check if kmer is present
            if self.sim_filter[index]:  # Complete hit - all children have
                complete_hits += 1
                if complete_hits >= absolute_threshold:  # Enough hits to return all children
                    return self.iter_children()
            elif self.rem_filter is not None and self.rem_filter[index]:  # Partial hit - some children have, some don't
                partial_hits += [index]
            else:  # Complete miss - no children have
                complete_misses += 1
                if complete_misses > len(filter_indices) - absolute_threshold:  # Stop since too many misses
                    return []
        # Search children since not enough hits but not enough misses only on kmer partial hits
        return self.left_child.faster_query_experiment(partial_hits, absolute_threshold - complete_hits) + \
            self.right_child.faster_query_experiment(partial_hits, absolute_threshold - complete_hits)

    def iter_children(self):
        if self.left_child is not None:
            return self.left_child.iter_children() + self.right_child.iter_children()
        return [self.experiment_name]

    # Print experiment name and the bits of the bloom filter, then iterate on children
    def print(self):
        print(self.experiment_name, '\t',
              ''.join(map(str, map(int, self.sim_filter))), '\n',
              ''.join(map(str, map(int, self.rem_filter))) if self.rem_filter is not None else '')
        if self.left_child is not None:
            self.left_child.print()
            self.right_child.print()

    # Obtain experiment name or bits of the bloom filter, then add to a graphviz Graph, then iterate on children
    def graphviz(self, graph, bits):
        graph.node(self.id,
                   ''.join(map(str, map(int, self.sim_filter))) + '\n' + (''.join(map(str, map(int, self.rem_filter)))
                                                                          if self.rem_filter is not None else '')
                   if bits else self.experiment_name)
        if self.left_child is not None:
            graph.edge(tail_name=self.id, head_name=self.left_child.id)
            graph.edge(tail_name=self.id, head_name=self.right_child.id)
            self.left_child.graphviz(graph, bits)
            self.right_child.graphviz(graph, bits)
