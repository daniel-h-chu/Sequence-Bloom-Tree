from bitarray import bitarray


class BaseNode(object):
    count = 0  # How many Nodes have been created

    def __init__(self, bloom_filter_length, hash_functions, similarity_function, experiment_name, bloom_filter=None):
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.similarity_function = similarity_function
        self.experiment_name = experiment_name
        # Initialize Bloom Filters to all 0s
        self.bloom_filter = bloom_filter if bloom_filter is not None else bitarray('0') * bloom_filter_length
        self.left_child = None
        self.right_child = None
        # Give node an id
        self.id = str(BaseNode.count)
        BaseNode.count += 1

    """ Creates a new parent Node whose children are left_child and right_child. The Node's filter(s) are set so that
    the SBT Topology/Relationship between nodes is maintained """
    @staticmethod
    def from_children(left_child, right_child):
        # Create new node
        node = left_child.copy()
        node.experiment_name = "I" + str(node.id)  # Label inner nodes
        # Set new node filters
        node.bloom_filter |= right_child.bloom_filter
        # Set new node's children
        node.left_child = left_child
        node.right_child = right_child
        return node

    """ Insert a kmer into the Node's bloom filter """
    def insert_kmer(self, kmer):
        for hash_function in self.hash_functions:
            self.bloom_filter[hash_function(kmer) % self.bloom_filter_length] = True

    """ Query a kmer from the Node's bloom filter """
    def query_kmer(self, kmer):
        for hash_function in self.hash_functions:  # Check if any bits are 0, if so return false
            if not self.bloom_filter[hash_function(kmer) % self.bloom_filter_length]:
                return False
        return True

    """ Return similarity between the first (# bits_to_check) bits of this Node's bloom filter and the first 
    (# bits_to_check) of another node's bloom filter """
    def similarity(self, node, bits_to_check=None):
        if bits_to_check is None:
            return self.similarity_function(self.bloom_filter, node.bloom_filter)
        return self.similarity_function(self.bloom_filter[:bits_to_check], node.bloom_filter[:bits_to_check])

    """ Deep copy fields of node (except for left and right children) """
    def copy(self):
        return BaseNode(self.bloom_filter_length, self.hash_functions, self.similarity_function, self.experiment_name,
                        self.bloom_filter.copy())

    """ Insert a single node to an existing SBT greedily by traversing down the most similar child starting from the 
    root """
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

    """ Query a list of kmers from a SBT by checking whether the respective bit is turned on in the bloom filter. If at
     least (# absolute_threshold) kmers are present, then the query proceeds to the children. If the current node is a 
     leaf then the node's name is returned. """
    def query_experiment(self, kmers: list, absolute_threshold):
        hits = []
        num_misses = 0
        for kmer in kmers:  # Check if kmer is present
            if self.query_kmer(kmer):
                hits.append(kmer)
            else:
                num_misses += 1
                if num_misses > len(kmers) - absolute_threshold:  # Stop since too many misses
                    return []
        # Passed threshold
        if self.left_child is None:  # Return since this node is a leaf/exp
            return [self.experiment_name]
        return self.left_child.query_experiment(hits, absolute_threshold) + \
            self.right_child.query_experiment(hits, absolute_threshold)

    """ Faster way to query a list of kmers from a SBT by only hashing the kmers once and then checking a list of 
    filter_indices that the kmers hash to """
    def fast_query_experiment(self, filter_indices, absolute_threshold):
        hits = []
        num_misses = 0
        for index in filter_indices:  # Check if each index is a hit
            if self.bloom_filter[index]:  # Hit
                hits.append(index)
            else:  # Complete miss - none of descendants have a hit at that index
                num_misses += 1
                if num_misses > len(filter_indices) - absolute_threshold:  # Stop since too many misses
                    return []
        # Passed threshold
        if self.left_child is None:  # Return since this node is a leaf
            return [self.experiment_name]
        return self.left_child.fast_query_experiment(hits, absolute_threshold) + \
            self.right_child.fast_query_experiment(hits, absolute_threshold)

    """ Print experiment name and the bits of the bloom filter, then call print on children """
    def print(self):
        print(self.experiment_name, '\t', ''.join(map(str, map(int, self.bloom_filter))))
        if self.left_child is not None:
            self.left_child.print()
            self.right_child.print()

    """ Obtain experiment name (bits=False) or bits (bits=True) of the bloom filter, then add to a graphviz Graph, then 
    iterate on children """
    def graphviz(self, graph, bits):
        graph.node(self.id, ''.join(map(str, map(int, self.bloom_filter))) if bits else self.experiment_name)
        if self.left_child is not None:
            graph.edge(tail_name=self.id, head_name=self.left_child.id)
            graph.edge(tail_name=self.id, head_name=self.right_child.id)
            self.left_child.graphviz(graph, bits)
            self.right_child.graphviz(graph, bits)
