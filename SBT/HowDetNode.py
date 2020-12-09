from bitarray import bitarray


class HowDetNode(object):
    count = 0  # How many Nodes have been created

    def __init__(self, bloom_filter_length, hash_functions, similarity_function, experiment_name, how_filter=None):
        self.bloom_filter_length = bloom_filter_length
        self.hash_function = hash_functions[0]  # Only need 1 hash function
        self.similarity_function = similarity_function
        self.experiment_name = experiment_name
        # Initialize Bloom Filters to all 0s
        self.how_filter = bitarray('0') * bloom_filter_length if how_filter is None else how_filter
        self.det_filter = None  # None if node is a leaf
        self.union_filter = None  # None if node is a leaf
        self.left_child = None
        self.right_child = None
        # Give node an id
        self.id = str(HowDetNode.count)
        HowDetNode.count += 1

    """ Creates a new parent Node whose children are left_child and right_child. The Node's filter(s) are set so that
    the SBT Topology/Relationship between nodes is maintained """
    @staticmethod
    def from_children(left_child, right_child):
        # Create new node
        node = HowDetNode(left_child.bloom_filter_length, [left_child.hash_function], left_child.similarity_function,
                          left_child.experiment_name)
        node.experiment_name = "I" + str(node.id)  # Label inner nodes
        # Set new node's filters
        if left_child.union_filter is not None:
            node.union_filter = left_child.union_filter.copy()
        else:
            node.union_filter = left_child.how_filter.copy()
        if right_child.union_filter is not None:
            node.union_filter |= right_child.union_filter
        else:
            node.union_filter |= right_child.how_filter
        node.how_filter = left_child.how_filter & right_child.how_filter
        node.det_filter = node.how_filter | ~node.union_filter
        # Set new node's children
        node.left_child = left_child
        node.right_child = right_child
        return node

    """ Insert a kmer into the Node's how filter """
    def insert_kmer(self, kmer):
        self.how_filter[self.hash_function(kmer) % self.bloom_filter_length] = True

    """ Query a kmer from the Node's determined filter """
    def query_kmer_det(self, kmer):
        return self.det_filter[self.hash_function(kmer) % self.bloom_filter_length]

    """ Query a kmer from the Node's how filter """
    def query_kmer_how(self, kmer):
        return self.how_filter[self.hash_function(kmer) % self.bloom_filter_length]

    """ Return similarity between the first (# bits_to_check) bits of this Node's how filter and the first 
    (# bits_to_check) of another node's how filter """
    def similarity(self, node, bits_to_check=None):
        if bits_to_check is None:
            return self.similarity_function(self.how_filter, node.how_filter)
        return self.similarity_function(self.how_filter[:bits_to_check], node.how_filter[:bits_to_check])

    """ Deep copy fields of node (except for left and right children) """
    def copy(self):
        return HowDetNode(self.bloom_filter_length, [self.hash_function], self.similarity_function, self.experiment_name,
                          self.how_filter.copy())

    """ Insert a single node to an existing SBT greedily by traversing down the most similar child starting from the 
        root """
    def insert_experiment(self, node):
        # 0 children - copy current node into left child and insert into right child
        if self.left_child is None:
            self.left_child = self.copy()
            self.experiment_name = "I" + str(self.id)  # Label inner nodes
            self.right_child = node
            self.union_filter = self.left_child.how_filter | self.right_child.how_filter
            self.how_filter &= self.right_child.how_filter
            self.det_filter = self.how_filter | ~self.union_filter
        # 2 children - iterate into the more similar child
        else:
            # Update filter before iterating onto similar child
            self.union_filter |= node.how_filter
            self.how_filter &= node.how_filter
            self.det_filter = self.how_filter | ~self.union_filter
            # Iterate into more similar child
            left_similarity = self.left_child.similarity(node)
            right_similarity = self.right_child.similarity(node)
            if left_similarity > right_similarity:
                self.left_child.insert_experiment(node)
            else:
                self.right_child.insert_experiment(node)

    """ Query a list of kmers from a SBT by checking whether the respective bit is turned on in the bloom filter. If at
     least (# absolute_threshold) kmers are present, then the query returns all children nodes. If at least |kmers| - 
      absolute_threshold kmers are not present, then the subtree at this node is pruned from search. Lastly, if neither
      of those two conditions are met, then the query proceeds to the children. """
    def query_experiment(self, kmers: list, absolute_threshold):
        partial_hits = []
        complete_hits = 0
        complete_misses = 0
        # If the node is a leaf then only check how filter
        if self.det_filter is None:
            for kmer in kmers:  # Check if each kmer is present
                if self.query_kmer_how(kmer):  # Complete Hit
                    complete_hits += 1
                    if complete_hits >= absolute_threshold:  # Enough hits to return all descendants
                        return self.iter_children()
                else:  # Complete Miss
                    complete_misses += 1
                    if complete_misses > len(kmers) - absolute_threshold:  # Stop since too many misses
                        return []
        # If the node is an inner node, then first check det before checking how
        for kmer in kmers:  # Check if each kmer is present, partially present, or absent
            if self.query_kmer_det(kmer):  # Complete hit or complete miss - all descendants have or have not
                if self.query_kmer_how(kmer):  # Complete Hit
                    complete_hits += 1
                    if complete_hits >= absolute_threshold:  # Enough hits to return all descendants
                        return self.iter_children()
                else:  # Complete Miss
                    complete_misses += self.query_kmer_how(kmer)
                    if complete_misses > len(kmers) - absolute_threshold:  # Stop since too many misses
                        return []
            else:  # Partial hit: some descendants have, some don't
                partial_hits.append(kmer)
        # Search children since not enough hits but not enough misses only on kmer partial hits
        return self.left_child.query_experiment(partial_hits, absolute_threshold - complete_hits) + \
            self.right_child.query_experiment(partial_hits, absolute_threshold - complete_hits)

    """ Faster way to query a list of kmers from a SBT by only hashing the kmers once and then checking a list of 
    filter_indices that the kmers hash to """
    def fast_query_experiment(self, filter_indices, absolute_threshold):
        partial_hits = []
        complete_hits = 0
        complete_misses = 0
        # If the node is a leaf then only check how filter
        if self.det_filter is None:
            for index in filter_indices:  # Check if each index is a hit
                if self.how_filter[index]:  # Complete Hit
                    complete_hits += 1
                    if complete_hits >= absolute_threshold:  # Enough hits to return all descendants
                        return self.iter_children()
                else:  # Complete Miss
                    complete_misses += 1
                    if complete_misses > len(filter_indices) - absolute_threshold:  # Stop since too many misses
                        return []
        for index in filter_indices:  # Check if each index is a hit
            if self.det_filter[index]:  # Complete hit or complete miss - all descendants have or don't have
                if self.how_filter[index]:  # Complete Hit
                    complete_hits += 1
                    if complete_hits >= absolute_threshold:  # Enough hits to return all descendants
                        return self.iter_children()
                else:  # Complete Miss
                    complete_misses += 1
                    if complete_misses > len(filter_indices) - absolute_threshold:  # Stop since too many misses
                        return []
            else:  # Partial hit - some descendants have, some don't
                partial_hits.append(index)
        # Search children since not enough hits but not enough misses only on kmer partial hits
        return self.left_child.fast_query_experiment(partial_hits, absolute_threshold - complete_hits) + \
            self.right_child.fast_query_experiment(partial_hits, absolute_threshold - complete_hits)

    """ Returns a list of the names of all descendant nodes """
    def iter_children(self):
        if self.left_child is not None:
            return self.left_child.iter_children() + self.right_child.iter_children()
        return [self.experiment_name]

    """ Print experiment name and the bits of the bloom filter, then call print on children """
    def print(self):
        print(self.experiment_name, '\t',
              ''.join(map(str, map(int, self.how_filter))), '\n',
              ''.join(map(str, map(int, self.det_filter))) if self.det_filter is not None else '', '\n',
              ''.join(map(str, map(int, self.union_filter))) if self.union_filter is not None else '')
        if self.left_child is not None:
            self.left_child.print()
            self.right_child.print()

    """ Obtain experiment name (bits=False) or bits (bits=True) of the bloom filter, then add to a graphviz Graph, then 
    iterate on children """
    def graphviz(self, graph, bits):
        graph.node(self.id,
                   ''.join(map(str, map(int, self.how_filter))) + '\n'
                   + (''.join(map(str, map(int, self.det_filter)))
                      if self.det_filter is not None else '') + '\n'
                   + (''.join(map(str, map(int, self.union_filter)))
                      if self.union_filter is not None else '')
                   if bits else self.experiment_name)
        if self.left_child is not None:
            graph.edge(tail_name=self.id, head_name=self.left_child.id)
            graph.edge(tail_name=self.id, head_name=self.right_child.id)
            self.left_child.graphviz(graph, bits)
            self.right_child.graphviz(graph, bits)
