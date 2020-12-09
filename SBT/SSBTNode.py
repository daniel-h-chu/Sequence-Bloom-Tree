""" Sequence Bloom Tree Node implementation based off of HowDe-SBT in Kingsford & Solomon (2018) """
from bitarray import bitarray


class SSBTNode(object):
    count = 0  # How many Nodes have been created

    def __init__(self, bloom_filter_length, hash_functions, similarity_function, experiment_name, sim_filter=None):
        self.bloom_filter_length = bloom_filter_length
        self.hash_function = hash_functions[0]  # Only need 1 hash function
        self.similarity_function = similarity_function
        self.experiment_name = experiment_name
        # Initialize Bloom Filters to all 0s
        self.sim_filter = bitarray('0') * bloom_filter_length if sim_filter is None else sim_filter
        self.rem_filter = None  # None if node is a leaf
        self.left_child = None
        self.right_child = None
        # Give node an id
        self.id = str(SSBTNode.count)
        SSBTNode.count += 1

    """ Creates a new parent Node whose children are left_child and right_child. The Node's filter(s) are set so that
    the SBT Topology/Relationship between nodes is maintained """
    @staticmethod
    def from_children(left_child, right_child):
        # Create new node
        node = left_child.copy()
        node.experiment_name = "I" + str(node.id)  # Label inner nodes
        # Set new node filters and update child node filters
        node.sim_filter &= right_child.sim_filter
        left_child.sim_filter &= ~node.sim_filter
        right_child.sim_filter &= ~node.sim_filter
        node.rem_filter = left_child.sim_filter | right_child.sim_filter
        if left_child.rem_filter is not None:
            node.rem_filter |= left_child.rem_filter
        if left_child.rem_filter is not None:
            node.rem_filter |= right_child.rem_filter
        # Set new node's children
        node.left_child = left_child
        node.right_child = right_child
        return node

    """ Insert a kmer into the Node's similarity filter """
    def insert_kmer(self, kmer):
        self.sim_filter[self.hash_function(kmer) % self.bloom_filter_length] = True

    """ Query a kmer from the Node's similarity filter """
    def query_kmer_sim(self, kmer):
        return self.sim_filter[self.hash_function(kmer) % self.bloom_filter_length]

    """ Query a kmer from the Node's remainder filter """
    def query_kmer_rem(self, kmer):
        return self.rem_filter[self.hash_function(kmer) % self.bloom_filter_length]

    """ Return similarity between the first (# bits_to_check) bits of this Node's sim filter and the first 
    (# bits_to_check) of another node's sim filter """
    def similarity(self, node, bits_to_check=None):
        if bits_to_check is None:
            return self.similarity_function(self.sim_filter, node.sim_filter)
        return self.similarity_function(self.sim_filter[:bits_to_check], node.sim_filter[:bits_to_check])

    """ Deep copy fields of node (except for left and right children) """
    def copy(self):
        return SSBTNode(self.bloom_filter_length, [self.hash_function], self.similarity_function, self.experiment_name,
                        self.sim_filter.copy())

    """ Insert a single node to an existing SBT greedily by traversing down the most similar child starting from the 
    root """
    def insert_experiment(self, node):
        # 0 children - copy current node into left child and insert into right child
        if self.left_child is None:
            self.left_child = self.copy()
            self.experiment_name = "I" + str(self.id)  # Label inner nodes
            self.right_child = node
            new_sim_filter = self.sim_filter & node.sim_filter
            self.rem_filter = (self.sim_filter & ~new_sim_filter) | (node.sim_filter & ~new_sim_filter)
            self.sim_filter = new_sim_filter
            # Drop bits that are already similar in a parent node
            self.left_child.sim_filter &= ~self.sim_filter
            self.right_child.sim_filter &= ~self.sim_filter
        # 2 children - iterate into the more similar child
        else:
            # Update filter before iterating onto similar child
            new_sim_filter = self.sim_filter & node.sim_filter  # If node is 1 then sim filter remains 1
            new_rem_filter = self.rem_filter | (self.sim_filter ^ node.sim_filter)  # Not all or nothing filter
            new_node_filter = ~self.sim_filter & node.sim_filter  # If sim filter is 1, then don't have to carry sim
            self.left_child.sim_filter |= (self.sim_filter & ~node.sim_filter)
            self.right_child.sim_filter |= (self.sim_filter & ~node.sim_filter)
            self.sim_filter = new_sim_filter
            self.rem_filter = new_rem_filter
            node.sim_filter = new_node_filter
            # Iterate onto more similar child
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
        for kmer in kmers:  # Check if kmer is present
            if self.query_kmer_sim(kmer):  # Complete hit - all descendants have
                complete_hits += 1
                if complete_hits >= absolute_threshold:  # Enough hits to return all descendants
                    return self.iter_children()
            elif self.rem_filter is not None and self.query_kmer_rem(kmer):  # Partial hit : some descendant have
                partial_hits.append(kmer)
            else:  # Complete miss - no descendants have
                complete_misses += 1
                if complete_misses > len(kmers) - absolute_threshold:  # Stop since too many misses
                    return []
        # Search children since not enough hits but not enough misses only on kmer partial hits
        return self.left_child.query_experiment(partial_hits, absolute_threshold - complete_hits) + \
            self.right_child.query_experiment(partial_hits, absolute_threshold - complete_hits)

    """ Faster way to query a list of kmers from a SBT by only hashing the kmers once and then checking a list of 
    filter_indices that the kmers hash to """
    def fast_query_experiment(self, filter_indices, absolute_threshold):
        partial_hits = []
        complete_hits = 0
        complete_misses = 0
        for index in filter_indices:  # Check if each index is a hit
            if self.sim_filter[index]:  # Complete hit - all descendants have
                complete_hits += 1
                if complete_hits >= absolute_threshold:  # Enough hits to return all descendants
                    return self.iter_children()
            elif self.rem_filter is not None and self.rem_filter[index]:  # Partial hit - some descendants have
                partial_hits.append(index)
            else:  # Complete miss - no descendants have
                complete_misses += 1
                if complete_misses > len(filter_indices) - absolute_threshold:  # Stop since too many misses
                    return []
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
        print(self.experiment_name, '\t', 'sim: ',
              ''.join(map(str, map(int, self.sim_filter))), '\n', 'rem: ',
              ''.join(map(str, map(int, self.rem_filter))) if self.rem_filter is not None else '')
        if self.left_child is not None:
            self.left_child.print()
            self.right_child.print()

    """ Obtain experiment name (bits=False) or bits (bits=True) of the bloom filter, then add to a graphviz Graph, then 
    iterate on children """
    def graphviz(self, graph, bits):
        graph.node(self.id,
                   'sim: ' + ''.join(map(str, map(int, self.sim_filter))) +
                   ('\nrem: ' + ''.join(map(str, map(int, self.rem_filter))) if self.rem_filter is not None else '')
                   if bits else self.experiment_name)
        if self.left_child is not None:
            graph.edge(tail_name=self.id, head_name=self.left_child.id)
            graph.edge(tail_name=self.id, head_name=self.right_child.id)
            self.left_child.graphviz(graph, bits)
            self.right_child.graphviz(graph, bits)
