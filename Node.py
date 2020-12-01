class Node(object):
    def __init__(self, bloom_filter_length, hash_functions, similarity_function, experiment_name=None,
                 bloom_filter=None):
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.similarity_function = similarity_function
        self.experiment_name = experiment_name
        self.bloom_filter = [0] * bloom_filter_length if bloom_filter is None else bloom_filter
        self.left_child = None
        self.right_child = None

    def copy(self):
        return Node(self.bloom_filter_length, self.hash_functions, self.similarity_function, self.experiment_name,
                    self.bloom_filter)

    def insert_kmer(self, kmer):
        for hash_function in self.hash_functions:
            self.bloom_filter[hash_function(kmer)] = 1

    def query_kmer(self, kmer):
        for hash_function in self.hash_functions:
            if self.bloom_filter[hash_function(kmer)] == 0:
                return False
        return True

    def union(self, node):
        for i in range(self.bloom_filter_length):
            self.bloom_filter[i] = self.bloom_filter[i] | node.bloom_filter[i]

    def similarity(self, node):
        return self.similarity_function(self.bloom_filter, node.bloom_filter)

    def insert_child_node(self, node):
        # 0 children - copy current node into left child and insert into right child
        if self.left_child is None and self.right_child is None:
            self.left_child = self.copy()
            self.experiment_name = None
            self.right_child = node
        # 1 child - insert into other child
        elif self.left_child is None:
            self.left_child = node
        elif self.right_child is None:
            self.right_child = node
        # 2 children - iterate into the more similar child
        elif self.left_child.similarity(node) >= self.right_child.similarity(node):
            self.left_child.insert_child_node(node)
        else:
            self.right_child.insert_child_node(node)
        # Union bloom filter
        self.union(node)

    def get_height(self):
        if self.left_child is None and self.right_child is None:
            return 0
        else:
            return max(self.left_child.get_height(), self.right_child.get_height())
